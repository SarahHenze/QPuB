####################################################################################################
# MCMCROUTINE.R FOR QPUB PACKAGE
# MCMC SCHEME BASED ON RANDOM WALK METROPOLIS HASTINGS ALGORITHM WITH GENERALIZED ADAPTIVE SCALING
# DEFINITION OF PROPOSAL FUNCTION
# DEFINITION OF PRIOR DISTRIBUTION
# DEFINITION OF POSTERIOR DISTRIBUTION
# DEFINITION OF LIKELIHOOD FUNCTION
####################################################################################################
{

suppressMessages(library(tictoc)) # timer
suppressMessages(library(sys)) # timer
suppressMessages(library(tmvtnorm)) # proposal
suppressMessages(library(corpcor)) # semidef matrices
suppressMessages(library(dqrng)) # dqRNGkind("Xoroshiro128+")
suppressMessages(dqRNGkind("Xoroshiro128+")) # dqrunif
suppressMessages(library(coda))     # creation of mcmc object
suppressMessages(library(base))
suppressMessages(library(matrixcalc))   # is.postive.definite
suppressMessages(library(ggdmc))   #rtnorm
#suppressMessages(library(TruncatedNormal)) # distribution function
####################################################################################################
######## MCMC ALGORITHM

run_metropolis_MCMC <- function(startvalue, minSample, para_min, para_max, mu, Cov){

      cat("\nRUNNING INFERENCE...\n")
      tic("Total running time")
      starttime <- Sys.time()
      cat(paste0("Starting time: ",starttime, "\n"))

      # keep track of total mass deviation
      if(plot_massdev==TRUE){ massdev <- c() }

      # initialize chain of parameter values
      chain <- array(dim = c(1,length(startvalue)))
      chain[1,] = startvalue
      post_chain = posterior(chain[1,], para_min, para_max)

      # plot with starting values to see improvement
      plotChain(chain[1:1,])
      # stop('Absicht')

      # keep track of acceptance rate, is it optimal?
      if(Track==TRUE){
            accCount = 0 # counter of proposals accepted
            accHist = numeric(1) # history of acceptance
      }

      stopFlag = 0
      iter = 1
      ESS = array(c(rep(NA,length(startvalue))))
      # run algorithm until chain converged
      MaxIterEnable <- FALSE
      while (stopFlag!=1){
	      if((iter/1000-iter%/%1000)==0){
                  cat(paste("Iteration:",iter,"\n"))
	            currtime <- Sys.time()
	            cat(paste("Time elapsed: ",currtime - starttime,"\n"))
	      }

            proposal = proposalfunction(chain[iter,], Cov, para_min, para_max, scaling)
            post_prop = posterior(proposal, para_min, para_max)
            accProb = min(1, exp(post_prop[[1]]   - post_chain[[1]]))

            temp_chain = c(rep(NA,length(startvalue)))
            chain<-rbind(chain,temp_chain)


            if (dqrunif(1) <= accProb){
                  # acceptance
                  if(exists('accCount')){ accCount = accCount+1 }
                  if(exists('accHist')){ accHist <- cbind(accHist, c(1)) }
                  chain[iter+1,] = proposal
                  if(exists('massdev')){ massdev = rbind(massdev, post_prop[[2]]) }
                  post_chain = post_prop
            } else {
                  # rejection
                  chain[iter+1,] = chain[iter,]
                  if(exists('massdev')){ massdev = rbind(massdev, post_chain[[2]]) }
            }

            if(exists('accCount')){
                  accRate = accCount/iter
                  cat(sprintf("Acceptance rate: %s\n", accRate)) # TODO still warning!
            }

            # adaptation


            Adap = adaptation(mu, Cov, scaling, chain, iter, accProb)
            mu = Adap[[1]]
            Cov = Adap[[2]]
            scaling = Adap[[3]]

            # back scaling the chain
            #chain_backscaled = chain
            #chain_backscaled[,1:numP] = chain_backscaled[,1:numP] * scalefac1


            if (TerminationCriterion==TRUE & (!MaxIterEnable)){
                  # removing 50% of the samples as burn-In. One can remove 10% or 25% depending on which produces the highest ESS.
                  if (iter >=1000){
                  mcmcObj = mcmc(chain, start=floor(iter/2), end=iter)
                  effSampleSize <- multiESS(mcmcObj)


                  if(Track==TRUE){ cat(paste("Effective sample size: ",effSampleSize, "\n")) }
                  if (effSampleSize >= minSample){
                        cat("\n\nStop run: Chain converged.\n\n")
                        stopFlag = 1

                        # finally
                        plot_chain = TRUE
                        plot_residuals = TRUE
                        plotChain(chain[1:iter,])
                        save(massdev,file='massdev.RData')
                        chain_backscaled = chain* scalefac1
                        #chain_backscaled[,1:numP] = chain_backscaled[,1:numP] * scalefac1
                        save(chain_backscaled,file='chain_backscaled.RData')
                        cat(paste0("Finishing time: ",Sys.time(), "\n"))
                        break
                     }

               }

            }
            else
            {

            if (iter == MaxIter)
               {
                   stopFlag = 1
                   plot_chain = TRUE
                   plot_residuals = TRUE
                   plotChain(chain[1:iter,])
                   save(massdev,file='massdev.RData')
                   chain_backscaled = chain* scalefac1
                   #chain_backscaled[,1:numP] = chain_backscaled[,1:numP] * scalefac1
                   save(chain_backscaled,file='chain_backscaled.RData')
                   cat(paste0("Finishing time: ",Sys.time(), "\n"))
                   break
               }


            }

            # create plots and save chain
            if (plotdet==TRUE){
                  if(iter<1000){
                        if((iter/100-iter%/%100)==0){
                              plotChain(chain[1:iter,])
                        }
                  } else if (1000<=iter && iter<10000) {
                        if((iter/1000-iter%/%1000)==0){
                              plotChain(chain[1:iter,])
                        }
                  }
            }
      	if((iter/freq_plot-iter%/%freq_plot)==0){
                  plotChain(chain[1:iter,])
      	}
            if((iter/freq_save-iter%/%freq_save)==0){
                  save(chain,file=sprintf('chain_%s.RData',ID))
                  save(massdev,file=sprintf('massdev.RData',ID))
            }

            iter = iter + 1

      } # end while

      toc()

      if(!exists('massdev')){ massdev=NULL }
      list(chain_backscaled,massdev)


}

####################################################################################################
######## PROPOSAL FUNCTION

if(uniformprior==TRUE){

      proposalfunction <- function(param, Cov, para_min, para_max, scaling)

          {

            Cov <- diag(scaling) * Cov * diag(scaling)
            posFlag=0

            if (!is.diagonal.matrix(Cov))
            {
                stop("Covariance matrix in the proposal is not a diagonal matrix!")

            }
            if (det(Cov)<=0)
            {

                Cov = make.positive.definite(Cov,tol=1e-3)
                posFlag=1

            }
            if (posFlag == 0)
            {
                rtmvnorm(1, mean = param, sigma = Cov, lower=para_min, upper=para_max, algorithm = 'gibbs', burn.in.samples=100, thinning=20)
            }
           else
            {
                 MaxIterEnable <<- TRUE
                 prop_param = c(rep(NA, length(param)))
                 for (i in 1:length(param))
                   {
                      prop_param[i] = rtnorm(n = 1,param[i], Cov[i,i], lower=para_min[i], upper=para_max[i])

                   }
                   return(prop_param)
            }
      }
}

    else {

      proposalfunction <- function(param, Cov, para_min, para_max, scaling){

            Cov <- diag(scaling) * Cov * diag(scaling)
            print(Cov)
            if (det(Cov)<=0){ Cov = make.positive.definite(Cov,tol=1e-3) }
            rtmvnorm(1, mean = param, sigma = Cov, lower=0, upper=Inf, algorithm="gibbs", burn.in.samples=burnGibbs, thinning=thinGibbs) #TODO lower/upper vector?

      }
}






####################################################################################################
######## LIKELIHOOD FUNCTION

likelihood <- function(param){

      cf = param[1:numP]
      sigma = param[numP+1]
      pun = param[numP+2]

      if(plot_massdev==TRUE){massdev_curr = 0}

      likelihood = 0

      for(tp in 2:(numCT+1)){ # over list of timepoints to compare
            for(aa in a0:numA){ # over amino acids
                  for (dg in 1:numD){ # over digestions

                        numPdg = dim(signalF[[dg]][[1]])[1] # number of products in this digestion

                        for (rp in 1:numR){ # over replicates
                                    amount = 0

                                    for (j in 1: numP)
                                    {
                                         amount = amount + cf[j]*signalF[[dg]][[rp]][j,tp]*ppm[[dg]][j,aa]
                                    }

                                    if (dg==1) { # first digest must be the longest

                                          if(signalP[[dg]][[rp]][tp] < amount){

                                                likelihood = likelihood + log(pun) + dnorm(amount-signalP[[dg]][[rp]][tp], 0.0, pun*sigma,log=TRUE)

                                          } else {

                                                likelihood = likelihood + dnorm(amount-signalP[[dg]][[rp]][tp], 0.0, sigma,log=TRUE)
                                          }

                                          if(exists('massdev_curr')){ massdev_curr = massdev_curr + abs(amount-signalP[[dg]][[rp]][tp]) }

                                    } else { # all other digests: respective substrate does not have convfac=1 but its parameter

                                          if(signalP[[dg]][[rp]][tp] < amount){

                                                likelihood = likelihood + log(pun) + dnorm(amount-signalP[[dg]][[rp]][tp]*cf[dg-1], 0.0, pun*sigma,log=TRUE)

                                          } else {

                                                likelihood = likelihood + dnorm(amount-signalP[[dg]][[rp]][tp]*cf[dg-1], 0.0, sigma,log=TRUE)
                                          }

                                          if(exists('massdev_curr')){ massdev_curr = massdev_curr + abs(amount-signalP[[dg]][[rp]][tp]) }
                                    }
                        } # end rp
                  } # end dg
            } # end aa
      } # end tp
      if(!exists('massdev_curr')){massdev_curr=NULL}
      list(likelihood,massdev_curr)
}

####################################################################################################
######## PRIOR DISTRIBUTION

if(uniformprior==TRUE){
      prior <- function(param, para_min, para_max){
            sum(dunif(param, min=para_min, max=para_max, log=TRUE))
      }
} else {
      prior <- function(param, para_min, para_max){ # note that para_min and para_max are not needed (set to NULL beforehand)
            sum(dlnorm(param, meanlog=param1, sdlog=param2, log=TRUE))
      }
}

####################################################################################################
######## POSTERIOR DISTRIBUTION

posterior <- function(param, para_min, para_max){

      list(likelihood(param)[[1]] + prior(param, para_min, para_max), # return posterior
           likelihood(param)[[2]]) # return massdev
}

####################################################################################################
######## ADAPTIVE MCMC ALGORITHM

if (RaoBlack==TRUE){
      # Rao-Blackwellised Adaptive-MH algorithm
      adaptation <- function(mu, Cov, scaling, chain, iter, accProb){

            #if (iter < 100)
            #    Cov <- diag(359)
            #else
            Cov <- Cov + gamma(iter+1) * (accProb * outer(chain[iter+1,]-mu,chain[iter+1,]-mu) + (1-accProb) * outer(chain[iter,]-mu,chain[iter,]-mu) - Cov)
            mu <- mu + gamma(iter+1) * (accProb*(chain[iter+1,]-mu) + (1-accProb)*(chain[iter,]-mu))
            scaling <- scaling * exp(gamma(iter+1)*(accProb-OptAccRate)) # Global adaptive scaling

            return(list(mu, Cov, scaling))
      }
} else {
      # Adaptive-MH algorithm (default)
      adaptation <- function(mu, Cov, scaling, chain, iter, accProb){
            Cov <- Cov + gamma(iter+1) * (outer(chain[iter+1,]-mu,chain[iter+1,]-mu) - Cov)
            mu <- mu + gamma(iter+1) * (chain[iter+1,]-mu)
            scaling <- scaling * exp(gamma(iter+1)*(accProb-OptAccRate)) # Global adaptive scaling

            return(list(mu, Cov, scaling))
      }
}

} # end file
