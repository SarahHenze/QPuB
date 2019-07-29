####################################################################################################
# MCMCROUTINE.R FOR QPUB PACKAGE
# MCMC SCHEME BASED ON RANDOM WALK METROPOLIS HASTINGS ALGORITHM WITH GENERALIZED ADAPTIVE SCALING
# DEFINITION OF PROPOSAL FUNCTION
# DEFINITION OF PRIOR DISTRIBUTION
# DEFINITION OF POSTERIOR DISTRIBUTION
# DEFINITION OF LIKELIHOOD FUNCTION
####################################################################################################
#
# (C) Developers Sarah Henze, Debdas Paul, Juliane Liepe
# Research group Quantitative and Systems Biology
# Max Planck Institute for Biophysical Chemistry
# Goettingen, Germany
#
# This project is licensed under the GNU General Public License v3.0 
#
####################################################################################################

{ # start file

suppressMessages(library(tictoc)) # timer
suppressMessages(library(sys)) # timer
suppressMessages(library(tmvtnorm)) # proposal
suppressMessages(library(corpcor)) # semidef matrices
suppressMessages(library(dqrng)) # dqRNGkind("Xoroshiro128+")
suppressMessages(dqRNGkind("Xoroshiro128+")) # dqrunif
suppressMessages(library(coda)) # creation of mcmc object
suppressMessages(library(base)) # Base R functions
suppressMessages(library(matrixcalc)) # is.positive.definite
suppressMessages(library(ggdmc)) #rtnorm

####################################################################################################
######## MCMC ALGORITHM

run_metropolis_MCMC <- function(startvalue, minSample, para_min, para_max, mu, Cov){

      cat("\nRUNNING INFERENCE...\n")
	cat("Diagnostic plots of the chain and the residuals can be found in folder 'plots_diagnostics'.\n")
      
      # set a timer
      tic("Total running time")
      starttime <- Sys.time()
      cat(paste0("Starting time: ",starttime, "\n"))

      # initialize chain of parameter values
      chain <- array(dim = c(1,length(startvalue)))
      chain[1,] = startvalue
      post_chain = posterior(chain[1,], para_min, para_max)

      # plot with starting values before QPuB
      plotChain(chain[1:1,])

      # print acceptance rate and effective sample size to track performance of the run
      if(Track==TRUE){
            accCount = 0 # counter of proposals accepted
            accHist = numeric(1) # history of acceptance
      }

      # run algorithm until Markov chain converged
      iter = 1
      ESS = array(c(rep(NA,length(startvalue))))
      stopFlag = 0
      MaxIterEnable <- FALSE
      while (stopFlag!=1){
	      if((iter/1000-iter%/%1000)==0){
                  cat(paste("Iteration:",iter,"\n"))
	      }

            # proposal step
            proposal = proposalfunction(chain[iter,], Cov, para_min, para_max, scaling)
            post_prop = posterior(proposal, para_min, para_max)
            accProb = min(1, exp(post_prop[[1]] - post_chain[[1]]))

            # assign space in chain for new iteration 
            temp_chain = c(rep(NA,length(startvalue)))
            chain <- rbind(chain,temp_chain)

            if (dqrunif(1) <= accProb){
                  # acceptance
                  if(exists('accCount')){ accCount = accCount+1 }
                  if(exists('accHist')){ accHist <- cbind(accHist, c(1)) }
                  chain[iter+1,] = proposal
                  post_chain = post_prop
            } else {
                  # rejection
                  chain[iter+1,] = chain[iter,]
            }

            if(exists('accCount') && ((iter/1000-iter%/%1000)==0)){
                  accRate = accCount/iter
                  cat(sprintf("Acceptance rate: %s\n", accRate))
            }

            # adaptation
            Adap = adaptation(mu, Cov, scaling, chain, iter, accProb)
            mu = Adap[[1]]
            Cov = Adap[[2]]
            scaling = Adap[[3]]

            # check for convergence
            if(TerminationCriterion==TRUE && (!MaxIterEnable)){ # using convergence criterion based on effective sample size
                  if((iter >= ESSiter) && ((iter/1000-iter%/%1000)==0)){
                        mcmcObj = mcmc(chain, start=floor(iter/2), end=iter) # Removing 50% of the samples as burn-In. One can remove 10% or 25% depending on which produces the highest ESS.

                        # Sometimes there is the following error: 
                              # Error in if (effSampleSize >= minSample) { : missing value where TRUE/FALSE needed
                              # Calls: run_metropolis_MCMC
                              # In addition: Warning message:
                              # In mcse.multi(chain, ...) : Not enough samples. Estimate is not positive definite.
				mcmcObj = tryCatch({ mcmc(chain, start=floor(iter/2), end=iter) },
							warning = function(cond){ 
                                                      message("Here's the original warning message:")
									message(cond)
									stop("If it says 'Not enough samples. Estimate is not positive definite.' then try setting ESSiter to a higher value and restart.")})

                        effSampleSize <- multiESS(mcmcObj)
                        if(Track==TRUE){ cat(paste("Effective sample size: ",effSampleSize, "\n")) }
                        
                        if(effSampleSize >= minSample){
                              cat("\nStop run: Chain converged.\n\n")
                              cat(paste0("Finishing time: ", Sys.time(), "\n"))
                              stopFlag = 1
            
                              # final plots
                              plot_chain = TRUE
                              plot_residuals = TRUE
                              plotChain(chain[1:iter,])
                              
                              # save the chain
                              save(chain,file='chain.RData')
                              
                              # backscaling the chain (needed if products are scaled to fit order of magnitude of substrate, see datapreparation.r)
                              chain_backscaled <- matrix(,nrow=0,ncol=numP)
                              chain_backscaled <- rbind(chain_backscaled, chain[,1:numP] * scale_prod) # backscale conversion factors
                              chain_backscaled <- cbind(chain_backscaled, chain[,c(numP+1,numP+2)]) # sigma and pun don't need to be backscaled
                              save(chain_backscaled,file='chain_backscaled.RData')
                              
                              break
                        }
                  }
            } else { # using maximum number of iterations
                  if(iter == MaxIter){
                        cat("\nStop run: Maximum number of iterations reached.\n\n")
                        cat(paste0("Finishing time: ", Sys.time(), "\n"))
                        stopFlag = 1
                        
                        # final plots
                        plot_chain = TRUE
                        plot_residuals = TRUE
                        plotChain(chain[1:iter,])
                        
                        # save the chain
                        save(chain,file='chain.RData')
                        
                        # backscaling the chain (needed if products are scaled to fit order of magnitude of substrate, see datapreparation.r)
                        chain_backscaled <- matrix(,nrow=0,ncol=numP)
                        chain_backscaled <- rbind(chain_backscaled, chain[,1:numP] * scale_prod) # backscale conversion factors
                        chain_backscaled <- cbind(chain_backscaled, chain[,c(numP+1,numP+2)]) # sigma and pun don't need to be backscaled
                        save(chain_backscaled,file='chain_backscaled.RData')
                        
                        break
                  }
            }

            # create diagnostic plots and save chain
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
                  save(chain,file='chain.RData')
            }

            iter = iter + 1

      } # end while

      toc() # stop timer

      return(chain_backscaled) # return final chain
}

####################################################################################################
######## PROPOSAL FUNCTION

proposalfunction <- function(param, Cov, para_min, para_max, scaling){

      Cov <- diag(scaling) * Cov * diag(scaling)
      posFlag = 0

      # covariance matrix must be diagonal
      if (!is.diagonal.matrix(Cov)){ stop("Covariance matrix in the proposal is not a diagonal matrix!") }
      
      # covariance matrix must be positive definite
      if (det(Cov)<=0){
            Cov = make.positive.definite(Cov,tol=1e-3)
            posFlag = 1
      }
      
      # sample from truncated multivariate normal distribution
      if (posFlag == 0){ # if positive definite
            rtmvnorm(1, mean=param, sigma=Cov, lower=para_min, upper=para_max, algorithm='gibbs', burn.in.samples=burnGibbs, thinning=thinGibbs)
      } else { # if forced to be positive definite
            MaxIterEnable <<- TRUE
            prop_param = c(rep(NA, length(param)))
            for (i in 1:length(param)){
                  prop_param[i] = rtnorm(n = 1,param[i], Cov[i,i], lower=para_min[i], upper=para_max[i])
            }
            return(prop_param)
      }
}

####################################################################################################
######## LIKELIHOOD FUNCTION

likelihood <- function(param){

      cf = param[1:numP]
      sigma = param[numP+1]
      pun = param[numP+2]

      likelihood = 0

      for(tp in 2:(numCT+1)){ # over list of timepoints to compare
            for(aa in a0:numA){ # over amino acids
		      for (rp in 1:numR){ # over replicates
                        amount = 0

                        for (j in 1:numP){ # over products
                             amount = amount + cf[j] * signalF[[rp]][j,tp] * ppm[j,aa]
                        }

                        if(signalP[[rp]][tp] < amount){

                              likelihood = likelihood + log(pun) + dnorm(amount-signalP[[rp]][tp], 0.0, pun*sigma,log=TRUE)
                              # The log(pun) needs to be added in the case of mass gain, because the multiplicative change in the standard deviation by the punishment factor pun induces a flattening of the curve.
                              # To make it continuous in 0 (signalP[[rp]][tp] = amount), the curve needs to be shifted upwards by log(pun).

                        } else {

                              likelihood = likelihood + dnorm(amount-signalP[[rp]][tp], 0.0, sigma,log=TRUE)
                        }
                  } # end rp
            } # end aa
      } # end tp
      
      likelihood # return likelihood
}

####################################################################################################
######## PRIOR DISTRIBUTION

prior <- function(param, para_min, para_max){
      
      sum(dunif(param, min=para_min, max=para_max, log=TRUE))
      
}

####################################################################################################
######## POSTERIOR DISTRIBUTION

posterior <- function(param, para_min, para_max){

      likelihood(param) + prior(param, para_min, para_max) # return posterior
}

####################################################################################################
######## ADAPTIVE MCMC ALGORITHM

# Rao-Blackwellised Adaptive-MH algorithm
adaptation <- function(mu, Cov, scaling, chain, iter, accProb){

      Cov <- Cov + gamma(iter+1) * (accProb * outer(chain[iter+1,]-mu,chain[iter+1,]-mu) + (1-accProb) * outer(chain[iter,]-mu,chain[iter,]-mu) - Cov)
      mu <- mu + gamma(iter+1) * (accProb*(chain[iter+1,]-mu) + (1-accProb)*(chain[iter,]-mu))
      scaling <- scaling * exp(gamma(iter+1)*(accProb-OptAccRate)) # Global adaptive scaling

      return(list(mu, Cov, scaling))
}

} # end file
