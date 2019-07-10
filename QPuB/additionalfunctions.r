####################################################################################################
# ADDITIONALFUNCTIONS.R FOR QPUB PACKAGE
# PLOT DISTRIBUTION AND HISTORY OF THE PARAMETERS
# PLOT RESIDUALS: CURRENT MASS DEVIATION
# CALCULATE ABSOLUTE CONCENTRATIONS (OPTIONAL: BASED ON TITRATION DATA)
# PLOT KINETICS BASED ON ABSOLUTE CONCENTRATIONS
# PLOT IMPROVEMENT OF TOTAL MASS DEVIATION OVER COURSE OF ROUTINE
####################################################################################################

{

suppressMessages(library(grDevices)) # cairo_pdf
suppressMessages(library(matrixStats)) # standard deviation of columns of a matrix (colSds)
      
####################################################################################################
######## PLOTTING ANALYTICS

chain_FLAG = FALSE

roundUpNice <- function(x, nice=c(1,2,3,4,5,6,7,8,9,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
# get the highest and lowest line to set ylim
ylow = list()
yhigh = list()
for(dg in 1:numD){
      maxline = max(unlist(lapply(signalP[[dg]], function(x) x[which.max(x)])))
      minline = min(unlist(lapply(signalP[[dg]], function(x) x[which.min(x)])))
      if(abs(maxline)>abs(minline)){
            ylow[[dg]] = -roundUpNice(abs(1.5*maxline))
            yhigh[[dg]] = roundUpNice(abs(1.5*maxline))
      } else {
            ylow[[dg]] = -roundUpNice(abs(1.5*minline))
            yhigh[[dg]] = roundUpNice(abs(1.5*minline))
      }
}

plotChain <- function(chain){
      
      Ncurr = dim(chain)[1] # length of current chain
      if(is.null(Ncurr)){ # before the run starts: plot residuals using startvalues
            Ncurr=1
            post = matrix(chain[-c(length(chain),length(chain)-1)], nrow=1) # only convfac (not sigma and pun)
            if(plot_chain==TRUE){ chain_FLAG = TRUE }
            plot_chain = FALSE
            
      } else {
            burnIn = round(1*Ncurr/2)
            post = chain[-(1:burnIn),-c(dim(chain)[2],dim(chain)[2]-1)] # chain without burnIn, only convfac (not sigma and pun)
      }
      
      ####################
      ### PLOTTING CHAIN
      
      if(plot_chain==TRUE){
            # cairo_pdf(normalizePath(paste('plots_diagnostics/chain_',Ncurr,'.pdf',sep=""), winslash="\\", mustWork=TRUE), width=10, height=15, onefile=TRUE) # raster graphics
 		pdf(paste('plots_diagnostics/chain_',Ncurr,'.pdf',sep=""), width=10, height=15, onefile=TRUE) # raster graphics
            # pdf(paste(sprintf('plots_%s/chain_',ID),Ncurr,".pdf",sep=""), width=10, height=15) # vector graphics
            par(mfrow = c(5,2), oma=c(0,2,0,0))

            for(para in 1:dim(chain)[2]){
                  hist(chain[-(1:burnIn),para], main=paste('parameter',para), xlab="parameter values (using chain excl. BurnIn N/2)", breaks=70, col='gray', las=1)
                  if(QMEold==TRUE && (para!=dim(chain)[2]) && (para!=dim(chain)[2]-1)){
                        abline(v=old[para], col='orange')
                  }
                  
                  plot(chain[,para], type = "l", xlab="iteration incl. BurnIn", ylab="parameter value", main = paste("parameter",para), las=1)
            }

            dev.off()
      }
      
      ####################
      ### PLOTTING RESIDUALS

      if(plot_residuals==TRUE){
            
            # pdf(normalizePath(paste('plots_diagnostics/residuals_',Ncurr,".pdf",sep=""), winslash="\\", mustWork=TRUE), width=15, height=5)
		pdf(paste('plots_diagnostics/residuals_',Ncurr,".pdf",sep=""), width=15, height=5)
            par(mfrow = c(1,5), oma=c(0,5,3.5,0), mar=c(5,1.5,6,1)) # bottom, left, top, and right

            COLORS = rainbow(numR) # numR=2: red cyan, numR=3: red green blue, numR=4: red green cyan purple

            for (dg in 1:numD){ # over digestions

                  numPdg = dim(signalF[[dg]][[1]])[1] # number of products in that digestion

                  for(tp in 1:(numCT+1)){

                        if(tp%%5==1){

                              plot(c(1:numA),rep(signalP[[dg]][[1]][tp],numA),
                                    type="l", col=COLORS[1], lwd=2,
                                    yaxt="n", ylim=c(ylow[[dg]],yhigh[[dg]]), ylab="", xlab="aa position", cex.lab=2,
                                    axes=FALSE)
                              if(tp==1){title("time point 0", cex.main=2)
                              }else{title(paste("time point",timepoints[tp-1,2],' - ',timepoints[tp-1,1]), cex.main=2)}

                              axis(side=1, at=seq(1,numA,2), labels=seq(1,numA,2))
                              axis(side=2, at=seq(ylow[[dg]],yhigh[[dg]],length.out=5), las=1) #TODO good for every dataset?
                              mtext("signal deviation", side=2, line=4, at=maxline*0.5, cex=1.3)

                              mtext(sprintf("Digestion %s",dg), side=3, line=-1, outer=TRUE, cex=2.3, font=2)

                        } else {

                              plot(c(1:numA),rep(signalP[[dg]][[1]][tp],numA),
                                    type="l", col=COLORS[1], lwd=2,
                                    yaxt="n", ylim=c(ylow[[dg]],yhigh[[dg]]), ylab="", xlab="aa position", cex.lab=2,
                                    axes=FALSE,
                                    main=paste("time point",timepoints[tp-1,2],' - ',timepoints[tp-1,1]), cex.main=2)
                              axis(1, at=seq(1,numA,2), labels=seq(1,numA,2))

                        }

                        for (rp in 1:numR){

                              lines(c(1:numA),rep(signalP[[dg]][[rp]][tp],numA), type="l", col=COLORS[rp], lwd=2)

                              for(aa in 1:numA){
                                    
                                    if(!sum(ppm[[dg]][,aa])==0){ # only plot those which are present in some product

                                          amount = apply(sweep(sweep(post, MARGIN=2, signalF[[dg]][[rp]][,tp], `*`,check.margin=FALSE), MARGIN=2, ppm[[dg]][,aa], '*',check.margin=FALSE), 1,sum)
      
                                          if(dg==1){
                                                Quant <- quantile(amount, probs = c(0.05, 0.5, 0.95))
                                                points(aa, Quant['50%'], col=COLORS[rp], pch=20)
                                                segments(aa, Quant['5%'], aa, Quant['95%'], col=adjustcolor(COLORS[rp],alpha.f=0.5), lwd=0.01)
      
                                          } else {
                                                boxplot(amount/mean(post[,dg-1]), at=aa, add=TRUE, col=COLORS[rp], outline=FALSE, axes=FALSE) #TODO
                                          }
                                    }
                              } # end aa
                        } # end rp
                  } # end tp

                  if(numD>1){plot.new()}

            }

            dev.off()

            if(chain_FLAG==TRUE){ plot_chain=TRUE }
            
      } # end plot residuals
} # end plotChain

####################################################################################################
######## SUMMARY OF PARAMETER DISTRIBUTION

summarize <- function(chain){
      
      ####################
      ### STATISTICS OF THE PARAMETERS
      
      Stats <- cbind(t(boxplot(chain, plot=FALSE)$stats), colMeans(chain), colSds(chain))
      Stats <- cbind(matrix(c(as.vector(data[[1]][[1]][-1,1]),'sigma','pun')), Stats)
      colnames(Stats) <- c('parameter','min','lower_quartile','median','upper_quartile','max', 'mean', 'standard_deviation')
      write.csv(Stats, 'statistics.csv', quote=FALSE, row.names=FALSE)
      cat("\tStatistical summary in statistics.csv\n")
      
      ####################
      ### PLOT SUMMARY OF PARAMETERS
      
      chain_pep = chain[, -c(dim(chain)[2],dim(chain)[2]-1)]
      
      pdf("boxplot_chain.pdf", width=10, height=15) #TODO split plot multipage for many products
      # postscript("boxplot_chain.ps", width=10, height=20) #TODO achse wieder gestaucht
      
      boxplot(chain_pep, at=seq(1,numP,1), col="blue", 
              xlim = c(0,numP+1), ylab="", ylim=c(0, max(chain_pep)), axes=FALSE, 
              horizontal=TRUE, varwidth=TRUE, outline=FALSE)
      box()
      axis(1, xlim=c(0, max(chain_pep)),col="blue",col.axis="blue", las=1)
      mtext("conversion factor",side=1,line=2.5, col="blue")
      title("Distribution of conversion factors of all products", line=1.5, cex=1.2, font=2)
      
      axis(2,seq(1,numP))
      mtext("products",side=2,col="black",line=2.5) 
      
      abline(h=1:numP, col="gray", lty=3)
      
      if(QMEold==TRUE){ #TODO
            points(old, seq(1,numP), col='green')
      }
      
      dev.off()
      
      cat("\tDistributions of conversion factors in boxplot_chain.pdf\n")
}

####################################################################################################
######## PLOT RELATION OF CONVERSION FACTOR AND PEPTIDE LENGTH

plot_length <- function(chain){ #TODO test
      len = seq(1:numA)
      seq = data[[1]][[1]][-1,1]
      # COLORS = rainbow()
      png("relation.png")
      par(mar=c(5,7,5,3))
      plot(0,0, col='white', xlim=c(0,numA), ylim=c(0,max(chain)), xlab="peptide length", ylab="", las=1)
      for(pep in 1:numP){
            Quant <- quantile(chain[,pep], probs = c(0.05, 0.5, 0.95))
            points(nchar(as.character(seq[pep])), Quant['50%'], col='blue', pch=20)
            segments(nchar(as.character(seq[pep])), Quant['5%'], nchar(as.character(seq[pep])), Quant['95%'], col=adjustcolor('blue',alpha.f=0.5), lwd=1)
      }
      title(main="Relation of conversion factor to peptide length")
      title(ylab="conversion factor", line = 4.5)
      dev.off()
}

####################################################################################################
######## PLOT MAXIMUM LIKELIHOOD

plotmaxlike <- function(Pmax, sigma_plotmaxlike, punplotmaxlike){
      
      cat("\nPLOTTING MAXIMUM LIKELIHOOD...\n")
      
      p1 = seq(0,Pmax,Pmax/20)
      p2 = seq(0,Pmax,Pmax/20)
      p3 = seq(0,Pmax,Pmax/20)
      p4 = seq(0,Pmax,Pmax/20)
      
      M = array(NA,c(length(p1),length(p2),length(p3),length(p4)))
      
      for (i in 1:length(p1)){
            if((i/10-i%/%10)==0){
                  print(length(p1)-i-1)
            }
            for (j in 1:length(p2)){
                  for (k in 1:length(p3)){
                        for (l in 1:length(p4)){
                              M[i,j,k,l] = likelihood(param=c(p1[i],p2[j],p3[k],p4[l],sigma_plotmaxlike, pun_plotmaxlike))[[1]]
                        }
                  }
            }
      }
      
      png("maxlike.png")
      scp <- scatterplot3d(which(M[,,,1]>=max(round(M)-0.5),arr.ind=TRUE), type='l', lwd=10, color=rainbow(1),
                           angle=60, main='4D maximal likelihood', box=TRUE, grid=TRUE, y.margin.add=2,
                           lab=c(7,7,1), lab.z=c(7,1), asp=1, scale.y=0.6, mar=c(3.5,3.5,4,2),
                           xlim=c(0,20.2), ylim=c(0,20.2), zlim=c(0,20.2),
                           x.ticklabs=seq(0,Pmax,Pmax/5), y.ticklabs=seq(0,Pmax,Pmax/5), z.ticklabs=seq(0,Pmax,Pmax/5),
                           xlab='param 1', ylab='', zlab='param 3', las=1)
      for (l in 2:length(p4)){
            scp_p <- scp$points3d(which(M[,,,l]>=max(round(M)-0.5),arr.ind=TRUE), type='l', lwd=10, col=rainbow(length(p4))[l])
      }
      legend.gradient(pnts=cbind(x =c(7.5,6.5,6.5,7.5), y =c(5,5,4,4)),cols=rainbow(length(p4)), limits=c(0,1), title='param 4')
      text(x=6.5, y=0.7, 'param 2', srt=60)
      
      dev.off()
      
} # end plotmaxlike

####################################################################################################
######## PROVIDE CONCENTRATIONS

# calculates the absolute concentrations resulting from input data using output conversion factors

concentrations <- function(chain){
      
      cat("\nCALCULATING ABSOLUTE CONCENTRATIONS...\n")
      
      ####################
      ### FITTING TITRATION DATA
      
      if(titr_given==TRUE){

            titr_list = list()
            for(rp in 2:dim(titration)[2]){
                  titr_list[[rp-1]] <- cbind(titration[1],titration[rp])
                  colnames(titr_list[[rp-1]]) = c('amount', 'intensity')
            }
            
            titr_coeff = numeric()
            for(rp in 1:length(titr_list)){
                  if(fit_sigmoidal==TRUE){
                        fit = nls(intensity ~ SSlogis(amount, Asym, xmid, scal), data = titr_list[[rp]]) #Asym/(1+exp((xmid-input)/scal))
                  } else {
                        fit = lm(amount ~ intensity, titr_list[[rp]])
                  }
                  titr_coeff = rbind(titr_coeff,as.numeric(unlist(coefficients(fit))))
            }
            titr_intercept = colMeans(titr_coeff)[1]
            titr_slope = colMeans(titr_coeff)[2]
            
            #TODO sigmoidal fit: get coeff - get means, sd - get concentrations
      }
      
      ####################
      ### CALCULATING ABSOLUTE CONCENTRATIONS
      
      samp = sample(dim(chain)[1], 1000, replace = FALSE) #TODO ohne burnin
                  
      conc_means <- list()
      conc_sd <- list()
      for(dg in 1:numD){ #TODO testen mit numD>1!!!
            conc_means[[dg]] = list()
            conc_sd[[dg]] = list()
            
            for(rp in 1:numR){
                  conc_means[[dg]][[rp]] = matrix(,nrow=0,ncol=dim(data[[dg]][[rp]])[2])
                  conc_sd[[dg]][[rp]] = matrix(,nrow=0,ncol=dim(data[[dg]][[rp]])[2])
                  
                  for(pep in 1:numP){
                        conc = matrix(,nrow=0,ncol=dim(data[[dg]][[rp]][pep,-1])[2])
                        
                        for(sam in 1:length(samp)){
                              conc_tmp = tryCatch( as.numeric(data[[dg]][[rp]][pep+1,-1] * chain[samp[sam],pep] * 1/titr_slope - titr_intercept), # concentrations with titration data
                                               error = function(e) { as.numeric(data[[dg]][[rp]][pep+1,-1] * chain[samp[sam],pep]) }) # concentrations without titration data
                              conc =  rbind(conc, conc_tmp)
                        }

                        conc_means[[dg]][[rp]] = rbind(conc_means[[dg]][[rp]], c(as.vector(data[[dg]][[rp]][pep+1,1]), colMeans(conc)))
                        conc_sd[[dg]][[rp]] = rbind(conc_sd[[dg]][[rp]], c(as.vector(data[[dg]][[rp]][pep+1,1]), colSds(conc)))
                        
                  } # end pep
            } # end rp
      } # end dg
      
      for(dg in 1:numD){ #TODO testen mit numD>1!!!
            for(rp in 1:numR){
                  colnames(conc_means[[dg]][[rp]]) = c('sequence', seq(0,numT-1,1))
                  write.csv(conc_means[[dg]][[rp]], sprintf('conc_means_%s.csv',rp), row.names=FALSE, quote=FALSE)
                  
                  colnames(conc_sd[[dg]][[rp]]) = c('sequence', seq(0,numT-1,1))
                  write.csv(conc_sd[[dg]][[rp]], sprintf('conc_sd_%s.csv',rp), row.names=FALSE, quote=FALSE)
            }
      }
      cat("Means and standard deviations of absolute concentrations in conc_means.csv and conc_sd.csv\n")

      ####################
      ### PLOT KINETICS
      
      if (!file.exists(sprintf('plots_kinetics_%s',ID))){
            dir.create(sprintf('plots_kinetics_%s',ID))
      }
      
      COLORS = rainbow(numR)
      for(dg in 1:numD){ #TODO testen mit numD>1!!!
            
            ### SUBSTRATE
            # png(normalizePath(paste(sprintf('plots_kinetics_%s',ID), sprintf('peptide%s.png',pep), sep='/'), winslash="\\", mustWork=TRUE))
            png(paste(sprintf('plots_kinetics_%s',ID), sprintf('substrate.png'), sep='/'))
            par(mar=c(5,7,5,3)) # bottom, left, top, and right
            
            ymin = min(data[[dg]][[1]][1,-1])
            ymax = max(data[[dg]][[1]][1,-1])
            for(rp in 1:numR){
                  if(min(data[[dg]][[rp]][1,-1]) < ymin){ ymin = min(data[[dg]][[rp]][1,-1]) }
                  if(max(data[[dg]][[rp]][1,-1]) > ymax){ ymax = max(data[[dg]][[rp]][1,-1]) }
            }
            for(rp in 1:numR){
                  plot(seq(1,numT,1), data[[dg]][[rp]][1,-1], ylim=c(ymin,ymax), xlab='time', ylab='', main='substrate', type='l', lwd=2, col=COLORS[rp], las=1)
                  par(new=TRUE)
            }
            title(ylab = "concentration", cex.lab = 1, line = 4.5)
            dev.off()
            
            ### PRODUCTS
            for(pep in 1:numP){
                  
                  ymin = min(as.numeric(conc_means[[dg]][[1]][pep,-1]))
                  ymax = max(as.numeric(conc_means[[dg]][[1]][pep,-1]))
                  for(rp in 1:numR){
                        # X = conc_means[[dg]][[rp]][pep,length(conc_means[[dg]][[rp]][pep,])] + conc_sd[[dg]][[rp]][pep,length(conc_sd[[dg]][[rp]][pep,])] #TODO not last timepoint but highest
                        X = as.numeric(conc_means[[dg]][[rp]][pep,-1]) + as.numeric(conc_sd[[dg]][[rp]][pep,-1]) #TODO not last timepoint but highest

                        if(min(X) < ymin){ ymin = min(X) }
                        if(max(X) > ymax){ ymax = max(X) }
                  }
            
                  # png(normalizePath(paste(sprintf('plots_kinetics_%s',ID), sprintf('peptide%s.png',pep), sep='/'), winslash="\\", mustWork=TRUE))
                  png(paste(sprintf('plots_kinetics_%s',ID), sprintf('peptide%s.png',pep), sep='/'))
                  par(mar=c(5,7,5,3))
                  for(rp in 1:numR){
                        x = seq(0,numT-1)
                        y = as.numeric(conc_means[[dg]][[rp]][pep,-1])
                        sd = as.numeric(conc_sd[[dg]][[rp]][pep,-1])
                        plot(x, y, ylim=c(ymin,ymax), xlab='time', ylab='', main=paste('peptide ',pep), type='o', col=COLORS[rp], las=1)
                        polygon(c(x,rev(x)), c(y+sd,rev(y-sd)), col=COLORS[rp], density=50, border=NA)
                        par(new=TRUE)
                  }
                  title(ylab="concentration", line=4.5)
                  dev.off()
            } # end pep
      } # end dg
      cat("Plots of the kinetics in folder plots_kinetics\n")
      
} # end concentrations

####################################################################################################
######## MASS DEVIATION PLOTS

massdevplot <- function(massdev){
      
      cat("\nPLOTTING IMPROVEMENT OF TOTAL MASS DEVIATION...\n")
      
      png('massdeviation.png')
      # postscript('massdeviation.png')
      plot(massdev,type='l', ylim=c(0, max(massdev)),
           main="total mass deviation", xlab='iterations', ylab='mass deviation', las=1)
      
      if(QMEold==TRUE){
            massdev_QME = 0
            for(aa in a0:numA){ # over amino acids
                  for (rp in 1:numR){ # over replicates
                        amount = sum(old * signalF[[dg]][[rp]][,numT] * scalefac1/scalefac2 * ppm[[dg]][,aa]) #TODO richtig so?
                        massdev_QME = massdev_QME + abs(amount-signalP[[dg]][[rp]][numT])

                  }
            }
            points(length(massdev), massdev_QME, pch=19)
            text(length(massdev), massdev_QME, labels='QMEold', cex= 1, pos=2)
      }
      dev.off()
      
      cat("Plot of total mass deviation in massdeviation.png\n")
      if(QMEold==TRUE){ cat(paste0("Total mass deviation of QME: ", massdev_QME, "\n")) }

} # end massdevplot

} # end file
