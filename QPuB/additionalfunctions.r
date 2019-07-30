####################################################################################################
# ADDITIONALFUNCTIONS.R FOR QPUB PACKAGE
# PLOT DISTRIBUTION AND TRACE OF THE PARAMETERS
# PLOT RESIDUALS: CURRENT MASS DEVIATION
# CALCULATE ABSOLUTE CONCENTRATIONS (OPTIONAL: BASED ON TITRATION DATA)
# PLOT KINETICS BASED ON ABSOLUTE CONCENTRATIONS
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
zoom <- 1.5
maxline = max(unlist(lapply(signalP, function(x) x[which.max(x)])))
minline = min(unlist(lapply(signalP, function(x) x[which.min(x)])))
if(abs(maxline)>abs(minline)){
      ylow = -roundUpNice(abs(zoom*maxline))
      yhigh = roundUpNice(abs(zoom*maxline))
} else {
      ylow = -roundUpNice(abs(zoom*minline))
      yhigh = roundUpNice(abs(zoom*minline))
}

plotChain <- function(chain){
      
      Ncurr = dim(chain)[1] # length of current chain
      if(is.null(Ncurr)){ # before the run starts: plot residuals using startvalues
            Ncurr = 0
            post = matrix(chain[-c(length(chain),length(chain)-1)], nrow=1) # only convfac (not sigma and pun)
            if(plot_chain==TRUE){ chain_FLAG = TRUE }
            plot_chain = FALSE
      } else if(Ncurr<10000){
            burnIn = round(Ncurr*burnFrac)
            post = chain[-(1:burnIn),-c(dim(chain)[2],dim(chain)[2]-1)] # chain without burnIn, only convfac (not sigma and pun)
      } else {
            burnIn = round(Ncurr*burnFrac)
            afterburnin = chain[-(1:burnIn),-c(dim(chain)[2],dim(chain)[2]-1)] # chain without burnIn, only convfac (not sigma and pun)
            samp = sample(Ncurr*burnFrac, 1000, replace = FALSE)
            post = afterburnin[samp,] # sample of chain without burnin, only convfac (not sigma and pun)
      }
      
      ####################
      ### PLOTTING CHAIN
      
      if(plot_chain==TRUE){
            pdf(file.path(paste('plots_diagnostics/chain_',Ncurr,".pdf",sep=""), fsep = .Platform$file.sep), width=8.3, height=11.7, compress=TRUE) # compressed: 852.7kB, uncompressed: 9.5MB (iteration 100000)
            par(mfrow = c(5,2), oma=c(0,2,0,0))

            # peptide products
            for(para in 1:(dim(chain)[2]-2)){
                  hist(chain[-(1:burnIn),para], main=paste('peptide',para), xlab="parameter values (using chain excl. BurnIn)", breaks=70, col='gray', las=1)
                  plot(chain[,para], type = "l", xlab="iteration incl. BurnIn", ylab="parameter value", main = paste("peptide",para), las=1)
            }
            # sigma
            hist(chain[-(1:burnIn),dim(chain)[2]-1], main='sigma', xlab="parameter values (using chain excl. BurnIn)", breaks=70, col='gray', las=1)
            plot(chain[,dim(chain)[2]-1], type = "l", xlab="iteration incl. BurnIn", ylab="parameter value", main = 'sigma', las=1)
            # punishment parameter
            hist(chain[-(1:burnIn),dim(chain)[2]], main='punishment parameter', xlab="parameter values (using chain excl. BurnIn)", breaks=70, col='gray', las=1)
            plot(chain[,dim(chain)[2]], type = "l", xlab="iteration incl. BurnIn", ylab="parameter value", main = 'punishment parameter', las=1)

            dev.off()
      }
      
      ####################
      ### PLOTTING RESIDUALS

      if(plot_residuals==TRUE){
            
            pdf(file.path(paste('plots_diagnostics/residuals_',Ncurr,".pdf",sep=""), fsep = .Platform$file.sep), width=15, height=5, compress=TRUE)
            par(mfrow = c(1,5), oma=c(0,5,3.5,0), mar=c(5,1.5,6,1)) # bottom, left, top, and right

            COLORS = rainbow(numR) # numR=2: red cyan, numR=3: red green blue, numR=4: red green cyan purple

            for(tp in 1:(numCT+1)){

                  if(tp%%5==1){ # put 5 plots in a row, first one with axes, others without

                        plot(c(1:numA),rep(signalP[[1]][tp],numA),
                              type="l", col=COLORS[1], lwd=2,
                              yaxt="n", ylim=c(ylow,yhigh), ylab="", xlab="aa position", cex.lab=2,
                              axes=FALSE)
                        if(tp==1){title("time point 0", cex.main=2)
                        } else {title(paste("time point",time_comp[tp-1,2],' - ',time_comp[tp-1,1]), cex.main=2)}

                        axis(side=1, at=seq(1,numA,2), labels=seq(1,numA,2))
                        axis(side=2, at=seq(ylow,yhigh,length.out=5), las=1)
                        mtext("signal deviation", side=2, line=4, at=maxline*0.5, cex=1.3)

                  } else {

                        plot(c(1:numA),rep(signalP[[1]][tp],numA),
                              type="l", col=COLORS[1], lwd=2,
                              yaxt="n", ylim=c(ylow,yhigh), ylab="", xlab="aa position", cex.lab=2,
                              axes=FALSE,
                              main=paste("time point",time_comp[tp-1,2],' - ',time_comp[tp-1,1]), cex.main=2)
                        axis(1, at=seq(1,numA,2), labels=seq(1,numA,2))

                  }

                  for (rp in 1:numR){ # over replicates

                        lines(c(1:numA),rep(signalP[[rp]][tp],numA), type="l", col=COLORS[rp], lwd=2)

                        for(aa in 1:numA){ # over amino acids
                              
                              if(!sum(ppm[,aa])==0){ # only plot those which are present in some product

                                    # amount is post * signalF * ppm, sum over all products
                                    amount = apply(sweep(sweep(post, MARGIN=2, signalF[[rp]][,tp], `*`,check.margin=FALSE), MARGIN=2, ppm[,aa], '*',check.margin=FALSE), MARGIN=1, sum)

                                      Quant <- quantile(amount, probs = c(0.05, 0.5, 0.95))
                                      points(aa, Quant['50%'], col=COLORS[rp], pch=20)
                                      segments(aa, Quant['5%'], aa, Quant['95%'], col=adjustcolor(COLORS[rp],alpha.f=0.5), lwd=0.01)
                              }
                        } # end aa
                  } # end rp
            } # end tp

            dev.off()

            if(chain_FLAG==TRUE){ plot_chain=TRUE }
            
      } # end plot residuals
} # end plotChain

####################################################################################################
######## SUMMARY OF PARAMETER DISTRIBUTION

summarize <- function(ch){
      
      ####################
      ### STATISTICS OF THE PARAMETERS
      
      Stats <- cbind(t(boxplot(ch, plot=FALSE)$stats), colMeans(ch), colSds(ch))
      Stats <- cbind(matrix(c(as.vector(dat[[1]][-1,1]),'sigma','pun')), Stats)
      colnames(Stats) <- c('parameter','min','lower_quartile','median','upper_quartile','max', 'mean', 'standard_deviation')
      write.csv(Stats, 'statistics.csv', quote=FALSE, row.names=FALSE)
      cat("\tStatistical summary in statistics.csv\n")
      
      ####################
      ### PLOT SUMMARY OF PARAMETERS
      
      chain_pep = ch[, -c(dim(ch)[2],dim(ch)[2]-1)]
      
      pdf("statistics.pdf", width=11.7, height=8.3)
      par(mar=c(5,7,5,3)) # bottom, left, top, and right
      
      if(dim(chain_pep)[2]<50){

            boxplot(chain_pep, at=seq(1,numP,1), 
                  xlim = c(0,numP+1), ylim=c(conv_lower, conv_upper), las=1, log="y",
                  col="blue", varwidth=TRUE, outline=FALSE)
            
            title(main = "Distribution of conversion factors of all products", line=2)     
            title(ylab = "conversion factor (log10-scale)", line=4.5)
            title(xlab = "peptide products", line=2.5) 
            
            abline(v=1:numP, col="gray", lty=3)
           
      } else {
            num_pages <- ceiling(numP/50)
            
            # all but last page
            for(page in 1:(num_pages-1)){ 
                  boxplot(chain_pep[,seq(50*(page-1)+1,50*page)], at=seq(50*(page-1)+1,50*page),
                        xlim = c(50*(page-1),50*page+1), ylim=c(conv_lower, conv_upper), las=1, log="y", xaxt="n",
                        col="blue", varwidth=TRUE, outline=FALSE)
                  
                  axis(1, at=seq(50*(page-1)+1,50*page), labels=seq(50*(page-1)+1,50*page))
                  title(main = "Distribution of conversion factors of all products", line=2)     
                  title(ylab = "conversion factor (log10-scale)", line=4.5)
                  title(xlab = "peptide products", line=2.5) 
                  
                  abline(v=seq(50*(page-1)+1,50*page), col="gray", lty=3)
            }
            # last page
            rest <- numP - 50*(num_pages-1)+1 -1
            boxplot(chain_pep[,seq(50*(num_pages-1)+1,numP)], at=seq(50*(num_pages-1)+1,numP),
                  xlim = c(50*(num_pages-1),numP+1), ylim=c(conv_lower, conv_upper), las=1, log="y", xaxt="n",
                  col="blue", varwidth=TRUE, outline=FALSE)

            axis(1, at=seq(50*(num_pages-1)+1,numP), labels=seq(50*(num_pages-1)+1,numP))
            title(main = "Distribution of conversion factors of all products", line=2)     
            title(ylab = "conversion factor (log10-scale)", line=4.5)
            title(xlab = "peptide products", line=2.5) 
            
            abline(v=seq(50*(num_pages-1)+1,numP), col="gray", lty=3)
      }
      
      dev.off()
      
      cat("\tDistributions of conversion factors in summary.pdf\n")
}


####################################################################################################
######## PROVIDE CONCENTRATIONS

# calculates the absolute concentrations resulting from input data using output conversion factors

concentrations <- function(ch){
      
      cat("\nCALCULATING ABSOLUTE CONCENTRATIONS...\n")
      
      if (!file.exists('concentrations')){
            dir.create('concentrations')
      }
      if (!file.exists('plots_concentrations')){
            dir.create('plots_concentrations')
      }
      
      ####################
      ### SUBSTRATE CONCENTRATION
      
      if(exists('titr')){
            cat(sprintf("... using substrate titration from file %s\n ", titr))
            # subs_titr defined in inputparser
            # fit titration before matching with kinetics
            subs_titr_means <- rowMeans(subs_titr[,-1])
            subs_titr_means <- cbind(subs_titr[1],subs_titr_means)
            colnames(subs_titr_means) = c('amount', 'intensity')
            
            COLORS = rainbow(numR)
            titr_coeff = numeric()
            plotfct('substratetitration_given', width=5, height=5)
			par(mar=c(5,7,5,3)) # bottom, left, top, and right
            ymax = max(subs_titr_means['intensity'], subs_titr_means['intensity'])
            fit = lm(intensity ~ amount, subs_titr_means)
            titr_coeff = rbind(titr_coeff,as.numeric(unlist(coefficients(fit))))
            cat("Linear fit of the substrate titration before normalization:\n")
            cat(sprintf("\tintercept: %s, slope: %s\n",titr_coeff[1],titr_coeff[2]))
            plot(as.numeric(unlist(subs_titr_means['amount'])), as.numeric(unlist(subs_titr_means['intensity'])), ylim=c(0,ymax), col=COLORS[rp], type='o', xlab='', ylab='', las=1)
            abline(fit)
			title(xlab = 'concentration', line=3)
			title(ylab = 'signal intensity', line=5)
			title(main = 'substrate titration as provided', line=2)
            dev.off()
            
            titr_intercept_before = colMeans(titr_coeff)[1]
            titr_slope_before = colMeans(titr_coeff)[2]
            signal_required <- titr_slope_before * loaded + titr_intercept_before
			cat(sprintf('Using this fit, the signal for amount loaded %s is %s\n',loaded,signal_required))
      
            # use signal_required to match titration data and signal kinetics, for every replicate
			subs_kin <- times
			for(rp in 1:numR){
				  subs_kin <- cbind(subs_kin, as.numeric(dat[[rp]][1,-1]))
			}
            normconst <- list()
            titration <- data.frame(subs_titr_means[,1])
            for(rp in 1:numR){
                  normconst[[rp]] <- signal_required/subs_kin[1,rp+1]
                  titration <- data.frame(cbind(titration, subs_titr_means[,2]/normconst[[rp]]))
            }
            
            # now fit titration after matching with kinetics
            titr_list = list()
            for(rp in 2:dim(titration)[2]){
                  titr_list[[rp-1]] <- cbind(titration[1],titration[rp])
                  colnames(titr_list[[rp-1]]) = c('amount', 'intensity')
            }
            
            fit_sigmoidal <- FALSE
            titr_coeff = numeric()
            plotfct('substratetitration_normalized', width=5, height=5)
			par(mar=c(5,7,5,3)) # bottom, left, top, and right
            ymax = max(titr_list[[1]]['intensity'], titr_list[[2]]['intensity'])
            for(rp in 1:length(titr_list)){
                  fit = lm(intensity ~ amount, titr_list[[rp]])
                  titr_coeff = rbind(titr_coeff,as.numeric(unlist(coefficients(fit))))
                  
                  plot(as.numeric(unlist(titr_list[[rp]]['amount'])), as.numeric(unlist(titr_list[[rp]]['intensity'])), ylim=c(0,ymax), col=COLORS[rp], type='o', xlab='', ylab='', las=1)
                  abline(fit)
                  par(new=TRUE)
            }
			title(xlab = 'concentration', line=3)
			title(ylab = 'signal intensity', line=5)
			title(main = 'normalized substrate titration', line=2)
            dev.off()
            titr_intercept = colMeans(titr_coeff)[1]
            titr_slope = colMeans(titr_coeff)[2]
            cat("Linear fit of the substrate titration after normalization, mean over replicates:\n")
			cat(sprintf("\tintercept: %s, slope: %s\n", titr_intercept,titr_slope))
            
            # amount of substrate loaded for every replicate, loaded given by userinput
            loaded_vec <- list()
            for(rp in 1:numR){
                  loaded_vec[[rp]] <- loaded
            }
      } else {
            cat("WARNING: No substrate titration provided. Only calculating normalized signals!\n")
            
            # "amount" of substrate loaded here is the initial intensity
            loaded_vec <- list()
            for(rp in 1:numR){
                  loaded_vec[[rp]] = dat[[rp]][1,2]
            }
      }
      
      # substrate concentration
      subs_conc <- list()
      for(rp in 1:numR){
            subs_conc[[rp]] = tryCatch( dat[[rp]][1,-1] * 1/titr_slope - titr_intercept/titr_slope, # concentrations with titration data
                                    error = function(e) { dat[[rp]][1,-1] }) # concentrations without titration data
            # save substrate concentration to csv
            write.csv(subs_conc[[rp]], file.path(paste('concentrations', sprintf('substrate_%s.csv',rp), sep='/'), fsep = .Platform$file.sep), row.names=FALSE, quote=FALSE)
      }
      
      # plot substrate degradation
      plotfct(file.path(paste('plots_concentrations','substrate', sep='/'), fsep = .Platform$file.sep), width=5, height=5)
      par(mar=c(5,7,5,3)) # bottom, left, top, and right
      ymax = max(unlist(lapply(subs_conc,max)))
      plot(times, subs_conc[[1]], type='o', col=COLORS[1], ylim=c(0,ymax), las=1, xlab='', ylab='')
      if(numR>1){
            for(rp in 2:numR){
                  par(new=TRUE)
                  plot(times, subs_conc[[rp]], type='o', col=COLORS[rp], ylim=c(0,ymax), las=1, xlab='', ylab='')
            }
      }
      title(xlab = 'time', line=3)
      title(ylab = 'concentration', line=5)
      title(main = 'substrate degradation', line=2)
      dev.off()
      
      # plot substrate degraded
      plotfct(file.path(paste('plots_concentrations','substrate_degraded', sep='/'), fsep = .Platform$file.sep), width=5, height=5)
      par(mar=c(5,7,5,3)) # bottom, left, top, and right
      ymax = max(unlist(lapply(subs_conc,max)))
      plot(times, abs(subs_conc[[1]]-loaded_vec[[1]]), type='o', col=COLORS[1], ylim=c(0,ymax), las=1, xlab='', ylab='')
      if(numR>1){
            for(rp in 2:numR){
                  par(new=TRUE)
                  plot(times, abs(subs_conc[[rp]]-loaded_vec[[rp]]), type='o', col=COLORS[rp], ylim=c(0,ymax), las=1, xlab='', ylab='')
            }
      }
      title(xlab = 'time', line=3)
      title(ylab = 'concentration', line=5)
      title(main = 'substrate degraded', line=2)
      dev.off()
      
      ####################
      ### PRODUCT CONCENTRATIONS
      
      samp = sample(dim(ch)[1], 1000, replace = FALSE) #TODO what if Niter < 1000 ?
      
      # calculate 0.05, 0.5, 0.95 quantiles based on parameter distributions and substrate titration
      conc_median = list()
      conc_five = list()
      conc_ninetyfive = list()
      
      for(rp in 1:numR){
            conc_median[[rp]] = matrix(,nrow=0,ncol=dim(dat[[rp]])[2])
            conc_five[[rp]] = matrix(,nrow=0,ncol=dim(dat[[rp]])[2])
            conc_ninetyfive[[rp]] = matrix(,nrow=0,ncol=dim(dat[[rp]])[2])
            
            for(pep in 1:numP){
                  conc = matrix(,nrow=0,ncol=dim(dat[[rp]][pep,-1])[2])
                  
                  for(sam in 1:length(samp)){
                        conc_tmp = tryCatch( as.numeric(dat[[rp]][pep+1,-1] * ch[samp[sam],pep] * 1/titr_slope - titr_intercept), # concentrations with titration data
                                          error = function(e) { as.numeric(dat[[rp]][pep+1,-1] * ch[samp[sam],pep]) }) # concentrations without titration data
                              
                        conc =  rbind(conc, conc_tmp)
                  }
                  Quant <- apply(conc, 2, function(x){ quantile(x, probs = c(0.05, 0.5, 0.95)) })
                  conc_median[[rp]] = rbind(conc_median[[rp]], c(as.vector(dat[[rp]][pep+1,1]), Quant['50%',]))
                  conc_five[[rp]] = rbind(conc_five[[rp]], c(as.vector(dat[[rp]][pep+1,1]), Quant['5%',]))
                  conc_ninetyfive[[rp]] = rbind(conc_ninetyfive[[rp]], c(as.vector(dat[[rp]][pep+1,1]), Quant['95%',]))
                  
            } # end pep
      } # end rp
      
      # move to zero
      for(rp in 1:numR){
            for(pep in 1:numP){
                  conc_median[[rp]][pep,-1] <- as.numeric(conc_median[[rp]][pep,-1]) - as.numeric(conc_median[[rp]][pep,-1][1])
                  conc_five[[rp]][pep,-1] <- as.numeric(conc_five[[rp]][pep,-1]) - as.numeric(conc_five[[rp]][pep,-1][1])
                  conc_ninetyfive[[rp]][pep,-1] <- as.numeric(conc_ninetyfive[[rp]][pep,-1]) - as.numeric(conc_ninetyfive[[rp]][pep,-1][1])
            }
      }
      
      # save products concentrations to csv
      for(rp in 1:numR){
            colnames(conc_median[[rp]]) = c('sequence', times)
            write.csv(conc_median[[rp]], file.path(paste('concentrations', sprintf('conc_median_%s.csv',rp), sep='/'), fsep = .Platform$file.sep), row.names=FALSE, quote=FALSE)
            
            colnames(conc_five[[rp]]) = c('sequence', times)
            write.csv(conc_five[[rp]], file.path(paste('concentrations', sprintf('conc_five_%s.csv',rp), sep='/'), fsep = .Platform$file.sep), row.names=FALSE, quote=FALSE)
            
            colnames(conc_ninetyfive[[rp]]) = c('sequence', times)
            write.csv(conc_ninetyfive[[rp]], file.path(paste('concentrations', sprintf('conc_ninetyfive_%s.csv',rp), sep='/'), fsep = .Platform$file.sep), row.names=FALSE, quote=FALSE)
      }
	cat("Substrate and product amounts:\n")
	cat("\tNumeric values in folder 'concentrations'.\n")
      
      # plot product concentration kinetics
      for(pep in 1:numP){
      
            plotfct(file.path(paste('plots_concentrations', sprintf('peptide%s',pep), sep='/'), fsep = .Platform$file.sep), width=5, height=5)
            par(mar=c(5,7,5,3))
            
            # get y limits
            ymin = min(as.numeric(conc_five[[1]][pep,-1]), na.rm=T)
            ymax = max(as.numeric(conc_ninetyfive[[1]][pep,-1]), na.rm=T)
            for(rp in 1:numR){
                  if(min(as.numeric(conc_five[[rp]][pep,-1]), na.rm=T) < ymin){ ymin = min(as.numeric(conc_five[[rp]][pep,-1])) }
                  if(max(as.numeric(conc_ninetyfive[[rp]][pep,-1]), na.rm=T) > ymax){ ymax = max(as.numeric(conc_ninetyfive[[rp]][pep,-1])) }
            }
            
            # plot median and 0.05-0.95 quantile shaded area
            for(rp in 1:numR){
                  y = as.numeric(conc_median[[rp]][pep,-1])
                  five = as.numeric(conc_five[[rp]][pep,-1])
                  ninetyfive = as.numeric(conc_ninetyfive[[rp]][pep,-1])
                  plot(times, y, ylim=c(as.numeric(ymin),as.numeric(ymax)), xlab='', ylab='', main='', type='o', pch=19, lwd=2, col=COLORS[rp], las=1)
                  polygon(c(times,rev(times)), c(ninetyfive,rev(five)), col=COLORS[rp], density=30, border=NULL)
                  par(new=TRUE)
            }
            title(xlab='time', line=2.5)
            title(ylab="concentration", line=5)
            title(main=paste('peptide ',pep), line=2)
            dev.off()
      } # end pep

	cat("\tKinetic plots in folder 'plots_concentrations'.\n")
      
} # end concentrations

} # end file
