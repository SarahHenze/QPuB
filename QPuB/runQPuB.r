####################################################################################################
# RUNQPUB.R FOR QPUB PACKAGE
####################################################################################################

{
path <- dirname(getwd())
cat("\nLOADING REQUIRED PACKAGES...\n\n")

suppressMessages(library(R.utils)) #commandArgs
suppressMessages(library(mcmcse))  # effective sample size (multiESS, minESS)

####################################################################################################
###### COMMAND LINE ARGUMENTS

# command line: Rscript runQPuB.r -fol inputfolder -pf inputfile -tim timepoints -titr titrationdata

args = commandArgs(trailingOnly=TRUE, asValue=TRUE) #TODO inputfile maybe not neccessary

if (length(args)<2) {
 	stop("At least two arguments must be supplied: inputfolder inputfile.\n -fol inputfolder: folder containing all of the following input files AND folder data with at least one csv file with input data.\n -pf inputfile: txt or csv with input parameters\n -tim timepoints: txt or csv with timepoints to compare (OPTIONAL, default: successive timepoints -titr titrationdata: txt or csv with titration data, at least two values (OPTIONAL)", call.=FALSE)
 }

keys <- attachLocally(args)

# for use in Rstudio:

#fol <- 'Kras'
#pf <- 'input.txt'
#tim <- 'timepoints.txt'
#titr <- ''


source(sprintf('%s/R/inputparser.r',path))
source(sprintf('%s/R/datapreparation.r',path))


####################################################################################################
####################################################################################################
###### INITIALIZE THE SAMPLER

# INFO: number of parameters numP+2 i.e. conversion factors of all products and sigma and punishment parameter

if(uniformprior==TRUE){
      para_min = c(rep(param_lower,numP), sigma_lower,0) # parameter minimum
      para_max = c(rep(param_upper,numP), sigma_upper,1) # parameter maximum
} else { para_min = NULL; para_max = NULL } # not needed for normal prior


para_start <-1 # parameter startvalue for conversion factors of products
startvalue = c(c(rep(para_start,numP)),Sigma,Punishment) # parameter startvalues of all products, sigma, pun


mu <- c(rep(0,numP+2)) # mean
Cov <- diag(numP+2) # covariance matrix

gamma <- function(iter){ return(1/iter**GammaExponent) } # sequence to realize vanishing adaptation
scaling <- c(rep(2.38/sqrt(numP+2),numP+2)) # startvalues for global scaling

source(sprintf('%s/R/additionalfunctions.r',path))
source(sprintf('%s/R/MCMCroutine.r',path))


#source(normalizePath(sprintf('%s/QPuB/MCMCroutine.r',path), winslash="\\", mustWork=TRUE))
#source(normalizePath(sprintf('%s/QPuB/additionalfunctions.r',path), winslash="\\", mustWork=TRUE))


####################################################################################################
###### RUN MCMC

# calculate the minimum effective sample size according to
# Vats, D., Flegal, J. M., and, Jones, G. L Multivariate Output Analysis for Markov chain Monte Carlo, arXiv preprint arXiv:1512.07713 (2015).
# p = number of parameters
# alpha = confidence level
# eps   = tolerance level.

# selbst implementieren
# https://www.johndcook.com/blog/2017/06/27/effective-sample-size-for-mcmc/

cf = ConfLevel  # 1 - confidence level
tol  = Tolerance # Monte Carlo error in percent

minSampleSize = minESS(p = numP+2, alpha = cf, eps = tol)
#browser()
cat(paste0("Convergence achieved when effective sample size > ", minSampleSize,'\n'))

output = run_metropolis_MCMC(startvalue=startvalue, minSample=minSampleSize, para_min=para_min, para_max=para_max, mu=mu, Cov=Cov)

chain_backscaled <- output[[1]]
massdev <- output[[2]]

####################################################################################################
####################################################################################################
##### OUTPUT OF FINAL RESULTS
#
cat("\nPRODUCING SUMMARY OF RESULTS...\n")
burnIn = round(dim(chain_backscaled)[1]*burnFrac)
# # # #
# # # # ##### SUMMARY OF PARAMETER DISTRIBUTION
summarize(chain_backscaled[-(1:burnIn),])

# # # # ##### PLOT RELATION OF CONVERSION FACTOR TO PEPTIDE LENGTH
if(plot_relation==TRUE){ plot_length(chain_backscaled[-(1:burnIn),]) }

# # # # ##### PLOT IMPROVEMENT OF MASS DEVIATION
if(plot_massdev==TRUE){ massdevplot(massdev) }

# # # # ##### PROVIDE CONCENTRATIONS (with or without titration data)
if(provideConc==TRUE){ concentrations(chain_backscaled[-(1:burnIn),]) }
# # # #
# # # # ##### PLOT MAXIMUM LIKELIHOOD
if(plot_maxlike==TRUE){
     suppressMessages(library(scatterplot3d))
     suppressMessages(library(SDMTools)) # legend.gradient
     pun_plotmaxlike = 0.1 #TODO
     plotmaxlike(paramax_plotmaxlike, sigma_plotmaxlike, pun_plotmaxlike)
}

} # end file

