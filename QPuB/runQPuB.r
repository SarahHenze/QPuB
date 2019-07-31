####################################################################################################
# RUNQPUB.R FOR QPUB PACKAGE
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

suppressMessages(library(R.utils)) # commandArgs
suppressMessages(library(mcmcse))  # effective sample size (multiESS, minESS)

####################################################################################################
###### COMMAND LINE ARGUMENTS

# command line: Rscript runQPuB.r -infol inputfolder -outfol outputfolder -tim timepoints file -titr titration file

args = commandArgs(trailingOnly=FALSE, asValue=TRUE)
keys <- attachLocally(args)

if (!exists('infol')){
stop("This is the QPuB package. \n
At least one argument must be supplied: -infol.\n 
-infol inputfolder (COMPULSORY): folder containing all of the following input files \n AND folder 'data' with at least one csv file with input data. \n 
-outfol name of outputfolder (OPTIONAL): if not provided, default name will be 'OUTPUT_inputfolder(i)'. \n
-tim timepoints (OPTIONAL): txt with timepoints to compare, default: successive timepoints. \n
-titr titrationdata (OPTIONAL): csv with titration data, at least two values \n
For documentation see: https://github.com/QuantSysBio/QPuB/manual.\n", call.=FALSE)
}

# for use in R:
# setwd("/home/shenze/Downloads/QPuB-master/QPuB")
# file <- 'runQPuB.r'
# infol <- 'examples/Prot_K386'
# outfol <- 'my_fancy_name'
# titr <- '190423_K386_titration_substrate_charge_3_K386.csv'

####################################################################################################
###### AQUIRING PATHS

# get path of QPuB
if(dirname(file)=='.'){
      path_qpub <- getwd() 
} else { 
      path_qpub <- normalizePath(dirname(file))
}

# on Windows: format path
infol <- gsub(infol, pattern="\\", replacement="/", fixed=TRUE)

# get path of input folder
if(dirname(infol)=='.'){ # if only folder name provided

      ID <- infol
      path <- dirname(path_qpub)
      if(file.exists(file.path(paste(path,ID,sep='/'), fsep = .Platform$file.sep))){ # try directory of QPuB folder
            infol <- file.path(paste(path,ID,sep='/'), fsep = .Platform$file.sep)
      } else {
            stop(sprintf("Cannot find input folder '%s'. Typo? Otherwise, please provide path of the folder.",infol))
      }
} else if(grepl('/', infol)){ # if short path is provided

	ID <- unlist(strsplit(infol, "/"))[2]
	path <- dirname(path_qpub)
	shortpath <- unlist(strsplit(infol, c("/")))[1]
	path <- file.path(paste(path,shortpath, sep='/'), fsep = .Platform$file.sep)
      if(file.exists(file.path(paste(path,ID,sep='/'), fsep = .Platform$file.sep))){ # try directory of QPuB folder
            infol <- file.path(paste(path,ID,sep='/'), fsep = .Platform$file.sep)
      } else {
            stop(sprintf("Cannot find input folder '%s'. Typo? Otherwise, please provide path of the folder.",infol))
      }

} else { # if full path is provided
     path <- dirname(infol)
     ID <- basename(infol)
}

cat(paste0("\nPath of QPuB: ", path_qpub, "\n"))
cat(paste0("Path of input folder: ", infol, "\n"))

####################################################################################################
###### INITIALIZE THE SAMPLER

# read user input
source(file.path(paste(path_qpub,'inputparser.r',sep='/'), fsep = .Platform$file.sep))

# run data preparation
source(file.path(paste(path_qpub,'datapreparation.r',sep='/'), fsep = .Platform$file.sep))

# INFO: number of parameters numP+2 i.e. conversion factors of all products and sigma and punishment parameter

# prior distribution
para_min = c(rep(conv_lower,numP), sigma_lower, 0) # parameter minimum
para_max = c(rep(conv_upper,numP), sigma_upper, 1) # parameter maximum

# starting values
startvalue = c(c(rep(conv_start,numP)), sigma_start, pun_start) # parameter startvalues for conversion factors of all products, sigma, pun

# initialisation of the adaptation
mu <- c(rep(0,numP+2)) # mean
Cov <- diag(numP+2) # covariance matrix
gamma <- function(iter){ return(1/iter**GammaExponent) } # sequence to realize vanishing adaptation
scaling <- c(rep(2.38/sqrt(numP+2),numP+2)) # startvalues for global scaling

# convergence criterion based on ESS
if(TerminationCriterion==TRUE){
      cf = ConfLevel  # 1 - confidence level
      tol = Tolerance # Monte Carlo error in percent
      minSampleSize = minESS(p = numP+2, alpha = cf, eps = tol)
      cat(paste0("\nConvergence achieved when effective sample size > ", minSampleSize,'\n'))
}

# load additional functions
source(file.path(paste(path_qpub,'additionalfunctions.r',sep='/'), fsep = .Platform$file.sep))
source(file.path(paste(path_qpub,'MCMCroutine.r',sep='/'), fsep = .Platform$file.sep))

####################################################################################################
###### RUN MCMC

chain_backscaled <- run_metropolis_MCMC(startvalue=startvalue, minSample=minSampleSize, para_min=para_min, para_max=para_max, mu=mu, Cov=Cov)

####################################################################################################
##### OUTPUT OF FINAL RESULTS

cat("\nPRODUCING SUMMARY OF RESULTS...\n")
burnIn = round(dim(chain_backscaled)[1]*burnFrac)

##### SUMMARY OF PARAMETER DISTRIBUTION
summarize(chain_backscaled[-(1:burnIn),])

##### PROVIDE CONCENTRATIONS (with or without titration data)
concentrations(chain_backscaled[-(1:burnIn),])

####################################################################################################

if(exists('joke')){ cat("") } #TODO

cat("\nTHE END.")

} # end file
