####################################################################################################
# INPUTPARSER.R FOR QPUB PACKAGE
# CREATES OUTPUT DIRECTORY
# READS IN INPUT PARAMETERS AND SETTINGS FROM INPUT FILE
# READS IN DATA FROM DATA FOLDER
# OPTIONAL: READS IN TIMEPOINTS TO COMPARE FROM TIMEPOINTS FILE (DEFAULT: SUCCESSIVE TIMEPOINTS)
# OPTIONAL: READS IN TITRATION DATA FROM TITRATION FILE
####################################################################################################

{

####################################################################################################
###### CREATE OUTPUT FOLDER

cat("CREATING OUTPUT DIRECTORY...\n")

ID <- fol

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

ID<-split_path(ID)[1]



out_try = sprintf("%s/Output/OUTPUT_%s", path, ID)
if (!file.exists(out_try)){
      dir.create(out_try)
      setwd(out_try)
      outname = paste0('OUTPUT_',ID)
      dir.create('plots_diagnostics')
} else {
      outnumber = 1
      out_flag=0
      while (out_flag!=1) {
		out_try = sprintf("%s/Output/OUTPUT_%s(%s)", path, ID, outnumber)
		#out_try = sprintf("%s/OUTPUT_%s", path, ID, outnumber)
            if (!file.exists(out_try)){
                  dir.create(out_try)
                  setwd(out_try)
                  outname = sprintf("OUTPUT_%s(%s)", ID, outnumber)
                  #outname = sprintf("OUTPUT_%s", ID, outnumber)
                  dir.create('plots_diagnostics')
                  out_flag = 1
            }
            outnumber = outnumber + 1
      }
}
cat(paste0("OUTPUT FOLDER: ", outname, "\n"))

# redirect console output to a file
cat("Terminal output redirected to runQPuB_ROUT.txt.\n")
sink('runQPuB_ROUT.txt', split=FALSE)

####################################################################################################
###### INPUT PARAMETERS

cat("READING INPUT PARAMETERS...\n")

input <- read.table(paste(path,fol,pf,sep='/'), sep=':', header=FALSE, blank.lines.skip=TRUE)

##### MANDATORY INPUT

# mechanism of underlying enzyme
if(any(input[,1]=="enzyme")){
      enzyme <- gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="enzyme"),2]))
      if(nchar(enzyme)==0){ stop("Mechanism of underlying enzyme not provided!\n") }
} else {
      stop("Mechanism of underlying enzyme not provided!\n")
}

# prior to use: uniform or normal
if(any(input[,1]=="uniformprior")){
      uniformprior <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="uniformprior"),2])))
      if(is.na(uniformprior)){
            uniformprior <- TRUE
            cat("NOTE: Uniform prior used.\n")
      }
} else {
      uniformprior <- TRUE
      cat("NOTE: Uniform prior used.\n")
}
# parameter 1 for prior
if(any(input[,1]=="param_lower")){
      param_lower <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="param_lower"),2])))
      if(is.na(param_lower)){
            if(uniformprior==TRUE){
                  # uniform prior
                  param_lower <- 10**(-2)
                  cat(sprintf("WARNING: Uniform prior: Lower bound set to default: param_lower = %s\n", param_lower))
            } else {
                  # normal prior
                  stop("Normal prior: Mean not provided!\n")
            }
      }
} else {
      if(uniformprior==TRUE){
            # uniform prior
            param_lower <- 10**(-2)
            cat(sprintf("NOTE: Uniform prior: Lower bound set to default: param_lower = %s\n", param_lower))
      } else {
            # normal prior
            stop("Normal prior: Mean not provided!\n")
      }
}
# parameter 2 for prior
if(any(input[,1]=="param_upper")){
      param_upper <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="param_upper"),2])))
      if(is.na(param_upper)){
            if(uniformprior==TRUE){
                  # uniform prior
                  param_upper <- 10**2
                  cat(sprintf("NOTE: Uniform prior: Upper bound set to default: param_upper = %s\n", param_upper))
            } else {
                  # normal prior
                  stop("Normal prior: Standard deviation not provided!\n")
            }
      }
} else {
      if(uniformprior==TRUE){
            # uniform prior
            param_upper <- 10**2
            cat(sprintf("NOTE: Uniform prior: Upper bound set to default: param2 = %s\n", param_upper))
      } else {
            # normal prior
            stop("Normal prior: Standard deviation not provided!\n")
      }
}
# scale products according to substrate
if(any(input[,1]=="scaleprod")){
      scaleprod <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="scaleprod"),2])))
      if(is.na(scaleprod)){ scaleprod <- TRUE }
} else {
      scaleprod <- TRUE
}

# how often to save the chain
if(any(input[,1]=="freq_save")){
      freq_save <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="freq_save"),2])))
       if(is.na(freq_save)){ freq_save <- 10000 }
} else {
      freq_save <- 10000
}

##### PLOTTING INPUT

# plot chain: distribution and history
if(any(input[,1]=="plot_chain")){
      plot_chain <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="plot_chain"),2])))
      if(is.na(plot_chain)){ plot_chain <- TRUE }
} else {
      plot_chain <- TRUE
}
# plot residuals
if(any(input[,1]=="plot_residuals")){
      plot_residuals <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="plot_residuals"),2])))
      if(is.na(plot_residuals)){ plot_residuals <- TRUE }
} else {
      plot_residuals <- TRUE
}
# how often to generate plots
if(any(input[,1]=="freq_plot")){
      freq_plot <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="freq_plot"),2])))
      if(is.na(freq_plot)){ freq_plot <- 10000 }
} else {
      freq_plot <- 10000
}
# plot detailed: every 100 under 1000, every 1000 under 10000  #TODO change?
if(any(input[,1]=="plotdet")){
      plotdet <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="plotdet"),2])))
      if(is.na(plotdet)){ plotdet <- FALSE }
} else {
      plotdet <- FALSE
}

##### OPTIONAL FINAL OUTPUT

# plot relation of conversion factor to peptide length
if(any(input[,1]=="plot_relation")){
      plot_relation <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="plot_relation"),2])))
      if(is.na(plot_relation)){ plot_relation <- TRUE }
} else {
      plot_relation <- TRUE
}

# compare to results of QMEold
if(any(input[,1]=="QMEold")){
      QMEold <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="QMEold"),2])))
      if(is.na(QMEold)){ QMEold <- FALSE }
} else {
      QMEold <- FALSE
}
# plot total mass deviation
if(any(input[,1]=="plot_massdev")){
      plot_massdev <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="plot_massdev"),2])))
      if(is.na(plot_massdev)){ plot_massdev <- TRUE }
} else {
      plot_massdev <- TRUE
}
# fit titration data using sigmoidal function instead of linear
if(any(input[,1]=="titr_sigmoidal")){ #TODO only when itr data provided
      titr_sigmoidal <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="titr_sigmoidal"),2])))
      if(is.na(titr_sigmoidal)){
            titr_sigmoidal <- FALSE
            cat("NOTE: Titration data fitted linearly.\n")
      }
} else {
      titr_sigmoidal <- FALSE
      cat("NOTE: Titration data fitted linearly.\n")
}
# provide absolute concentrations
if(any(input[,1]=="provideConc")){
      provideConc <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="provideConc"),2])))
      if(is.na(provideConc)){ provideConc <- TRUE }
} else {
      provideConc <- TRUE
}
# plot maximum likelihood
if(any(input[,1]=="plot_maxlike")){
      plot_maxlike <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="plot_maxlike"),2])))
      if(plot_maxlike==TRUE){
            paramax_plotmaxlike <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="paramax_plotmaxlike"),2])))
            if(is.na(paramax_plotmaxlike)){
                  stop("Maximum likelihood: Upper bound of parameter range not provided!\n")
            }
            sigma_plotmaxlike <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="sigma_plotmaxlike"),2])))
            if(is.na(sigma_plotmaxlike)){
                  stop("Maximum likelihood: Sigma not provided!\n")
            }
      } else if(!is.na(plot_maxlike)){
            plot_maxlike <- FALSE
      }
} else {
      plot_maxlike <- FALSE
}

##### ADVANCED ALGORITHM PARAMETERS

# optimal acceptance rate
if(any(input[,1]=="OptAccRate")){
      OptAccRate <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="OptAccRate"),2])))
      if(is.na(OptAccRate)){
            OptAccRate <- 0.234
            cat(sprintf("NOTE: Optimal acceptance rate set to default: OptAccRate = %s\n",OptAccRate))
      }
} else {
      OptAccRate <- 0.234
      cat(sprintf("NOTE: Optimal acceptance rate set to default: OptAccRate = %s\n",OptAccRate))
}
# track acceptance rate
if(any(input[,1]=="Track")){
      Track <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="Track"),2])))
      if(is.na(Track)){ Track <- FALSE }
} else {
      Track <- FALSE
}
# choose adaptive MCMC algorithm to use
if(any(input[,1]=="RaoBlack")){
      RaoBlack <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="RaoBlack"),2])))
      if(is.na(RaoBlack)){ RaoBlack <- TRUE }
} else {
      RaoBlack <- TRUE
}

# choose the termination criterion
if(any(input[,1]=="TerminationCriterion")){
      TerminationCriterion <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="TerminationCriterion"),2])))
      if(is.na(TerminationCriterion)){ TerminationCriterion <- TRUE }
} else {
      TerminationCriterion <- TRUE
}

if(any(input[,1]=="MaxIter")){
      MaxIter <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="MaxIter"),2])))
      if(is.na(MaxIter)){

      MaxIter <- 100000
      cat(sprintf("NOTE: Maximum iteration for MCMC is set to default: MaxIter= %s\n",MaxIter))

      }
      } else {
           MaxIter <- 100000
           cat(sprintf("NOTE: Maximum iteration for MCMC is set to default: MaxIter= %s\n",MaxIter))
     }


# Exponent for the vanishing factor
if(any(input[,1]=="GammaExponent")){
      GammaExponent <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="GammaExponent"),2])))
      if(is.na(GammaExponent)){
            OptAccRate <- 0.5
            cat(sprintf("NOTE: GammaExponent is set to default: GammaExponent = %s\n",GammaExponent))
      }
} else {
      GammaExponent <- 0.5
      cat(sprintf("NOTE: GammaExponent is set to default: GammaExponent = %s\n",GammaExponent))
}


# confidence level for multESS routine
if(any(input[,1]=="ConfLevel")){
      ConfLevel <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="ConfLevel"),2])))
      if(is.na(ConfLevel)){
            ConfLevel <- 0.05
            cat(sprintf("NOTE: ConfLevel is set to default: ConfLevel = %s\n",ConfLevel))
      }
} else {
      ConfLevel <- 0.05
      cat(sprintf("NOTE: ConfLevel is set to default: ConfLevel = %s\n",ConfLevel))
}

# tolerance of multESS routine
if(any(input[,1]=="Tolerance")){
      Tolerance <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="Tolerance"),2])))
      if(is.na(Tolerance)){
            Tolerance <- 0.05
            cat(sprintf("NOTE: Tolerance is set to default: Tolerance = %s\n",Tolerance))
      }
} else {
      Tolerance <- 0.05
      cat(sprintf("NOTE: Tolerance is set to default: Tolerance = %s\n",Tolerance))
}

####################################################################################################
#### STARTING PARAMETERS

# starting values for the conversion factors
if(any(input[,1]=="ConvStartValue")){
      ConvStartValue <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="ConvStartValue"),2])))
      if(is.na(ConvStartValue)){
            ConvStartValue <- 1
            cat(sprintf("NOTE: ConvStartValue is set to default for all the peptides: ConvStartValue = %s\n",ConvStartValue))
      }
} else {
      ConvStartValue <- 1
      cat(sprintf("NOTE: ConvStartValue is set to default for all the peptides: ConvStartValue = %s\n",ConvStartValue))
}

# initial standard deviation
if(any(input[,1]=="Sigma")){
      Sigma <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="Sigma"),2])))
      if(is.na(Sigma)){
           Sigma <- 10
            cat(sprintf("NOTE: Sigma is set to default: Sigma = %s\n",Sigma))
      }
} else {
      Sigma <- 10
      cat(sprintf("NOTE: Sigma is set to default: Tolerance = %s\n",Sigma))
}

# upper and lower bounds for the standard deviation
if(any(input[,1]=="sigma_upper")){
      sigma_upper <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="sigma_upper"),2])))
      if(is.na(sigma_upper)){
           sigma_upper <- 10**(4)
            cat(sprintf("NOTE: Upper bound for Sigma is set to default: sigma_upper = %s\n",sigma_upper))
      }
} else {
      sigma_upper <- 10**(4)
      cat(sprintf("NOTE: Upper bound for Sigma is set to default: sigma_upper = %s\n",sigma_upper))
}

if(any(input[,1]=="sigma_lower")){
      sigma_lower=as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="sigma_lower"),2])))
      if(is.na(sigma_lower)){
           sigma_upper <- 10**(-4)
            cat(sprintf("NOTE: Lower bound for Sigma is set to default: sigma_lower = %s\n",sigma_lower))
      }
} else {
      sigma_lower <- 10**(-4)
      cat(sprintf("NOTE: Lower bound for Sigma is set to default: sigma_upper = %s\n",sigma_lower))
}



# Initial value of the punishment parameter
if(any(input[,1]=="Punishment")){
      Punishment <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="Punishment"),2])))
      if(is.na(Punishment)){
           Punishment <- 0.1
            cat(sprintf("NOTE: Punishment parameter is set to default: Punishment  = %s\n",Punishment))
      }
} else {
      Punishment <- 0.1
      cat(sprintf("NOTE: Punishment parameter is set to default: Punishment = %s\n",Punishment))
}

# burnIn fraction: 0.25 means 25% burnIn
if(any(input[,1]=="burnFrac")){
      burnFrac <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="burnFrac"),2])))
      if(is.na(burnFrac)){
           burnFrac <- 0.5
            cat(sprintf("NOTE: burnIn fraction is set to default: burnFrac  = %s\n",burnFrac))
      }
} else {
      burnFrac <- 0.5
      cat(sprintf("NOTE: burnIn fraction is set to default: burnFrac = %s\n",burnFrac))
}

# number of samples as burnIn for the Gibbs sampler in rtmvnorm
if(any(input[,1]=="burnGibbs")){
      burnGibbs <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="burnGibbs"),2])))
      if(is.na(burnGibbs)){
           burnFrac <- 1000
            cat(sprintf("NOTE: burnIn samples for the Gibbs sampler in rtmvnorm is set to default: burnGibbs  = %s\n",burnGibbs))
      }
} else {
      burnFrac <- 1000
      cat(sprintf("NOTE: burnIn samples for the Gibbs sampler in rtmvnorm is set to default: burnGibbs = %s\n",burnGibbs))
}


# thinning parameter for the Gibbs sampler in rtmvnorm
if(any(input[,1]=="thinGibbs")){
      thinGibbs <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="thinGibbs"),2])))
      if(is.na(thinGibbs)){
            thinGibbs  <- 20
            cat(sprintf("NOTE: Thinning is set to default: Thinning  = %s\n",thinGibbs))
      }
} else {
      thinGibbs <- 20
      cat(sprintf("NOTE: Thinning  is set to default: Thinning = %s\n",thinGibbs))
}


####################################################################################################
#### FILTER DATA

# filtering some of the products
if(any(input[,1]=="FILTER")){
      FILTER <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="FILTER"),2])))
      if(is.na(FILTER)){
            FILTER  <- FALSE
            cat(sprintf("NOTE: FILTER is set to default: FILTER  = %s\n", FILTER))
      }
} else {
      thinGibbs <- FALSE
      cat(sprintf("NOTE: FILTER  is set to default: FILTER = %s\n", FILTER))
}

if (FILTER==TRUE)
{

   if(any(input[,1]=="Threshold")){
      Threshold <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="Threshold"),2])))
      if(is.na(Threshold)){
            thinGibbs  <- 10**6
            cat(sprintf("NOTE: Threshold for product filtering is set to default: Threshold  = %s\n",Threshold))
      }
   } else {
      Threshold <- 10**6
      cat(sprintf("NOTE: Threshold for product filtering is set to default: Threshold = %s\n",Threshold))
    }

}

print(input)




####################################################################################################
#### INPUT TITRATION

if(exists('titr')){
      titration <- read.table(paste(path, fol, titr, sep='/'), sep=',', header=TRUE)
      titr_given = TRUE
} else {
      titr_given = FALSE
}

####################################################################################################
#### INPUT QMEold

if(QMEold==TRUE){
      old <- read.table(paste(path, fol, 'convfac', sep='/'), header=FALSE)[2][[1]]
}

####################################################################################################
#### INPUT DATA

filenames <- list.files(paste(path, fol, 'data', sep='/'), pattern="*.csv", full.names=TRUE)
data <- list()
data[[1]] <- lapply(filenames, function(x){ read.csv(file=x, header=F) })
# TODO how to define dg of the files? e.g erap --> subfolders per digestion

####################################################################################################
#### INPUT TIMEPOINTS

if(exists('tim')){
      cat("\nReading in timepoint file:\n")
      timepoints <- read.table(paste(path, fol, tim, sep='/'), header=FALSE)
} else {
      # set default for timepoints to use: successive timepoints
      cat("Setting timepoints to default: comparing successive timepoints.\n") #TODO timepoints aus data
	timepoints <- matrix(ncol=2, byrow=TRUE,
                     data=c(0,1,
                            1,2,
                            2,3,
                            3,4))
}
print(timepoints)

} # end file
