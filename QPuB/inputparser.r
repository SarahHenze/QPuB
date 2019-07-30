####################################################################################################
# INPUTPARSER.R FOR QPUB PACKAGE
# CREATE OUTPUT DIRECTORY
# READ IN INPUT PARAMETERS AND SETTINGS FROM INPUT FILE input.txt
# READ IN DATA FROM DATA FOLDER
# OPTIONAL: READ IN TIMEPOINTS TO COMPARE FROM TIMEPOINTS FILE (DEFAULT: SUCCESSIVE TIMEPOINTS)
# OPTIONAL: READ IN TITRATION DATA FROM TITRATION FILE
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

####################################################################################################
###### CREATE OUTPUT FOLDER

cat("\nCREATING OUTPUT DIRECTORY...\n")

if(!file.exists(file.path(sprintf("%s/QPuB-Output", path), fsep = .Platform$file.sep))){
      dir.create(file.path(sprintf("%s/QPuB-Output", path), fsep = .Platform$file.sep))
}

if(exists('outfol')){ # userdefined name
      out_try = file.path(sprintf("%s/QPuB-Output/%s", path, outfol), fsep = .Platform$file.sep)
      if(!file.exists(out_try)){
            dir.create(out_try)
            setwd(out_try)
            outname = outfol
            dir.create('plots_diagnostics')
      } else {
		do.call(file.remove, list(list.files(file.path(sprintf("%s/QPuB-Output/%s", path, outfol), fsep = .Platform$file.sep), full.names = TRUE))) # empty folder
            dir.create('plots_diagnostics')
      }
} else { # no userdefined name
      out_try = file.path(sprintf("%s/QPuB-Output/OUT_%s", path, ID), fsep = .Platform$file.sep)
      if (!file.exists(out_try)){
            dir.create(out_try)
            setwd(out_try)
            outname = paste0('OUT_',ID)
            dir.create('plots_diagnostics')
      } else {
            outnumber = 1
            out_flag=0
            while (out_flag!=1) {
      		out_try = file.path(sprintf("%s/QPuB-Output/OUT_%s(%s)", path, ID, outnumber), fsep = .Platform$file.sep)
                  if (!file.exists(out_try)){
                        dir.create(out_try)
                        setwd(out_try)
                        outname = sprintf("OUT_%s(%s)", ID, outnumber)
                        dir.create('plots_diagnostics')
                        out_flag = 1
                  }
                  outnumber = outnumber + 1
            }
      }
}
cat(paste0("Path of output folder: ", file.path(paste(path,"QPuB-Output", outname, sep='/'), fsep = .Platform$file.sep), "\nDo not change folder name during run! \n"))

# redirect console output to a file
cat("\nTerminal output redirected to runQPuB_ROUT.txt.\n")
cat("RUNNING...")
sink('runQPuB_ROUT.txt', split=FALSE)

####################################################################################################
###### INPUT PARAMETERS

cat("READING INPUT PARAMETERS...\n")

if(file.exists(file.path(paste(infol,'input.txt',sep='/'), fsep = .Platform$file.sep))){
      input <- read.table(file.path(paste(infol,'input.txt',sep='/'), fsep = .Platform$file.sep), sep=':', header=FALSE, blank.lines.skip=TRUE)
} else {
      stop("Cannot find input file. Note that the file has to be named 'input.txt'. Abort.")
}

##################################################
# LOW LEVEL INPUT
##################################################

########## MANDATORY INPUT

# mechanism of underlying enzyme
if(any(input[,1]=="enzyme")){
      enzyme <- gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="enzyme"),2]))
      if(enzyme!='endopep' && enzyme!='exopep'){ stop("Mechanism of underlying enzyme not provided! Maybe typo? Abort.\n") }
} else {
      stop("Mechanism of underlying enzyme not provided! Abort.\n")
}

# amount of substrate loaded for kinetics
if(exists('titr')){
      if(any(input[,1]=="loaded")){
            loaded <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="loaded"),2])))
            if(is.na(loaded)){
                  stop("Amount of substrate loaded for the kinetics not provided! Abort.\n")
            }
      } else {
		stop("Amount of substrate loaded for the kinetics not provided! Abort.\n")
	}
}

########## OPTIONAL INPUT WITH DEFAULTS

### DIAGNOSTIC PLOTS
# plot chain: distribution and history
if(any(input[,1]=="plot_chain")){
      plot_chain <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="plot_chain"),2])))
      if(is.na(plot_chain)){
            plot_chain <- TRUE             
            cat("WARNING: Typo in inputfile? Set plot_chain to TRUE by default.\n")
      }
} else {
      plot_chain <- TRUE
}
# plot residuals
if(any(input[,1]=="plot_residuals")){
      plot_residuals <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="plot_residuals"),2])))
      if(is.na(plot_residuals)){
            plot_residuals <- TRUE
            cat("WARNING: Typo in inputfile? Set plot_residuals to TRUE by default.\n")
      }
} else {
      plot_residuals <- TRUE
}
# how often to generate plots
if(any(input[,1]=="freq_plot")){
      freq_plot <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="freq_plot"),2])))
      if(is.na(freq_plot)){
            freq_plot <- 10000
            cat("WARNING: Typo in inputfile? Set freq_plot to 10000 by default.\n")
      }
} else {
      freq_plot <- 10000
}
# plot detailed: every 100 under 1000, every 1000 under 10000
if(any(input[,1]=="plotdet")){
      plotdet <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="plotdet"),2])))
      if(is.na(plotdet)){
            plotdet <- FALSE 
            cat("WARNING: Typo in inputfile? Set plotdet to FALSE by default.\n")
      }
} else {
      plotdet <- FALSE
}
# generate vector graphics instead of flat png
if(any(input[,1]=="vecplot")){
      vecplot <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="vecplot"),2])))
      if(is.na(vecplot)){
            vecplot <- FALSE
            cat("WARNING: Typo in inputfile? Set vecplot to TRUE by default.\n")
      }
} else {
      vecplot <- FALSE
}
if(vecplot==TRUE){
      plotfct <- function(x, width, height){ postscript(paste0(x,'.eps'), width=width, height=height, paper='special', horizontal=F) }
} else {
      plotfct <- function(x, width, height){ png(paste0(x,'.png'), width=width, height=height, units = "in", res=72) }
}

### BACKUP
# how often to save the chain
if(any(input[,1]=="freq_save")){
      freq_save <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="freq_save"),2])))
      if(is.na(freq_save)){
            freq_save <- 10000
            cat("WARNING: Typo in inputfile? Set freq_save to 10000 by default.\n")
      }
} else {
      freq_save <- 10000
}

##################################################
# ADVANCED INPUT
##################################################

### DATA MANIPULATION
# scale products to match order of magnitude of substrate
if(any(input[,1]=="scaleprod")){
      scaleprod <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="scaleprod"),2])))
      if(is.na(scaleprod)){ 
            scaleprod <- TRUE
            cat("WARNING: Typo in inputfile? Set scaleprod to TRUE by default. Signal intensities of products will be scaled to match order of magnitude of substrate!\n")
      }
} else {
      scaleprod <- TRUE
      cat("WARNING: Signal intensities of products will be scaled to match order of magnitude of substrate!\n")
}
# delete all products which have a signal below the threshold and do not contribute much
if(any(input[,1]=="filterprod")){
      filterprod <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="filterprod"),2])))
      if(is.na(filterprod)){
            filterprod  <- FALSE
            cat(sprintf("WARNING: Typo in inputfile? Products with low intensities are not filtered out by default: filterprod  = %s\n", filterprod))
      }
} else {
      filterprod <- FALSE
      cat(sprintf("WARNING: Products with low intensities are not filtered out by default: filterprod = %s\n", filterprod))
}
if(filterprod==TRUE){
      if(any(input[,1]=="threshold")){
            threshold <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="threshold"),2])))
            if(is.na(threshold)){
                  threshold  <- 10**6
                  cat(sprintf("NOTE: Typo in inputfile? Threshold for product filtering is set to default: threshold  = %s\n",threshold))
            }
      } else {
            threshold <- 10**6
            cat(sprintf("NOTE: Threshold for product filtering is set to default: threshold = %s\n",threshold))
      }
}

### UNIFORM PRIOR DISTRIBUTION
# lower bound of uniform prior range of conversion factors for all products
if(any(input[,1]=="conv_lower")){
      conv_lower <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="conv_lower"),2])))
      if(is.na(conv_lower)){
            conv_lower <- 10**(-4)
            cat(sprintf("NOTE: Typo in inputfile? Uniform prior: Lower bound for conversion factors set to default: conv_lower = %s\n", conv_lower))
      }
} else {
      conv_lower <- 10**(-4)
      cat(sprintf("NOTE: Uniform prior: Lower bound for conversion factors set to default: conv_lower = %s\n", conv_lower))
}
# upper bound of uniform prior range of conversion factors for all products
if(any(input[,1]=="conv_upper")){
      conv_upper <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="conv_upper"),2])))
      if(is.na(conv_upper)){
            conv_upper <- 10**4
            cat(sprintf("NOTE: Typo in inputfile? Uniform prior: Upper bound for conversion factors set to default: conv_upper = %s\n", conv_upper))
      }
} else {
      param_upper <- 10**4
      cat(sprintf("NOTE: Uniform prior: Upper bound for conversion factors set to default: param2 = %s\n", param_upper))
}
# lower bound of uniform prior range of the standard deviation in the data
if(any(input[,1]=="sigma_lower")){
      sigma_lower=as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="sigma_lower"),2])))
      if(is.na(sigma_lower)){
            sigma_lower <- 10**(-4)
            cat(sprintf("NOTE: Typo in inputfile? Uniform prior: Lower bound for sigma set to default: sigma_lower = %s\n",sigma_lower))
      }
} else {
      sigma_lower <- 10**(-4)
      cat(sprintf("NOTE: Uniform prior: Lower bound for sigma set to default: sigma_lower = %s\n",sigma_lower))
}
# upper bound of uniform prior range of the standard deviation in the data
if(any(input[,1]=="sigma_upper")){
      sigma_upper <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="sigma_upper"),2])))
      if(is.na(sigma_upper)){
            sigma_upper <- 10**(4)
            cat(sprintf("NOTE: Typo in inputfile? Uniform prior: Upper bound for sigma set to default: sigma_upper = %s\n",sigma_upper))
      }
} else {
      sigma_upper <- 10**(4)
      cat(sprintf("NOTE: Uniform prior: Upper bound for sigma set to default: sigma_upper = %s\n",sigma_upper))
}

### STARTING VALUES OF THE PARAMETERS
# initial values for the conversion factors for all products
if(any(input[,1]=="conv_start")){
      conv_start <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="conv_start"),2])))
      if(is.na(conv_start)){
            conv_start <- 1
            cat(sprintf("NOTE: Typo in inputfile? Starting value for conversion factors for all products is set to default: conv_start = %s\n",conv_start))
      }
	if((conv_start<conv_lower) || (conv_start>conv_upper)){ stop("Conversion factors: Starting value not in prior range!") }
} else {
      conv_start <- 1
      cat(sprintf("NOTE: Starting value for conversion factors for all products is set to default for all the peptides: conv_start = %s\n",conv_start))
}
# initial value for the standard deviation sigma
if(any(input[,1]=="sigma_start")){
      sigma_start <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="sigma_start"),2])))
      if(is.na(sigma_start)){
            sigma_start <- 10
            cat(sprintf("NOTE: Typo in inputfile? Starting value for sigma is set to default: sigma_start = %s\n",sigma_start))
      }
	if((sigma_start<sigma_lower) || (sigma_start>sigma_upper)){ stop("Sigma: Starting value not in prior range!") }
} else {
      sigma_start <- 10
      cat(sprintf("NOTE: Starting value for sigma is set to default: sigma_start = %s\n",sigma_start))
}
# initial value of the punishment parameter
if(any(input[,1]=="pun_start")){
      pun_start <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="pun_start"),2])))
      if(is.na(pun_start)){
            pun_start <- 0.1
            cat(sprintf("NOTE: Typo in inputfile? Starting value for the punishment parameter is set to default: pun_start  = %s\n",pun_start))
      }
} else {
      pun_start <- 0.1
      cat(sprintf("NOTE: Starting value for the punishment parameter is set to default: pun_start = %s\n",pun_start))
}

##### ADVANCED ALGORITHM PARAMETERS
# optimal acceptance rate
if(any(input[,1]=="OptAccRate")){
      OptAccRate <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="OptAccRate"),2])))
      if(is.na(OptAccRate)){
            OptAccRate <- 0.234
            cat(sprintf("NOTE: Typo in inputfile? Optimal acceptance rate set to default: OptAccRate = %s\n",OptAccRate))
      }
} else {
      OptAccRate <- 0.234
      cat(sprintf("NOTE: Optimal acceptance rate set to default: OptAccRate = %s\n",OptAccRate))
}
# track acceptance rate
if(any(input[,1]=="Track")){
      Track <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="Track"),2])))
      if(is.na(Track)){
            Track <- FALSE 
            cat(sprintf("NOTE: Typo in inputfile? Setting tracking of acceptance rate to default: Track = %s\n",Track))
      }
} else {
      Track <- FALSE
      cat(sprintf("NOTE: Setting tracking of acceptance rate to default: Track = %s\n",Track))
}
# Exponent for the sequence to realize vanishing adaptation
if(any(input[,1]=="GammaExponent")){
      GammaExponent <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="GammaExponent"),2])))
      if(is.na(GammaExponent)){
            GammaExponent <- 0.5
            cat(sprintf("NOTE: Typo in inputfile? Exponent of Gamma function is set to default: GammaExponent = %s\n",GammaExponent))
      }
} else {
      GammaExponent <- 0.5
      cat(sprintf("NOTE: Exponent of Gamma function is set to default: GammaExponent = %s\n",GammaExponent))
}
# Fraction of the Markov chain to consider as burn-in
if(any(input[,1]=="burnFrac")){
      burnFrac <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="burnFrac"),2])))
      if(is.na(burnFrac)){
            burnFrac <- 0.5
            cat(sprintf("NOTE: Typo in inputfile? BurnIn fraction is set to default: burnFrac  = %s\n",burnFrac))
      }
} else {
      burnFrac <- 0.5
      cat(sprintf("NOTE: BurnIn fraction is set to default: burnFrac = %s\n",burnFrac))
}

### PROPOSAL FUNCTION
# burn-in of the internal Gibbs sampler (rtmvnorm)
if(any(input[,1]=="burnGibbs")){
      burnGibbs <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="burnGibbs"),2])))
      if(is.na(burnGibbs)){
            burnGibbs <- 1000
            cat(sprintf("NOTE: Typo in inputfile? Proposal: burnIn for the internal Gibbs sampler in rtmvnorm is set to default: burnGibbs  = %s\n",burnGibbs))
      }
} else {
      burnGibbs <- 1000
      cat(sprintf("NOTE: Proposal: burnIn for the internal Gibbs sampler in rtmvnorm is set to default: burnGibbs = %s\n",burnGibbs))
}
# thinning parameter of the internal Gibbs sampler (rtmvnorm)
if(any(input[,1]=="thinGibbs")){
      thinGibbs <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="thinGibbs"),2])))
      if(is.na(thinGibbs)){
            thinGibbs  <- 20
            cat(sprintf("NOTE: Typo in inputfile? Proposal: Thinning for the internal Gibbs sampler in rtmvnorm is set to default: Thinning  = %s\n",thinGibbs))
      }
} else {
      thinGibbs <- 20
      cat(sprintf("NOTE: Proposal: Thinning for the internal Gibbs sampler in rtmvnorm is set to default: Thinning = %s\n",thinGibbs))
}

### TERMINATION CRITERION
# termination criterion: ESS
if(any(input[,1]=="TerminationCriterion")){
      TerminationCriterion <- as.logical(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="TerminationCriterion"),2])))
      if(is.na(TerminationCriterion)){
            TerminationCriterion <- TRUE
            cat("NOTE: Typo in inputfile? Termination of QPuB run based on effective sample size.\n")
      }
} else {
      TerminationCriterion <- TRUE
      cat("NOTE: Termination of QPuB run based on effective sample size.\n")
}
# Iteration when to start checking for convergence
if(TerminationCriterion==TRUE){
      if(any(input[,1]=="ESSiter")){
            ESSiter <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="ESSiter"),2])))
            if(is.na(ESSiter)){
                  ESSiter <- 10000
                  cat(sprintf("NOTE: Typo in inputfile? Iteration when to start checking for convergence is set to default: ESSiter= %s\n",ESSiter))
            }
      } else {
            ESSiter <- 10000
            cat(sprintf("NOTE: Iteration when to start checking for convergence is set to default: ESSiter= %s\n",ESSiter))
      }
}
# confidence level for multESS routine
if(TerminationCriterion==TRUE){
      if(any(input[,1]=="ConfLevel")){
            ConfLevel <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="ConfLevel"),2])))
            if(is.na(ConfLevel)){
                  ConfLevel <- 0.05
                  cat(sprintf("NOTE: Typo in inputfile? ESS Confidence level is set to default: ConfLevel = %s\n",ConfLevel))
            }
      } else {
            ConfLevel <- 0.05
            cat(sprintf("NOTE: ESS Confidence level is set to default: ConfLevel = %s\n",ConfLevel))
      }
}
# tolerance of multESS routine
if(TerminationCriterion==TRUE){
      if(any(input[,1]=="Tolerance")){
            Tolerance <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="Tolerance"),2])))
            if(is.na(Tolerance)){
                  Tolerance <- 0.05
                  cat(sprintf("NOTE: Typo in inputfile? ESS Tolerance is set to default: Tolerance = %s\n",Tolerance))
            }
      } else {
            Tolerance <- 0.05
            cat(sprintf("NOTE: ESS Tolerance is set to default: Tolerance = %s\n",Tolerance))
      }
}
# termination criterion: maximum number of iterations reached
if(TerminationCriterion==FALSE){
      if(any(input[,1]=="MaxIter")){
            MaxIter <- as.numeric(gsub(pattern="[[:space:]]", replacement="", x=as.vector(input[which(input[,1]=="MaxIter"),2])))
            if(is.na(MaxIter)){
                  MaxIter <- 10**6
                  cat(sprintf("NOTE: Typo in inputfile? Maximum iteration number for MCMC is set to default: MaxIter= %s\n",MaxIter))
            }
      } else {
            MaxIter <- 10**6
            cat(sprintf("NOTE: Maximum iteration number for MCMC is set to default: MaxIter= %s\n",MaxIter))
      }
}

print(input)

####################################################################################################
#### INPUT DATA

cat("\nREADING INPUT DATA...\n")

if(file.exists(file.path(paste(infol, 'data', sep='/'), fsep = .Platform$file.sep))){
      filenames <- list.files(file.path(paste(infol, 'data', sep='/'), fsep = .Platform$file.sep), pattern="*.csv", full.names=TRUE)
      if(length(filenames)==0){
            stop("Cannot find files in folder 'data'. The signal intensities have to be given as .csv files in the folder 'data'. Abort.")
      } else {
		dat <- lapply(filenames, function(x){ read.csv(file=x, header=T, check.names=FALSE) })
		times <- as.numeric(colnames(dat[[1]])[-1])
      }
} else {
      stop("Cannot find folder 'data'. The signal intensities have to be given as .csv files in the folder 'data'. Abort.")
}
print(dat)

# check data format
tmp_nrow <- dim(dat[[1]])[1]
tmp_ncol <- dim(dat[[1]])[2]
for(rp in 1:length(dat)){
      if(any(is.na(dat[[rp]]))){ stop(sprintf("No NA allowed in the data. Please check file %s!",rp)) }
      if(tmp_nrow != dim(dat[[rp]])[1]){ stop("Replicates don't have matching number of rows!") }
      if(tmp_ncol != dim(dat[[rp]])[2]){ stop("Replicates don't have matching number of columns!") }
      if(class(dat[[rp]][,1])!="factor"){ stop("First column must contain peptide sequences!") }
      for(clm in 2:tmp_ncol){
            if(!is.numeric(dat[[rp]][,clm])){ stop("All columns but the first must contain numeric values!") }
      }
}

####################################################################################################
#### INPUT TIMEPOINTS

cat("\nSETTING TIMEPOINTS FOR COMPARISON...\n")

if(exists('tim')){
      if(file.exists(file.path(paste(infol, tim, sep='/'), fsep = .Platform$file.sep))){
            cat("\nReading in timepoint file:\n")
            time_comp <- read.table(file.path(paste(infol, tim, sep='/'), fsep = .Platform$file.sep), header=FALSE)
      } else {
            stop(sprintf("Cannot find file %s. Abort.",tim))
      }
} else {
      # set default for timepoints to use: successive timepoints
      cat("Setting timepoints to default: comparing successive timepoints.\n")
      numT <- length(times) # number of timepoints
      TIM <- as.vector(rbind(seq(0,numT-1),seq(1,numT)))
      time_comp <- matrix(ncol=2, byrow=TRUE,
                     data=head(TIM,n=length(TIM)-2))
}
print(time_comp)

####################################################################################################
#### INPUT TITRATION

cat("\nREADING INPUT TITRATION...\n")

if(exists('titr')){
      if(file.exists(file.path(paste(infol,titr,sep='/'), fsep = .Platform$file.sep))){
            subs_titr <- read.table(file.path(paste(infol, titr, sep='/'), fsep = .Platform$file.sep), sep=',', header=TRUE)
			cat(sprintf("...from file %s\n", titr))
            print(subs_titr)
      } else {
            stop(sprintf("Cannot find file %s. Abort.",titr))
      }
} else {
      cat("No titration file given.\n")
}

} # end file
