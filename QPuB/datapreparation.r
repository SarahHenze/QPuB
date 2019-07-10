####################################################################################################
# DATAPREPARATION.R FOR QPUB PACKAGE
# SCALES DATA
# ASSIGNS CLEAVAGE AND SPLICE SITES OF PEPIDES
# CREATES POSITION PROBABILITY MATRIX
# EXTRACTS MASS SPECTROMETRY SIGNALS FROM DATA
####################################################################################################

{

####################################################################################################

cat("\nPREPARING DATA...\n")

####################################################################################################
# DEFINE GLOBAL VARIABLES
numD <- length(data) # number of digestions
numR <- length(data[[1]]) # number of replicates
numT <- dim(data[[1]][[1]])[2]-1  # number of timepoints
numCT <- dim(timepoints)[1] # number of comparisons of timepoints
numP <- dim(data[[1]][[1]])[1]-1 # number of products
substrate <- as.vector(unlist(data[[1]][[1]][1,1])) # sequence of substrate
numA <- nchar(substrate) # number of amino acids in the substrate
numPlen=c(rep(NA,numP))
for (i in 1:numP)
{
    numPlen[i] = nchar(as.vector(unlist(data[[1]][[1]][i,1])))

}

################################################################################
# FILTER DATA

if (FILTER == TRUE)
{

threshold <- Threshold # delete all products which have a signal below the threshold and do not contribute much
for(dg in 1:numD){
    minmax <- matrix(0,nrow=numR,ncol=numP)
    for(rp in 1:numR){
        minmax[rp,] <- apply(data[[dg]][[rp]][-1,-1], MARGIN=1, max)
    }
    throw <- c(which(apply(minmax, MARGIN=2, min) <= threshold))
    cat("The following peptides are excluded: ",throw,"\n")
    for(rp in 1:numR){
        data[[dg]][[rp]] <- data[[dg]][[rp]][-throw,]
    }
}
numP <- dim(data[[1]][[1]])[1]-1 # updated number of product

}

####################################################################################################
# PLOT INPUT DATA

if (!file.exists(sprintf('plots_input_%s',ID))){
      dir.create(sprintf('plots_input_%s',ID))
}

COLORS = rainbow(numR)
for(dg in 1:numD){ #TODO testen mit numD>1!!!

      ### SUBSTRATE
      png(paste(sprintf('plots_input_%s',ID), 'substrate.png', sep='/'))
      par(mar=c(5,7,5,3)) # bottom, left, top, and right

      ymin = min(data[[dg]][[1]][1,-1])
      ymax = max(data[[dg]][[1]][1,-1])
      for(rp in 1:numR){
            if(min(data[[dg]][[rp]][1,-1]) < ymin){ ymin = min(data[[dg]][[rp]][1,-1]) }
            if(max(data[[dg]][[rp]][1,-1]) > ymax){ ymax = max(data[[dg]][[rp]][1,-1]) }
      }
      for(rp in 1:numR){
            plot(seq(0,numT-1,1), data[[dg]][[rp]][1,-1], ylim=c(ymin,ymax), xlab='time', ylab='', main='substate', type='l', lwd=2, col=COLORS[rp], las=1)
            par(new=TRUE)
      }
      title(ylab = "signal intensity", cex.lab = 1, line = 5)
      dev.off()

      ### PRODUCTS
      for(pep in 1:numP){

            # png(normalizePath(paste(sprintf('plots_kinetics_%s',ID), sprintf('peptide%s.png',pep), sep='/'), winslash="\\", mustWork=TRUE))
            png(paste(sprintf('plots_input_%s',ID), sprintf('peptide%s.png',pep), sep='/'))
            par(mar=c(5,7,5,3)) # bottom, left, top, and right

            ymin = min(data[[dg]][[1]][pep+1,-1])
            ymax = max(data[[dg]][[1]][pep+1,-1])
            for(rp in 1:numR){
                  if(min(data[[dg]][[rp]][pep+1,-1]) < ymin){ ymin = min(data[[dg]][[rp]][pep+1,-1]) }
                  if(max(data[[dg]][[rp]][pep+1,-1]) > ymax){ ymax = max(data[[dg]][[rp]][pep+1,-1]) }
            }
            for(rp in 1:numR){
                  plot(seq(0,numT-1,1), data[[dg]][[rp]][pep+1,-1], ylim=c(ymin,ymax), xlab='time', ylab='', main=paste('peptide ',pep), type='l', lwd=2, col=COLORS[rp], las=1) # products
                  par(new=TRUE)
            }
            title(ylab = "signal intensity", cex.lab = 1, line = 5)
            dev.off()
      } # end pep
} # end dg
# stop('absicht')

####################################################################################################
#### DEFINE MECHANISM OF ENZYME

# define amino acid position where to start fitting
if(enzyme=='exopep'){
      a0 <- numA - min(nchar(levels(data[[1]][[1]][,1]))) +1
} else if(enzyme=='endopep'){
      a0 <- 1
}

####################################################################################################
# ASSIGN PEPTIDE ORIGIN IN SUBSTRATE
checkForsubstrings <- function(Str,substrate,k){

    k2 = numeric()
    for(i in 1:(length(k))){
        shorterSubstrate = paste(strsplit(substrate,"")[[1]][-c(1:k[i])],sep="",collapse="")
        tmp = gregexpr(pattern=Str,shorterSubstrate,perl=TRUE)[[1]]
        if(tmp[1]>0){k2 = c(k2,tmp+k[i])}
    }

    k2 = unique(k2)
    return(k2)
}

positions <- list()
for(dg in 1:numD){

      seq = as.vector(data[[dg]][[1]][-1,1])
      positions[[dg]] <- list()

      for(i in 1:length(seq)){

          l = nchar(seq[i])
          k = grep(seq[i],substrate)

          # PCP found
          if(length(k)>0){

            ind = gregexpr(pattern=seq[i],substrate)
            additionalIndices = checkForsubstrings(seq[i],substrate,ind[[1]])
            ind = sort(unique(c(ind[[1]],additionalIndices)))

            positions[[dg]][[i]] = matrix(NA,length(ind),2)

            for(j in 1:length(ind)){
            p1 = ind[j]
            p4 = p1+nchar(seq[i])-1
            positions[[dg]][[i]][j,] = c(p1,p4)
            }

          }
          else{ # no PCP found, check for all possible PSP with one splicing event

             # split peptide to P vector
             P = strsplit(seq[i],split="")[[1]]

             # get permutations of length l
             x = c(1:l)
             y = c(1:l)

             z = as.vector(outer(x,y,paste,sep="_"))
             q = matrix(NA,length(z),2)

             for(j in 1:length(z)){
                 q[j,] = as.numeric(strsplit(z[j],split="_")[[1]])
             }

             qs = apply(q,1,sum)
             k = which(qs==l)
             q = q[k,]

             # generate all strings for searches
             S = matrix(NA,dim(q)[1],2)
             for(j in 1:dim(q)[1]){
                 S[j,1] = paste(P[1:q[j,1]],sep="",collapse="")
                 S[j,2] = paste(P[(q[j,1]+1):l],sep="",collapse="")
             }

             # search each entry in prot for the two corresponding fragments and extract acc and positions
             positions[[dg]][[i]] = numeric()

             for(j in 1:dim(S)[1]){

                 psp <- list()
                 res1 = which((grepl(pattern=S[j,1],x=substrate)==TRUE))
                 res2 = which((grepl(pattern=S[j,2],x=substrate)==TRUE))

                 if((length(res1)>0)&(length(res2)>0)){

                     ind = gregexpr(pattern=S[j,1],substrate,perl=TRUE)
                     additionalIndices = checkForsubstrings(S[j,1],substrate,ind[[1]])
                     ind = sort(unique(c(ind[[1]],additionalIndices)))

                     n1 = rep(NA,length(ind))
                     n2 = rep(NA,length(ind))

                     for(r in 1:length(ind)){
                         n1[r] = ind[r]
                         n2[r] = n1[r]+nchar(S[j,1])-1
                     }

                     ind = gregexpr(pattern=S[j,2],substrate,perl=TRUE)
                     additionalIndices = checkForsubstrings(S[j,2],substrate,ind[[1]])
                     ind = sort(unique(c(ind[[1]],additionalIndices)))

                     n3 = rep(NA,length(ind))
                     n4 = rep(NA,length(ind))

                     for(r in 1:length(ind)){
                         n3[r] = ind[r]
                         n4[r] = n3[r]+nchar(S[j,2])-1
                     }

                     # get all internal combinations
                     counter = 0
                     z = as.vector(outer(n2,n3,paste,sep="_"))
                     tmp = matrix(NA,length(z),4)
                     for(n in 1:length(n1)){
                         for(m in 1:length(n3)){
                             p1 = n1[n]
                             p2 = n2[n]
                             p3 = n3[m]
                             p4 = n4[m]
                             counter = counter+1
                             tmp[counter,] = c(p1,p2,p3,p4)
                         }
                     }
                     positions[[dg]][[i]] = rbind(positions[[dg]][[i]],tmp)

                 } # end combination found
             } # end j
          } # end PSP
      } # end i
} # end dg

####################################################################################################
# POSITION PROBABILITY MATRIX
ppm = list()
for(dg in 1:numD){ # over digestions
      ppm[[dg]] = matrix(0,dim(data[[dg]][[1]])[1]-1,numA)
      for(pep in 1:numP){ # over products
            for(aa in 1:numA){ # over amino acids

                  if(dim(positions[[dg]][[pep]])[2]==2){ # PCP
                        for(poss in 1:dim(positions[[dg]][[pep]])[1]){ # over all possibilities
                              if((positions[[dg]][[pep]][poss,1]<=aa)&(aa<=positions[[dg]][[pep]][poss,2])){
                                    possib = dim(positions[[dg]][[pep]])[1]
                                    fraction = 1/possib
                                    ppm[[dg]][pep,aa] = fraction
                              }
                        }
                  } else if(dim(positions[[dg]][[pep]])[2]==4){ # PSP

                        SR1 = matrix(unique(positions[[dg]][[pep]][,c(1,2)]),ncol=2) # different splice reactants 1
                        SR2 = matrix(unique(positions[[dg]][[pep]][,c(3,4)]),ncol=2) # different splice reactants 2

                        frac1 = 1/dim(SR1)[1]
                        frac2 = 1/dim(SR2)[1]

                        for(poss_sr1 in 1:dim(SR1)[1]){ # over all possibilities of splice reactants 1
                                    if((SR1[poss_sr1,1]<=aa)&(aa<=SR1[poss_sr1,2])){
                                          ppm[[dg]][pep,aa] = ppm[[dg]][pep,aa]+frac1
                                    }
                        }
                        for(poss_sr2 in 1:dim(SR2)[1]){ # over all possibilities of splice reactants 2
                                    if((SR2[poss_sr2,1]<=aa)&(aa<=SR2[poss_sr2,2])){
                                          ppm[[dg]][pep,aa] = ppm[[dg]][pep,aa]+frac2
                                    }
                        }
                  }
            } # end aa
      } # end pep
} # end dg

####################################################################################################
# EXTRACT SIGNALS FROM INPUT DATA
signalP <- list()
signalF <- list()
for(dg in 1:numD){ # over digestions
      signalP[[dg]] <- list()
      signalF[[dg]] <- list()
      for(rp in 1:numR){ # over replicates

            signalS <- as.numeric(data[[dg]][[rp]][1,2:dim(data[[dg]][[rp]])[2]]) # signals substrate
            signalF1 <- data[[dg]][[rp]][2:dim(data[[dg]][[rp]])[1],2:dim(data[[dg]][[rp]])[2]] # signals products

            # signalP signals sum of products
            # signalF signal differences of products over time
            signalP[[dg]][[rp]] <- signalS[1]-signalS[1]
            signalF[[dg]][[rp]] = matrix(signalF1[,1]-signalF1[,1],ncol=1)

            for(tp in 1:numCT){
                  signalP[[dg]][[rp]] <- cbind(signalP[[dg]][[rp]], signalS[timepoints[tp,1]+1] - signalS[timepoints[tp,2]+1])
                  signalF[[dg]][[rp]] <- cbind(signalF[[dg]][[rp]], matrix(signalF1[,timepoints[tp,2]+1]-signalF1[,timepoints[tp,1]+1],ncol=1))
            }

      } # end rp
} # end dg

####################################################################################################
# SCALE SIGNALS

log10_ceiling <- function(x){ 10^(ceiling(log10(x))) } # get order of magnitude, ceiling

# general scaling to avoid large values
order_data_new = mean(signalF[[1]][[1]])
if(order_data_new>1000){
	scalefac2 = 1/log10_ceiling(order_data_new)*100
      for(dg in 1:numD){
            for(rp in 1:numR){
                  signalP[[dg]][[rp]] = signalP[[dg]][[rp]] * scalefac2
                  signalF[[dg]][[rp]] = signalF[[dg]][[rp]] * scalefac2
            }
      }
      cat(sprintf("Data scaled by %s\n",scalefac2))
}

# scale products to fit substrate
if(scaleprod==TRUE){

      sig = list()
      for(tp in 1:(numCT+1)){ # over timepoints
            sig[[tp]] = rep(0,numA)
            for(aa in a0:numA){ # over amino acids
                  siga = sum((signalF[[1]][[1]][,tp]) * ppm[[1]][,aa])
                  sig[[tp]][aa] = sig[[tp]][aa] + siga
            }
      }
      sig = unlist(lapply(sig, mean))
      quot = sig/as.numeric(abs(signalP[[1]][[1]]))
      scalefac1 = 1/max(quot[!is.infinite(quot)], na.rm=TRUE)

      for(dg in 1:numD){
            for(rp in 1:numR){
                  signalF[[dg]][[rp]] = signalF[[dg]][[rp]] * scalefac1
            }
      }
      cat(sprintf("Products scaled by %s\n",scalefac1))

} else { scalefac1 = 1 }


} # end file
