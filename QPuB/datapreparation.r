####################################################################################################
# DATAPREPARATION.R FOR QPUB PACKAGE
# DEFINE GLOBAL VARIABLES
# DEFINE MECHANISM OF ENZYME
# FILTER DATA
# PLOT INPUT SIGNAL KINETICS
# ASSIGN CLEAVAGE AND SPLICE SITES OF PEPIDES
# CREATE POSITION PROBABILITY MATRIX
# EXTRACT MASS SPECTROMETRY SIGNALS FROM DATA
# SCALE DATA
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

cat("\nPREPARING DATA...\n")

####################################################################################################    
### DEFINE GLOBAL VARIABLES
numR <- length(dat) # number of replicates
numT <- length(times) # number of timepoints
numCT <- dim(time_comp)[1] # number of comparisons of timepoints
numP <- dim(dat[[1]])[1]-1 # number of products
substrate <- as.vector(unlist(dat[[1]][1,1])) # sequence of substrate
numA <- nchar(substrate) # number of amino acids in the substrate

####################################################################################################
# MECHANISM OF ENZYME

# define amino acid position where to start fitting
if(enzyme=='exopep'){
      a0 <- numA - min(nchar(levels(dat[[1]][,1]))) +1
} else if(enzyme=='endopep'){
      a0 <- 1
}

################################################################################
# FILTER DATA

# delete all products which have a signal intensity below the threshold and do not contribute much

if (filterprod == TRUE){

      minmax <- matrix(0,nrow=numR,ncol=numP)
      for(rp in 1:numR){
            minmax[rp,] <- apply(dat[[rp]][-1,-1], MARGIN=1, max)
      }
      throw <- c(which(apply(minmax, MARGIN=2, min) <= threshold))
      if(length(throw)!=0){
            cat("The following peptides are excluded: ",throw,"\n")
            for(rp in 1:numR){
              dat[[rp]] <- dat[[rp]][-(throw+1),]
            }
      }
      numP <- dim(dat[[1]])[1]-1 # updated number of product
      
      # save new data to file
      FILES = lapply(basename(filenames),function(x){unlist(strsplit(x,'[.]'))[1]})
      for(rp in 1:numR){
            write.csv(dat[[rp]], sprintf('%s_filtered.csv',FILES[rp]), row.names=FALSE, quote=FALSE)
      }
}

####################################################################################################
# PLOT INPUT SIGNAL KINETICS

if (!file.exists('plots_inputsignals')){
      dir.create('plots_inputsignals')
}

COLORS = rainbow(numR)

### SUBSTRATE
plotfct(file.path(paste('plots_inputsignals', 'substrate', sep='/'), fsep = .Platform$file.sep), width=5, height=5)
par(mar=c(5,7,5,3)) # bottom, left, top, and right

ymin = min(dat[[1]][1,-1], na.rm=T)
ymax = max(dat[[1]][1,-1], na.rm=T)
for(rp in 1:numR){
      if(min(dat[[rp]][1,-1], na.rm=T) < ymin){ ymin = min(dat[[rp]][1,-1], na.rm=T) }
      if(max(dat[[rp]][1,-1], na.rm=T) > ymax){ ymax = max(dat[[rp]][1,-1], na.rm=T) }
}
for(rp in 1:numR){
      plot(times, dat[[rp]][1,-1], ylim=c(ymin,ymax), xlab='', ylab='', main='', type='l', lwd=2, col=COLORS[rp], las=1)
      par(new=TRUE)
}
title(xlab ='time', cex.lab = 1, line = 3)
title(ylab = 'signal intensity', cex.lab = 1, line = 5)
title(main = 'substrate', cex.lab = 1, line = 2)
dev.off()

### SUBSTRATE DEGRADED 
plotfct(file.path(paste('plots_inputsignals', 'substrate_degraded', sep='/'), fsep = .Platform$file.sep), width=5, height=5)
par(mar=c(5,7,5,3)) # bottom, left, top, and right

ymin = min(dat[[1]][1,-1], na.rm=T)
ymax = max(dat[[1]][1,-1], na.rm=T)
for(rp in 1:numR){
      if(min(dat[[rp]][1,2]-dat[[rp]][1,-1], na.rm=T) < ymin){ ymin = min(dat[[rp]][1,2]-dat[[rp]][1,-1], na.rm=T) }
      if(max(dat[[rp]][1,2]-dat[[rp]][1,-1], na.rm=T) > ymax){ ymax = max(dat[[rp]][1,2]-dat[[rp]][1,-1], na.rm=T) }
}
for(rp in 1:numR){
      plot(times, dat[[rp]][1,2]-dat[[rp]][1,-1], ylim=c(ymin,ymax), xlab='', ylab='', main='', type='o', lwd=2, col=COLORS[rp], las=1)
      par(new=TRUE)
}
title(xlab = 'time', cex.lab = 1, line = 3)
title(ylab = 'signal intensity', cex.lab = 1, line = 5)
title(main = 'substate degraded', cex.lab = 1, line = 2)
dev.off()

### PRODUCTS   
for(pep in 1:numP){
      
      plotfct(file.path(paste('plots_inputsignals', sprintf('peptide%s',pep), sep='/'), fsep = .Platform$file.sep), width=5, height=5)
      par(mar=c(5,7,5,3)) # bottom, left, top, and right

      ymin = min(dat[[1]][pep+1,-1], na.rm=T)
      ymax = max(dat[[1]][pep+1,-1], na.rm=T)
      for(rp in 1:numR){
            if(min(dat[[rp]][pep+1,-1], na.rm=T) < ymin){ ymin = min(dat[[rp]][pep+1,-1]) }
            if(max(dat[[rp]][pep+1,-1], na.rm=T) > ymax){ ymax = max(dat[[rp]][pep+1,-1]) }
      }
      for(rp in 1:numR){
            plot(times, dat[[rp]][pep+1,-1], ylim=c(ymin,ymax), xlab='', ylab='', main='', type='o', lwd=2, col=COLORS[rp], las=1)
            par(new=TRUE)
      }
      title(xlab = 'time', cex.lab = 1, line = 3)
      title(ylab = 'signal intensity', cex.lab = 1, line = 5)
      title(main = paste('peptide ',pep), cex.lab = 1, line = 2)
      dev.off()
} # end pep

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

seq = as.vector(dat[[1]][-1,1])
positions <- list()

for(i in 1:length(seq)){

    l = nchar(seq[i])
    k = grep(seq[i],substrate)

    # PCP found
    if(length(k)>0){

      ind = gregexpr(pattern=seq[i],substrate)
      additionalIndices = checkForsubstrings(seq[i],substrate,ind[[1]])
      ind = sort(unique(c(ind[[1]],additionalIndices)))

      positions[[i]] = matrix(NA,length(ind),2)

      for(j in 1:length(ind)){
      p1 = ind[j]
      p4 = p1+nchar(seq[i])-1
      positions[[i]][j,] = c(p1,p4)
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
       positions[[i]] = numeric()

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
               positions[[i]] = rbind(positions[[i]],tmp)

           } # end combination found
       } # end j
    } # end PSP
} # end i

# save peptide number, sequence and positions in csv
identifier <- list()
for(prod in 1:numP){
      # unique assignment
      var=1
      A = cbind(c(prod),as.character(dat[[1]][prod+1,1]),paste(positions[[prod]][var,], collapse='_'))
      # non-unique assignment
      VAR = dim(positions[[prod]])[1]
      if(VAR>1){
            for(var in 2:VAR){
                  A = cbind(A,paste(positions[[prod]][var,], collapse='_'))
            }
      }
      identifier[prod] <- list(A)
}
lapply(identifier, function(x) write.table( data.frame(unlist(x),check.names=F), 'identifier.csv', row.names=FALSE, col.names=FALSE, append= T, sep=','))

####################################################################################################
# POSITION PROBABILITY MATRIX
ppm = list()
ppm = matrix(0,dim(dat[[1]])[1]-1,numA)
for(pep in 1:numP){ # over products
      for(aa in 1:numA){ # over amino acids

            if(dim(positions[[pep]])[2]==2){ # PCP
                  for(poss in 1:dim(positions[[pep]])[1]){ # over all possibilities
                        if((positions[[pep]][poss,1]<=aa)&(aa<=positions[[pep]][poss,2])){
                              possib = dim(positions[[pep]])[1]
                              fraction = 1/possib
                              ppm[pep,aa] = fraction
                        }
                  }
            } else if(dim(positions[[pep]])[2]==4){ # PSP

                  SR1 = matrix(unique(positions[[pep]][,c(1,2)]),ncol=2) # different splice reactants 1
                  SR2 = matrix(unique(positions[[pep]][,c(3,4)]),ncol=2) # different splice reactants 2

                  frac1 = 1/dim(SR1)[1]
                  frac2 = 1/dim(SR2)[1]

                  for(poss_sr1 in 1:dim(SR1)[1]){ # over all possibilities of splice reactants 1
                              if((SR1[poss_sr1,1]<=aa)&(aa<=SR1[poss_sr1,2])){
                                    ppm[pep,aa] = ppm[pep,aa]+frac1
                              }
                  }
                  for(poss_sr2 in 1:dim(SR2)[1]){ # over all possibilities of splice reactants 2
                              if((SR2[poss_sr2,1]<=aa)&(aa<=SR2[poss_sr2,2])){
                                    ppm[pep,aa] = ppm[pep,aa]+frac2
                              }
                  }
            }
      } # end aa
} # end pep

####################################################################################################
# EXTRACT SIGNALS FROM INPUT DATA
signalP <- list()
signalF <- list()
for(rp in 1:numR){ # over replicates

      signalS <- as.numeric(dat[[rp]][1,2:dim(dat[[rp]])[2]]) # signals substrate
      signalF1 <- dat[[rp]][2:dim(dat[[rp]])[1],2:dim(dat[[rp]])[2]] # signals products

      # signalP signals sum of products
      # signalF signal differences of products over time
      signalP[[rp]] <- signalS[1]-signalS[1]
      signalF[[rp]] = matrix(signalF1[,1]-signalF1[,1],ncol=1)

      for(tp in 1:numCT){
            signalP[[rp]] <- cbind(signalP[[rp]], signalS[time_comp[tp,1]+1] - signalS[time_comp[tp,2]+1])
            signalF[[rp]] <- cbind(signalF[[rp]], matrix(signalF1[,time_comp[tp,2]+1]-signalF1[,time_comp[tp,1]+1],ncol=1))
      }

} # end rp

####################################################################################################
# SCALE SIGNALS

log10_ceiling <- function(x){ 10^(ceiling(log10(x))) } # get order of magnitude, ceiling

# general scaling to avoid large values
ordermag = mean(c(signalF[[1]][[1]], signalS))
if(ordermag>1000){
	scale_mag = 1/log10_ceiling(ordermag)*100
        for(rp in 1:numR){
              signalP[[rp]] = signalP[[rp]] * scale_mag
              signalF[[rp]] = signalF[[rp]] * scale_mag
        }
      cat(sprintf("Data scaled by %s to avoid large orders of magnitude.\n",scale_mag))

      # sigma_start <- sigma_start * 0.2 #TODO mean(signalP)*0.5?
}

# scale products to fit order of magnitude of substrate
if(scaleprod==TRUE){

      sig = list()
      for(tp in 1:(numCT+1)){ # over timepoints
            sig[[tp]] = rep(0,numA)
            for(aa in a0:numA){ # over amino acids
                  siga = sum((signalF[[1]][,tp]) * ppm[,aa])
                  sig[[tp]][aa] = sig[[tp]][aa] + siga
            }
      }
      sig = unlist(lapply(sig, mean))
      quot = sig/as.numeric(abs(signalP[[1]]))
      scale_prod = 1/max(quot[!is.infinite(quot)], na.rm=TRUE)

        for(rp in 1:numR){
              signalF[[rp]] = signalF[[rp]] * scale_prod
        }
      cat(sprintf("Products scaled by %s to match order of magnitude of the substrate.\n",scale_prod))

} else { scale_prod = 1 }

} # end file
