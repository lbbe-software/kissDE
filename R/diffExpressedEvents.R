### R code from vignette source 'KissCS_Sknsh.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: parameters
###################################################

diffExpressedEvents <- function(n,nr,file) {

option.Plot <- FALSE
if (option.Plot == TRUE) {
system("mkdir Figures") 
}

dataName <- "data name"

# n <- 2 #number of biological conditions
# nr <- c(2,2) #number of replicates for each condition
# cond <- rep(c("cond1","cond2"),nr)
cond <- rep(sapply(1:n, FUN=function(i){paste("cond",i,sep="")}), nr)

# k_kissplice <- 41
# data.file <-"results_wgEncodeCshlLongRnaSeqSknshCellPapFastqRd1Rep3_10M_wgEncodeCshlLongRnaSeqSknshCellPapFastqRd1Rep4_10M_wgEncodeCshlLongRnaSeqSknshraCellPapFastqRd1Rep1_10M_wgEncodeCshlLongRnaSeqSknshraCellPapFastqRd1Rk41_coherents_type_1.txt" # output of pretraitement.py


###################################################
### code chunk number 2: load libraries
###################################################
# if("DESeq"%in%installed.packages()){
#   library(DESeq)
# } else {
#   source("http://bioconductor.org/biocLite.R")
#   biocLite("DESeq")
#   library(DESeq)
# }

# options("repos" = c(CRAN = "http://cran.r-project.org/"))
# ## set a default repo. avoid error if no cran mirror by default in the user session
# if("aod"%in%installed.packages()){
#   library(aod)
# } else{
# install.packages("aod")
# library(aod)
# }
# # aod: negbin function
# if("xtable"%in%installed.packages()){
#   library(xtable)
# } else{
# install.packages("xtable")
# library(xtable)
# }
# if("DSS"%in%installed.packages()){
#   library(DSS)
# } else {
#   source("http://bioconductor.org/biocLite.R")
#   biocLite("DSS")
#   library(DSS)
# }

# if("glmnet"%in%installed.packages()){
#   library(glmnet)
# } else{
# install.packages("glmnet")
# library(glmnet)
# }

set.seed(10) # to help comparison between vignette
options(warn=-1) # supress the warning for the users


###################################################
### code chunk number 3: Read data
###################################################
# TODO: forbid runnning the vignette if only one condition
# on personal computer
dataCounts <- read.table(file,h=F,sep="\t")

# namesData <- c("ID", "Length")
# for(i in 1:n ) {
#   for( j in 1:nr[i] ) {
#       namesData <- c( namesData,paste("Cond",i,"_","R",j,sep="",collapse=""))
#     }
# }

# names(dataCounts) <- namesData
# # add path information
# dataCounts$Path <- gl( 2, 1, dim(dataCounts)[1], labels = c("UP", "LP") )




#####################################################################################
# TODO : test if a data or .fa (KisSplice output) is given and apply or not convertToCounts
namesData <- c("ID", "Length")
for (i in 1:n){
  for (j in 1:nr[i]){
      namesData <- c( namesData,paste("Cond",i,"_","R",j,sep = "",collapse = ""))
    }
}

convertToCounts <- function(fileName,namesData) {

  toConvert <- file(fileName, open = "r")
  line <- readLines(toConvert)
  index <- 1
  if (substr(line[index], start = 0, stop = 1) == '>') {
    len <- length(strsplit(line[index], "|", fixed = TRUE)[[1]])-6
  }
  lineCounts <- strsplit(line[index], "|", fixed = TRUE)[[1]][5:(5 + len)] # computes the right number of conditions and their counts
  eventCounts <- as.vector(as.numeric(lapply(lineCounts, function(x) strsplit(x, "_")[[1]][2])))
  # eventName <- c(strsplit(substr(line[index],start = 2,stop = length(strsplit(line[index],"|")[[1]])),"|T",fixed = TRUE)[[1]][1])
  eventName <- gsub("|","_",c(strsplit(substr(line[index],start = 2,stop = length(strsplit(line[index],"|")[[1]])),"|T",fixed = TRUE)[[1]][1]),fixed=TRUE)
  eventLength <- as.vector(as.numeric(strsplit(strsplit(substr(line[index],start = 2,stop = length(strsplit(line[index],"|")[[1]])),"|",fixed = TRUE)[[1]][4],"_")[[1]][4])) # parsing to get information to build the dataframe : bcc id, path length and counts /!\ I assumed that types are always noted as "|Type_X|" for the sake of the parsing simplicity, which is true in KisSplice outputs so far."
  events.df <- data.frame(eventName,eventLength,t(eventCounts))
  names(events.df) <- namesData
  index <- index+1
  while (index <= length(line)) {
     if (substr(line[index], start = 0, stop = 1) == '>') {lineCounts <- strsplit(line[index], "|", fixed=TRUE)[[1]][5:(5+len)] # computes the right number of conditions and their counts
       eventCounts <- as.vector(as.numeric(lapply(lineCounts, function(x) strsplit(x, "_")[[1]][2])))
       # eventName <- c(strsplit(substr(line[index],start = 2,stop = length(strsplit(line[index],"|")[[1]])),"|T",fixed = TRUE)[[1]][1])
       eventName <- gsub("|","_",c(strsplit(substr(line[index],start = 2,stop = length(strsplit(line[index],"|")[[1]])),"|T",fixed = TRUE)[[1]][1]),fixed=TRUE)
       eventLength <- as.vector(as.numeric(strsplit(strsplit(substr(line[index],start = 2,stop = length(strsplit(line[index],"|")[[1]])),"|",fixed = TRUE)[[1]][4],"_")[[1]][4])) # parsing to get information to build the dataframe : bcc id, path length and counts /!\ I assumed that types are always noted as "|Type_X|" for the sake of the parsing simplicity, which is true in KisSplice outputs so far."
       events.df.temp <- data.frame(eventName,eventLength,t(eventCounts))
       names(events.df.temp) <- namesData
       events.df <- rbind(events.df,events.df.temp)
     }
     index <- index + 1
  }
  return (events.df)
 }
dataCounts <- convertToCounts(file,namesData)
dataCounts$Path <- gl( 2, 1, dim(dataCounts)[1], labels = c("UP", "LP") )
#############################################################################################

# return(dataCounts)


###################################################
### code chunk number 4: Normalisation
###################################################
# Normalisation with DESeq
conds <- c()
for( i in 1:n ) {
  for( j in 1:nr[i] ) {
      conds <- c( conds,paste( "Cond", i, sep = "",collapse = "") )
    }
} 
cds <- newCountDataSet( dataCounts[ ,3:(3+length(conds)-1)], conds ) # create object
cdsSF <- estimateSizeFactors( cds )
## see : http://seqanswers.com/forums/showpost.php?p=16468&postcount=13
sizeFactors( cdsSF )
#Vincent: In case there are too few data, the sizeFactors are not estimated, and the normalisation fails. I introduce this variable: shouldWeNormalise. I think that it could be put a general parameter at some point.
shouldWeNormalise=sum(is.na(sizeFactors(cdsSF))) < 1
dim <- dim(dataCounts)[2]
dataCounts[ ,(dim+1):(dim+length(conds)) ] <- round(counts(cdsSF, normalized=shouldWeNormalise))
colnames(dataCounts)[(dim+1):(dim+length(conds))] <- paste(namesData[3:(3+sum(nr)-1)],"_Norm",sep="")


###################################################
### code chunk number 5: fig_hclust_raw
###################################################
if (option.Plot == TRUE){
plot( hclust(as.dist(1-cor(dataCounts[ ,3:(3+length(conds)-1)])),"ward") )  
}



###################################################
### code chunk number 6: fig_hclust_norm
###################################################
if (option.Plot == TRUE) {
plot( hclust(as.dist(1-cor(dataCounts[ ,(dim+1):(dim+length(conds))])),"ward") )
}

###################################################
### code chunk number 7: replicates
###################################################
if (option.Plot==TRUE) {
par( mfrow=c(n,1) )
for( i in 1:n ) {
  if ( nr[i] > 1 ) {
      plot( dataCounts[ ,paste("Cond",i,"_R1_Norm",sep="",collapse="")]+1, dataCounts[ ,paste("Cond",i,"_R2_Norm",sep="",collapse="")]+1, log="xy",
           xlim=c( 1, max( dataCounts[ ,paste( "Cond",i,"_R1_Norm",sep="",collapse="" ) ], dataCounts[ ,paste("Cond", i, "_R2_Norm", sep = "", collapse = "") ] ) ),
           xlab = paste( "Cond",i,"_R1_Norm",sep="",collapse="" ),
           ylim=c(1, max(dataCounts[ ,paste("Cond",i,"_R1_Norm",sep="",collapse="")], dataCounts[ ,paste("Cond",i,"_R2_Norm",sep="",collapse="")])),main = paste("Replicate 1 vs Replicate 2, ","Cond",i,sep="",collapse=""),
           ylab = paste("Cond",i,"_R2_Norm",sep="",collapse="")
           )
      #end plot
                                        #adding ev. with zero counts in on replicates
      zeroC.R1R2 <- dataCounts[which(dataCounts[ ,paste("Cond",i,"_R1_Norm",sep="",collapse="")] == 0 | dataCounts[ ,paste("Cond",i,"_R2_Norm",sep="",collapse="")]==0), ]
      
      points(zeroC.R1R2[ , paste("Cond",i,"_R1_Norm",sep="",collapse="") ]+1, zeroC.R1R2[ ,paste("Cond",i,"_R2_Norm",sep="",collapse="")]+1, col=3)
      abline(a=0, b=1, col=2, lty=2, lwd=2)
    }
}
}

###################################################
### code chunk number 8: MA-plot
###################################################
if (option.Plot == TRUE) {
dataM <- apply(dataCounts[ ,grepl("Norm",names(dataCounts))],1,sum) / sum(nr)
dataA <- (apply( dataCounts[ ,grepl("Norm",names(dataCounts)) & grepl( "Cond2",names(dataCounts)) ] ,1, sum ) / sum(nr[2])) / ( apply( dataCounts [ ,grepl("Norm",names(dataCounts)) & grepl("Cond1",names(dataCounts)) ],1,sum ) / sum( nr[2] ) )

plot(dataM, dataA, 
     log="xy", pch=".", xlab="Average", ylab="Minus (LOG ratio)", cex=3, 
     ylim=c(0.001, 1000), las=1,main=paste("MA-plot",unique(cond)[1],"vs.",unique(cond)[2],sep=" ",collapse=" "))
abline(h=1, col=2)
legend("topleft",c(unique(cond)[2]))
legend("bottomleft",c(unique(cond)[1]))
}



###################################################
### code chunk number 9: intra-group and inter-group-variance
###################################################
## Mean and variance over all conditions and replicates (normalized counts!) 


dataCounts$mn <- apply(dataCounts[ ,(dim+1):(dim+length(conds))], 1, mean)
dataCounts$var <- apply(dataCounts[ ,(dim+1):(dim+length(conds))], 1, var)
## correction term
nbAll <- sum(nr) # number of all observations in all groups
dataCounts$ct  <- apply(dataCounts[ ,(dim+1):(dim+length(conds))], 1, sum)^2/nbAll
#### BE CAREFULL!!!! Here only if two conditions computations ,
### not a big deal beacause only for the plot


## sum of squares between groups
dataCounts$ss <- apply(dataCounts[ ,(dim+1):(dim+length(conds)/n)],1,sum)^2/nr[1] + 
  apply(dataCounts[ ,((dim+1)+length(conds)/n):(dim+length(conds))],1,sum)^2/nr[2] 
## substract the correction term from the SS and divide by the degrees of 
df <- 1 # freedom(groups); here: df=2-1=1
dataCounts$varInter <- (dataCounts$ss - dataCounts$ct)/df
# intra-variability 
dataCounts$varC1 <- apply( dataCounts[ ,(dim+1):(dim+nr[1])], 1, var )
dataCounts$varC2 <- apply( dataCounts[ ,((dim+1)+nr[1]):(dim+nr[2]+nr[1])], 1, var )
dataCounts$varIntra <- apply( data.frame(dataCounts$varC1, dataCounts$varC2), 1, mean )
head(dataCounts)# here wrong if various conditions



###################################################
### code chunk number 10: intra-vs-inter
###################################################
if (option.Plot == TRUE) {
plot( x = dataCounts$varIntra, y = dataCounts$varInter, xlab = "Intra-variability", ylab = "Inter-variability", las = 1, log = "xy")
abline( a = 0, b = 1, col = 2, lty = 2, lwd = 2 )
}

###################################################
### code chunk number 11: ecdf
###################################################
if (option.Plot == TRUE) {
dataCounts$allNormCounts <- apply(dataCounts[ ,grepl("Norm",names(data))], 1, sum)
plot(ecdf(dataCounts$allNormCounts), verticals=TRUE)
}



###################################################
### code chunk number 12: snp-list
###################################################
# reduce data frame to the interesting columns
dataPart <- dataCounts[ ,c(1:2,which(grepl("_Norm",names(dataCounts)))) ] # -> when not want to filter data 
dataPart$Path <- gl( 2,1,dim(dataCounts)[1],labels=c("UP","LP") )#***********************

# modify the dataPart data.frame such that a row corresponds to a bubble i.e. the two alleles  
dataPart2 <- cbind( dataPart[ seq( 1, dim( dataPart )[1], 2 ), ], dataPart[ seq( 2, dim(dataPart)[1] , 2 ), grepl("Norm", names( dataPart) ) ] )
names(dataPart2)[3:(3+nbAll-1)] <- paste( "UP",names(dataPart2)[3:(3+nbAll-1)], sep="_" )
names(dataPart2)[(3+nbAll+1):(3+2*nbAll+1-1)] <- paste( "LP", names(dataPart2)[(3+nbAll+1):(3+2*nbAll+1-1)], sep="_" )

dataPart2 <- dataPart2[ ,c(-(3+nbAll))]
if (anyDuplicated(dataPart2[ ,1]) > 0) {
  dataPart2 <- dataPart2[!duplicated(as.character(dataPart2[ ,1])), ]
}

rownames(dataPart2)=as.character(dataPart2[ ,1])
# print(head(dataPart2))
#---------------------------------------------
# each list element is a table itself with columns:
# ID - biological condition - counts - path
# Parameters for function:
# df - name of data frame
# startPosColumn4Counts - number of column where the counts start (e.g. 4 for the 4th column)
# endPosCol4Counts - number of column where the counts end (function assumes that all columns between start and end are counts)
#---------------------------------------------
snptable <- function(df,startPosColumn4Counts, endPosCol4Counts){
    
  snpTab = data.frame(ID=rep(as.factor(df['ID']), endPosCol4Counts-startPosColumn4Counts+1),
    cond=as.factor(unlist(lapply(strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts],"_"),FUN=function(d){d[2]}))),
    counts=as.numeric(df[startPosColumn4Counts:endPosCol4Counts]),
    path=as.factor(unlist(lapply(strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts],"_"),FUN=function(d){d[1]}))),
    row.names = NULL)
  
  return(snpTab)
}

#---------------------------------------------

# snptable(snp.list[[1]], noRep=3, cond=c("Cond1","Cond2"), 4, 9)

# create list for the complete data set
allSNPtables  <- apply(dataPart2,1,snptable, startPosColumn4Counts = which(grepl("UP",names(dataPart2)))[1],endPosCol4Counts = ncol(dataPart2))

#example of the first element of the SNPList
allSNPtables[[1]]


###################################################
### code chunk number 13: DSS dispersion estimation
###################################################

dataNormCountsSNP <- as.matrix(dataPart2[ ,3:ncol(dataPart2)]) # the counts matrix
colnames(dataNormCountsSNP) <- 1:ncol(dataNormCountsSNP)
designs <- rep(c(1:(n*2)),c(nr,nr)) # the design matrix

dispData <- newSeqCountSet(dataNormCountsSNP, as.data.frame(designs))
dispData <-estDispersion(dispData)

dispDataMeanCond <- newSeqCountSet(dataNormCountsSNP, as.data.frame(designs))
#dispDataMeanCond <- estDispersion(dispData,trend=T)
#Vincent: I removed the option trend=T, since it was creating trouble with Gaurav's data
dispDataMeanCond <- estDispersion(dispData)


###################################################
### code chunk number 14: variance - mean - SNP level1
###################################################
#-----------------------------------------------------------
# compute mean and variance per SNP (instead of per allele)
#-----------------------------------------------------------
snp.mean.variance.df <- as.data.frame(cbind(apply(dataPart2[ ,which(grepl("_Norm",names(dataPart2)))],1,mean),apply(dataPart2[ ,which(grepl("_Norm",names(dataPart2)))],1,var)))
names(snp.mean.variance.df) <- c("Mean", "Variance")
rownames(snp.mean.variance.df) <- as.character(dataPart2[ ,1])

# estimate the dispersion parameter D of the Quasi-Poisson distribution
# by fitting equation (7) to the data
lm.D <- lm(snp.mean.variance.df$Variance ~ snp.mean.variance.df$Mean-1)
coef(lm.D)

## estimate the overdispersion parameter theta of the NB distritubion
## by fitting equation (10) to the data

modelNB <- Variance ~ Mean + 1/theta * Mean^2
nls.modelNB <- nls(modelNB, data=snp.mean.variance.df, start=list(theta=100))
coef(nls.modelNB)
phi <- 1/coef(nls.modelNB) # to be used as fixed parameter later on


#compute model fit
# log plot, so exclude 0 from the x values
x <- c(seq(0.1,1,0.1), seq(2,5000,1) )
yQP <- x*coef(lm.D)
yNB <- x + 1/coef(nls.modelNB) * x^2


###################################################
### code chunk number 15: SNPlevel3
###################################################
if (option.Plot == TRUE) {
plot(snp.mean.variance.df$Mean, snp.mean.variance.df$Variance, 
     xlab="Mean SNP count", 
     ylab="Variance SNP count",
     log="xy", las=1)
abline(a=0, b=1, col=2, lwd=2)
lines(x,yQP,col=3, lwd=2)
lines(x,yNB,col=6, lwd=2)
legend("topleft", c("Poisson","Quasi-Poisson", "Negative Binomial"), 
   text.col=c(2,3,6), box.lty=0);
}



###################################################
### code chunk number 16: function fitNBglmModelsDSSPhi
###################################################
fitNBglmModelsDSSPhi <- function(snpdata, phiDSS, phiDSScond, phiGlobal){
# S: simple, A: additive, I : interaction models
  # modèle de Poisson   
  nbglmS0 <- negbin(counts~cond, data=snpdata, random=~1, fixpar=list(3,0))
  nbglmA0 <- negbin(counts~cond + path, data=snpdata, random=~1, fixpar=list(4,0))
  nbglmI0 <- negbin(counts~cond * path, data=snpdata, random=~1, fixpar=list(5,0))  

  nbAnov0 <- anova(nbglmS0, nbglmA0, nbglmI0)
  nbAIC0 <- c(AIC(nbglmS0,k = log(nbAll))@istats$AIC, AIC(nbglmA0,k = log(nbAll))@istats$AIC, AIC(nbglmI0,k = log(nbAll))@istats$AIC)
  # singular.hessian:  true when fitting provided a singular hessian, indicating an overparamaterized model.
  nbSingHes0 <- c(nbglmS0@singular.hessian, nbglmA0@singular.hessian, nbglmI0@singular.hessian)
  # code: ‘code’ An integer (returned by ‘optim’) indicating why the optimization process terminated.
  nbCode0 <- c(nbglmS0@code, nbglmA0@code, nbglmI0@code)  
  
  # modèle binomial négatif, avec un phi global
  nbglmSgb <- negbin(counts~cond, data=snpdata, random=~1, fixpar=list(3,phiGlobal))
  nbglmAgb <- negbin(counts~cond + path, data=snpdata, random=~1, fixpar=list(4,phiGlobal))
  nbglmIgb <- negbin(counts~cond * path, data=snpdata, random=~1, fixpar=list(5,phiGlobal))	
  
  nbAnovgb <- anova(nbglmSgb, nbglmAgb, nbglmIgb)
  nbAICgb <- c(AIC(nbglmSgb,k = log(nbAll))@istats$AIC, AIC(nbglmAgb,k = log(nbAll))@istats$AIC, AIC(nbglmIgb,k = log(nbAll))@istats$AIC)
  # the BIC in fact, since we use k = log(nobs)
  nbSingHesgb <- c(nbglmSgb@singular.hessian, nbglmAgb@singular.hessian, nbglmIgb@singular.hessian)
  nbCodegb <- c(nbglmSgb@code, nbglmAgb@code, nbglmIgb@code)
  	
  # modèle binomial négatif, avec un phi DSS
  nbglmS <- negbin(counts~cond, data=snpdata, random=~1, fixpar=list(3, phiDSS))
  nbglmA <- negbin(counts~cond + path, data=snpdata, random=~1, fixpar=list(4, phiDSS))
  nbglmI <- negbin(counts~cond * path, data=snpdata, random=~1, fixpar=list(5, phiDSS))
  
  nbAnov <- anova(nbglmS, nbglmA, nbglmI)
  nbAIC <- c(AIC(nbglmS,k = log(nbAll))@istats$AIC, AIC(nbglmA,k = log(nbAll))@istats$AIC, AIC(nbglmI,k = log(nbAll))@istats$AIC)
  nbSingHes <- c(nbglmS@singular.hessian, nbglmA@singular.hessian, nbglmI@singular.hessian)
  nbCode <- c(nbglmS@code, nbglmA@code, nbglmI@code)
  
  # modèle binomial négatif, avec un phi DDS, conditionnellement par rapport à la moyenne d'expression  
  nbglmScond <- negbin(counts~cond, data=snpdata, random=~1, fixpar=list(3, phiDSScond))
  nbglmAcond <- negbin(counts~cond + path, data=snpdata, random=~1, fixpar=list(4, phiDSScond))
  nbglmIcond <- negbin(counts~cond * path, data=snpdata, random=~1, fixpar=list(5, phiDSScond))

  nbAnovcond <- anova(nbglmScond, nbglmAcond, nbglmIcond)
  nbAICcond <- c(AIC(nbglmScond,k = log(nbAll))@istats$AIC, AIC(nbglmAcond,k = log(nbAll))@istats$AIC, AIC(nbglmIcond,k = log(nbAll))@istats$AIC)
  nbSingHescond <- c(nbglmScond@singular.hessian, nbglmAcond@singular.hessian, nbglmIcond@singular.hessian) 
  nbCodecond <- c(nbglmScond@code, nbglmAcond@code, nbglmIcond@code) 
  
  rslts <- c(nbAnov0@anova.table$'P(> Chi2)'[2:3],
  nbAnovgb@anova.table$'P(> Chi2)'[2:3],
  nbAnov@anova.table$'P(> Chi2)'[2:3],
  nbAnovcond@anova.table$'P(> Chi2)'[2:3],
             nbAIC0,
             nbAICgb, 
             nbAIC,
             nbAICcond,
             
             nbCode0,
             nbCodegb,
             nbCode,
             nbCodecond,
             
             nbSingHes0,
             nbSingHesgb,
             nbSingHes,
             nbSingHescond)
 ## NB: here useless, colnames re-used after 
  ## names(rslts) <- c("(0)A vs S","(0)I vs A",
  ##                   "(gb)A vs S","(gb)I vs A",
  ##                   "A vs S","I vs A",
  ##                   "(c)A vs S","(c)I vs A",
                    
  ##                   "(0)bicS","(0)bicA","(0)bicI",
  ##                   "(gb)bicS","(gb)bicA","(gb)bicI",
  ##                   "bicS","bicA","bicI",
  ##                   "(c)bicS","(c)bicA","(c)bicI",
                    
  ##                   "(0)codeS","(0)codeA","(0)codeI",
  ##                   "(gb)codeS","(gb)codeA","(gb)codeI", 
  ##                   "codeS","codeA","codeI",
  ##                   "(c)codeS","(c)codeA","(c)codeI",
                    
  ##                   "(0)shS","(0)shA","(0)shI",
  ##                   "(gb)shS","(gb)shA","(gb)shI", 
  ##                   "shS","shA","shI",
  ##                   "(c)shS","(c)shA","(c)shI")
  
  return(rslts)  
}


###################################################
### code chunk number 17: pALLGlobalPhi.glm.nb
###################################################
pALLGlobalPhi.glm.nb=data.frame(t(rep(NA,44)))

for (i in 1:length(allSNPtables)) {
  pALLGlobalPhi.glm.nb[i, ] = try(fitNBglmModelsDSSPhi(allSNPtables[[i]],dispersion(dispData)[i],dispersion(dispDataMeanCond)[i], phi) ,silent=T)
}



###################################################
### code chunk number 18: errors
# ###################################################
# print(pALLGlobalPhi.glm.nb[grepl("Error",pALLGlobalPhi.glm.nb[ ,1]),1])


###################################################
### code chunk number 19: excl_errors
###################################################
nonsing.snps = which(!grepl("Error",pALLGlobalPhi.glm.nb[ ,1]))

#Here I think that we could use the new name of for pALLGlobalPhi.glm.nb  all along. Or better rename it with a clearer name. I think that it is not GlobalPhi only. For now, I just modified it locally, to avoid computing twice the glm for non singular hessian snps. 

pALLGlobalPhi.glm.nb.nonsing=data.frame(t(rep(NA,44)))
j=1
for (i in nonsing.snps) {
  pALLGlobalPhi.glm.nb.nonsing[j, ] = as.numeric(pALLGlobalPhi.glm.nb[i, ])
  j=j+1
}

pALLGlobalPhi.glm.nb=pALLGlobalPhi.glm.nb.nonsing


colnames(pALLGlobalPhi.glm.nb) <- c("(0)A vs S","(0)I vs A",
                                    "(gb)A vs S","(gb)I vs A",
                                    "A vs S","I vs A",
                                    "(c)A vs S","(c)I vs A",
                                    
                                    "(0)bicS","(0)bicA","(0)bicI",
                                    "(gb)bicS","(gb)bicA","(gb)bicI",
                                    "bicS","bicA","bicI",
                                    "(c)bicS","(c)bicA","(c)bicI",
                                    
                                    "(0)codeS","(0)codeA","(0)codeI",
                                    "(gb)codeS","(gb)codeA","(gb)codeI", 
                                    "codeS","codeA","codeI",
                                    "(c)codeS","(c)codeA","(c)codeI",
                                    
                                    "(0)shS","(0)shA","(0)shI",
                                    "(gb)shS","(gb)shA","(gb)shI", 
                                    "shS","shA","shI",
                                    "(c)shS","(c)shA","(c)shI")

rownames(pALLGlobalPhi.glm.nb) <- dataPart2[nonsing.snps,1]



###################################################
### code chunk number 20: repeats
###################################################

tableCount <- function(df){
  df_pen <- df
  tt <- table(df[ ,c(2,4)]) 
  for (i in 1:nrow(tt)) {
		for (j in 1:ncol(tt))
			tt[i,j] <- sum(df[df[ ,2]==rownames(tt)[i]&df[ ,4]==colnames(tt)[j],3])
	}
  return(tt)
}

repeats <- which( apply(pALLGlobalPhi.glm.nb[ ,c(33:34,36:37,39:40,42:43)],1,sum)==8)
# here repeats is the 


###################################################
### code chunk number 21: reapeats_counts
###################################################
# for(i in repeats)
# {
#   print(as.character(dataPart2[i,1]))
#   print(tableCount(allSNPtables[[i]]))
# }
pALLGlobalPhi.glm.nb = pALLGlobalPhi.glm.nb[!is.na(pALLGlobalPhi.glm.nb[ ,1]), ]


###################################################
### code chunk number 22: best model
###################################################
bestmodel.table.n = apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)
bestmodel.table = bestmodel.table.n
bestmodel.table[bestmodel.table==1]="Poisson"
bestmodel.table[bestmodel.table==2]="NB, global phi"
bestmodel.table[bestmodel.table==3]="NB, DSS phi"
bestmodel.table[bestmodel.table==4]="NB, cond DSS phi"
bestmodel.singhes = c()
for (i in 1:length(bestmodel.table.n)) {
  bestmodel.singhes[i] = c(pALLGlobalPhi.glm.nb[i,c(35,38,41,44)])[bestmodel.table.n[i]]
}
bestmodel.singhes = unlist(bestmodel.singhes)
bestmodel = table(bestmodel.table)

bestmodel2 = table(bestmodel.table,bestmodel.singhes)
colnames(bestmodel2) = c("F","T")

bestmodel3 = table(bestmodel.table,apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum),dnn=c("best model","singular hessian"))

colnames(bestmodel3) = paste(colnames(bestmodel3), "Models")


###################################################
### code chunk number 23: print bestmodel
###################################################
# print(bestmodel)


###################################################
### code chunk number 24: print bestmodel2
###################################################
# print(bestmodel2)
# print(bestmodel3)


###################################################
### code chunk number 25: function addOneCount
###################################################
addOneCount <- function(df)
{
  df.pen <- df #? why copy?
  tt <- table(df[ ,c(2,4)]) 
	for (i in 1:nrow(tt)) {
		for (j in 1:ncol(tt))
			tt[i,j] <- sum(df[df[ ,2]==rownames(tt)[i]&df[ ,4]==colnames(tt)[j],3])
	}
	tt = as.data.frame(tt)
  
  for (j in unique(as.character(tt$cond))) {
    
    for (k in unique(as.character(tt$path))) {
      l = sample(which(df$cond==j & df$path==k ),1,replace=F)
      df.pen[l,3] = df.pen[l,3]+1
      
    }
  }
	return(df.pen)
}


###################################################
### code chunk number 26: singhes
###################################################
# singhes = which(bestmodel.table_n >1 & bestmodel.singhes !=0)
# singhes_n = names(singhes) #SNPs pour lesquels le modèle de Poisson n'est pas le plus adapté et pour lesquels l'hessienne est singulière lors de l'ajustement du modèle binomial négatif
# 
# pALLGlobalPhi.glm.nb_pen = pALLGlobalPhi.glm.nb
# for(i in singhes){
#   
#   pALLGlobalPhi.glm.nb_pen[i, ] = try(fitNBglmModelsDSSPhi(addOneCount(allSNPtables[[which(rownames(dataPart2)==singhes_n[i])]]),dispersion(dispData)[which(rownames(dataPart2)==singhes_n[i])],dispersion(dispDataMeanCond)[which(rownames(dataPart2)==singhes_n[i])], phi),silent=T)
# }


###################################################
### code chunk number 27: glmnet
###################################################

pALLGlobalPhi.glm.nb.glmnet = pALLGlobalPhi.glm.nb
pALLGlobalPhi.glm.nb.glmnet$glmnet.pval = 1
pALLGlobalPhi.glm.nb.glmnet$glmnet.code = 0
singhes0 = which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min) == 1)# SNP pour lesquels le modèle poissonien est plus adapté
for (i in singhes0) {
  Xinter   = model.matrix(~cond*path,data= allSNPtables[[which(rownames(dataPart2)==rownames(pALLGlobalPhi.glm.nb)[i])]]); 
  outinter = glmnet(Xinter,allSNPtables[[which(rownames(dataPart2)==rownames(pALLGlobalPhi.glm.nb)[i])]]$counts,family="poisson",lambda=1e-4,alpha=0)
	Xprinc   = model.matrix(~path+cond,data= allSNPtables[[which(rownames(dataPart2)==rownames(pALLGlobalPhi.glm.nb)[i])]]); 
	outprinc = glmnet(Xprinc,allSNPtables[[which(rownames(dataPart2)==rownames(pALLGlobalPhi.glm.nb)[i])]]$counts,family="poisson",lambda=1e-4,alpha=0)
	Pv       = 1-pchisq(deviance(outprinc) - deviance(outinter),df=1)
	pALLGlobalPhi.glm.nb.glmnet$glmnet.pval[i] = Pv
	pALLGlobalPhi.glm.nb.glmnet$glmnet.code[i] = outinter$jerr
}


###################################################
### code chunk number 28: print singhes2
###################################################
# for (i in singhes[apply(pALLGlobalPhi.glm.nb_pen[singhes,33:44],1,sum)!=0]){
#     if(!i%in%repeats){ 
#     print(as.character(dataPart2[i,1]));
#     print(tableCount(allSNPtables[[i]]));
#     print(paste("Best model: ",
#    c("Poisson","NB,global phi","NB, DSS phi","NB, cond DSS phi")[which.min(pALLGlobalPhi.glm.nb_pen[i,c(11,14,17,20)])],sep = "",collapse=""))
#    print("-------------------------")
# }
# }


###################################################
### code chunk number 29: final_pval
###################################################

pALLGlobalPhi.glm.nb$final.pval.a.ia = 1

i <- 1 #modèle de Poisson
  
#	li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)==0)
#	pALLGlobalPhi.glm.nb$final.pval.a.ia[li] = pALLGlobalPhi.glm.nb[li,2]
	
	li.singhes <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i)
	pALLGlobalPhi.glm.nb$final.pval.a.ia[li.singhes] = pALLGlobalPhi.glm.nb.glmnet$glmnet.pval[li.singhes]
	
i  <- 2 # modèle binomial négatif, avec un phi global
	
	li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)==0)
	pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li,4]
	
# 	li.singhes = which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)!=0)
# 	pALLGlobalPhi.glm.nb$final.pval.a.ia[li.singhes] = pALLGlobalPhi.glm.nb_pen[li.singhes,4]

i  <- 3 # modèle binomial négatif, avec les phi estimés avec DSS
	
	li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)==0)
	pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li,6]
	
# 	li.singhes <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)!=0)
# 	pALLGlobalPhi.glm.nb$final.pval.a.ia[li.singhes] = pALLGlobalPhi.glm.nb_pen[li.singhes,6]

i  <- 4 # modèle binomial négatif, avec les phi estimés avec DSS, conditionnellement à la moyenne de l'expression
	
	li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)==0)
	pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li,8]
	
# 	li.singhes <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)!=0)
# 	pALLGlobalPhi.glm.nb$final.pval.a.ia[li.singhes] <- pALLGlobalPhi.glm.nb_pen[li.singhes,8]
	
pALLGlobalPhi.glm.nb$final.padj.a.ia <- p.adjust(pALLGlobalPhi.glm.nb$final.pval.a.ia, method="fdr")

# print(table(pALLGlobalPhi.glm.nb$final.padj.a.ia <= 0.05,p.adjust(pALLGlobalPhi.glm.nb[ ,4],method="fdr")<=0.05, dnn = c("Analysis2 - New analysis","Analysis1 - NegBin, global phi")))	






##################### Modified Vincent & Alice 29 / 11/ 13
signifEvents <- cbind(dataPart2[nonsing.snps, ],pALLGlobalPhi.glm.nb$final.padj.a.ia )[ pALLGlobalPhi.glm.nb$final.padj.a.ia <= 0.05, ]
######### sorting by deltaPSI / deltaF
finalDelta <- NA 
for ( i in 1:(n-1)) {#n =  nb conditions
# 1 er condition
  PSI.replicat1 <- c()
  for (numReplicatI in 1:nr[i]) {          
      nameUp <- paste('UP_Cond', i, '_R', numReplicatI, '_Norm', sep='')
      nameLow <- paste('LP_Cond', i, '_R', numReplicatI, '_Norm', sep='') 
      tmp.psi <- signifEvents[ , nameUp] / (signifEvents[ ,nameUp] + signifEvents[ ,nameLow])
      PSI.replicat1 <- cbind( PSI.replicat1, tmp.psi)
    }
  PSIcondI <- apply(PSI.replicat1, MARGIN=1, mean, na.rm = T)
  for (j in (i+1):n) {     #j = 2nd condition      
      PSI.replicat2 <- c()
      for (numReplicatI in 1:nr[j]) {          
          nameUp <- paste('UP_Cond', j, '_R', numReplicatI, '_Norm', sep='')
          nameLow <- paste('LP_Cond', j, '_R', numReplicatI, '_Norm', sep='') 
          tmp.psi <- signifEvents[ ,nameUp] / (signifEvents[ ,nameUp] + signifEvents[ ,nameLow])
          PSI.replicat2 <- cbind( PSI.replicat2, tmp.psi)
        }
      PSIcondJ <- apply(PSI.replicat2, MARGIN=1, mean, na.rm = T)
      
      
      # deltaPSI for cond i,j
      deltaPSIij <- abs( PSIcondJ- PSIcondI ) 
      
      finalDelta <- apply( cbind( finalDelta, deltaPSIij) , MARGIN = 1, max, na.rm = T)
    }
}

signifEvents <- cbind(signifEvents, finalDelta)# adding DeltaPsi/f to the final table
colnames(signifEvents)[length(colnames(signifEvents))] <- 'Deltaf/DeltaPSI'# renaming last columns
signifEvents.sorted <- signifEvents[ order( finalDelta, decreasing = T), ]

########################################## Writing final file  with significant events #########################################################################



# write.table( signifEvents_sorted , paste(dataName, "Analysis.txt",sep="",collapse="" ), sep="\t", quote=F, row.names=F, col.names=T )

####################################################################
##################### End modified  Vincent & Alice 29 / 11/ 13


###################################################
### code chunk number 30: plotVariance_vs_Mean_SNP2
###################################################
if (option.Plot == TRUE) {
plot( snp.mean.variance.df$Mean, snp.mean.variance.df$Variance, 
     xlab="Mean SNP count", 
     ylab="Variance SNP count",
     log="xy", las=1)
points( snp.mean.variance.df$Mean[pALLGlobalPhi.glm.nb$final.padj.a.ia <= 0.05 & p.adjust(pALLGlobalPhi.glm.nb[ ,4],method="fdr")<=0.05], snp.mean.variance.df$Variance[pALLGlobalPhi.glm.nb$final.padj.a.ia <= 0.05 & p.adjust(pALLGlobalPhi.glm.nb[ ,4],method="fdr")<=0.05], col="green")     

points(snp.mean.variance.df$Mean[pALLGlobalPhi.glm.nb$final.padj.a.ia <= 0.05 & p.adjust(pALLGlobalPhi.glm.nb[ ,4],method="fdr")>0.05], snp.mean.variance.df$Variance[pALLGlobalPhi.glm.nb$final.padj.a.ia <= 0.05 & p.adjust(pALLGlobalPhi.glm.nb[ ,4],method="fdr")>0.05], col="blue")     

points(snp.mean.variance.df$Mean[pALLGlobalPhi.glm.nb$final.padj.a.ia > 0.05 & p.adjust(pALLGlobalPhi.glm.nb[ ,4],method="fdr")<=0.05], snp.mean.variance.df$Variance[pALLGlobalPhi.glm.nb$final.padj.a.ia > 0.05 & p.adjust(pALLGlobalPhi.glm.nb[ ,4],method="fdr")<=0.05], col="red")     

abline(a=0, b=1, col=2, lwd=2)
lines(x,yQP,col=3, lwd=2)
lines(x,yNB,col=6, lwd=2)
legend("topleft", c("Poisson","Quasi-Poisson", "Negative Binomial"), text.col=c(2,3,6), box.lty=0)
legend("bottomright", c("Signif Analysis1 & Signif Analysis2","Signif Analysis1 & No Signif Analysis2", "No Signif Analysis1 & Signif Analysis2"), text.col=c("green","red","blue"), box.lty=0)  
}



###################################################
### code chunk number 31: histogram_NB_a_vs_ia
###################################################
if (option.Plot == TRUE) {
hist(pALLGlobalPhi.glm.nb$final.pval.a.ia, breaks=20)
}

###################################################
### code chunk number 32: sessioninfo
###################################################
# sessionInfo()

return(signifEvents.sorted)
}


