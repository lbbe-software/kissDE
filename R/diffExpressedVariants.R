.lineParse <- function(line, indexStart, isQuality) {
  options(warn=-1)
  beginningLineToWrite <- ""

  splitElements <- strsplit(line, "|", fixed=TRUE)[[1]] # splits the line
  if (indexStart == 6) {
    for (k in 1:(indexStart-2)) {
    beginningLineToWrite <- paste (beginningLineToWrite, splitElements[k], sep="|") #writes the firsts elements of the line : bcc, cycle... but NOT branching_nodes
    }
  } else {
    for (k in 1:(indexStart-1)) {
    beginningLineToWrite <- paste (beginningLineToWrite, splitElements[k], sep="|") #writes the firsts elements of the line : bcc, cycle... 
    }
  }

  ElementsNb <- length(splitElements) #number of elements in the line
  splitCounts <- splitElements[indexStart : (length(splitElements) - 1)] # avoids the name of the bcc ... and the rank (last one) to get only the counts
  s <- sapply(splitCounts, function(splitCounts) regmatches(splitCounts[[1]], gregexpr(pattern = "[0-9]+",splitCounts[[1]]))) #gets the junctions id (ex 1 in AS1) and the count (ex 6 in AS1_6)
  

  if (isQuality == TRUE ) { #if there is a quality information  (for SNPs)
    s<- s[regexpr('Q',names(s)) == -1] # we discard "Q_" information as they are not counts
  } 
  return(list(beginning=beginningLineToWrite,countsperCond=s))
}

.countsSet <- function(line, indexStart, counts=0, pairedEnd=FALSE, order=NULL, exonicReads=TRUE, isQuality) {
  resultParsing <- .lineParse(line, indexStart, isQuality)
  beginningLineInfo <- resultParsing$beginning
  countsperCond<- resultParsing$countsperCond
  nbVec <- rep(0, length(countsperCond))
  countsVec <- rep(0, length(countsperCond))
  psiVec <- rep(0, length(countsperCond))
  for (i in 1:length(countsperCond)) {
    nbVec[i] <- as.numeric(countsperCond[[i]][1])
    countsVec[i] <- as.numeric(countsperCond[[i]][2])
    if (counts > 1) { #specific issues linked with --counts option
      if (grepl("ASSB", names(countsperCond)[i]) == TRUE) { #so that counts on ASSB junction are not counted twice.
        psiVec[i] <- countsVec[i]
        countsVec[i] <- - countsVec[i]
      }
      if ((counts == 2) & (exonicReads == FALSE)) {
        if (grepl("^S[0-9]+", names(countsperCond)[i]) == TRUE) { #when exonic reads are not wanted we must discard reads counted in S_X
          countsVec[i] <- 0
        }
      }
    }
  }
  if (counts > 1) {
    d <- data.frame(nbVec,countsVec)
    names(d) <- c("NB", "COUNTS")
    sums <- aggregate(d$COUNTS, by=list(d$NB), sum) #sums the counts for each junction that belongs to the same event
    dpsi <- data.frame(nbVec,psiVec)#dpsi will store counts of ASSB counts 
    names(dpsi) <- c("NB", "ASSB")
    assbPsi <- aggregate(dpsi$ASSB, by=list(dpsi$NB), sum)
    if (pairedEnd == TRUE) {
      if (is.null(order)) {
        order <- rep(1:((dim(sums)[1])/2), rep(2,((dim(sums)[1])/2)))
      } else {
        if (! is.vector(order)) {
          print("Error, order vector seems to be in a wrong format.")
        }
      }
      d2 <- data.frame(order, sums)
      names(d2)[3] <- 'sums'
      sums2 <- aggregate(d2$sums, by=list(d2$order), sum) # in case data is paired-end, there is one more sum to do, for each part of the pair
      sums <- sums2
      dpsi2 <- data.frame(order, assbPsi)
      names(dpsi2) [3] <- 'sums'
      assbPsi2 <- aggregate(dpsi2$sums, by=list(dpsi2$order), sum)
      assbPsi <- assbPsi2
    } 
    listASSB <- t(assbPsi)[2,]
  } else { ### counts == 0 
    if (pairedEnd == TRUE) {
      if (is.null(order)) {
        order <- rep(1:(length(countsperCond)/2), rep(2,(length(countsperCond)/2))) # for length(s)=8, will create a vector c(1,1,2,2,3,3,4,4) (assuming data is ordered)
      } else {
        if (! is.vector(order)) {
          print("Error, order vector seems to be in a wrong format.")
        }
      }
    } else {
        order <-c(1:length(countsperCond))
    }
    d <- data.frame(order,countsVec)
    names(d) <- c("ORDER", "COUNTS")
    sums <- aggregate(d$COUNTS, by=list(d$ORDER), sum)
    listASSB <- NULL
  }
    listCounts <- t(sums)[2,]
  return(list(firstPart=beginningLineInfo, vCounts=listCounts, psiCounts= listASSB))
}

.addOneCount <- function(df)
{
  df$counts <- unlist(lapply(df[,'counts'],function(x){x+1}))
  return (df)
}

.getInfoLine <- function(line, counts=0, pairedEnd=FALSE, order=NULL, exonicReads=TRUE, isQuality) {
  if ( grepl("branching_nodes", line) ) { 
    indexStart <- 6 
    } else {
    indexStart <- 5
  }
  resultCountsSet <- .countsSet(line, indexStart, counts, pairedEnd, order, exonicReads, isQuality)
  lineFirstPart <- resultCountsSet$firstPart
  lineFirstPartSplit <- strsplit(lineFirstPart,"|",fixed=TRUE)[[1]]
  name <- paste(lineFirstPartSplit[2],lineFirstPartSplit[3],sep="|")
  name <- substr(name, start = 2, stop = nchar(name))
  length <- strsplit(lineFirstPartSplit[5],"_")[[1]][4]
  vCounts <- resultCountsSet$vCounts
  return (list(eventName=name,variantLength=length,variantCounts=vCounts, psiInfo=resultCountsSet$psiCounts))
}


.readAndPrepareData <- function(countsData,conditions) {
  ###################################################
  ### code chunk number 1: Read data
  ###################################################
  if (is.null(countsData$psiInfo)){
    countsEvents <- countsData #count table provided by the user
    psiInfo <- NULL
  } else {
    countsEvents <- countsData$countsEvents #provided by kissplice2counts
    if ( dim(countsData$psiInfo)[2] >1 ){
      psiInfo <- countsData$psiInfo#info about ASSB counts
    } else {
      psiInfo <- NULL
    }
  }
  sortedconditions <- sort(conditions)
  n <- length(unique(sortedconditions))
  nr <- rle(sortedconditions)$lengths
  sortedindex <- order(conditions)+2
  namesData <- c("ID","Length",rep(NA,length(conditions)))
  for (k in 1:nr[1]){
    namesData[2+k] <- paste(sortedconditions[k],"_r",k,sep="",collapse="")
  }
  for (i in 2:n) {
    for (j in 1:nr[n]) {
      namesData[2+cumsum(nr)[i-1]+j] <- paste(sortedconditions[cumsum(nr)[i-1]+j],"_r",j,sep="",collapse="")
    }
  }#proper names for conditionsXrelicates
  countsEvents[,-(1:2)] <- countsEvents[,sortedindex]
  colnames(countsEvents) <- namesData
  if(! is.null(psiInfo)){
    psiInfo[,-1] <- psiInfo[,sortedindex-1]
    colnames(psiInfo) <- c("events.names",namesData[c(-1,-2)])
    ASSBinfo <- data.frame(psiInfo)
    colnames(ASSBinfo) <- c("events.names",namesData[c(-1,-2)])
    colnames(psiInfo) <- namesData[c(-1,-2)]
  } else {
    ASSBinfo <- NULL
  }
  countsEvents$Path <- gl( 2, 1, dim(countsEvents)[1], labels = c("UP", "LP"))

  ###################################################
  ### code chunk number 2: Normalisation
  ###################################################
  # Normalisation with DESeq
  conds <- c()
  for( i in 1:n ) {
    for( j in 1:nr[i] ) {
      conds <- c( conds,paste( "Cond", i, sep = "",collapse = "") )
    }
  } 
  cds <- newCountDataSet( countsEvents[ ,3:(3+length(conds)-1)], conds ) # create object
  cdsSF <- estimateSizeFactors(cds)
  sizeFactors( cdsSF )
  shouldWeNormalise=sum(is.na(sizeFactors(cdsSF))) < 1
  dim <- dim(countsEvents)[2]
  countsEvents[ ,(dim+1):(dim+length(conds)) ] <- round(counts(cdsSF, normalized=shouldWeNormalise))
  colnames(countsEvents)[(dim+1):(dim+length(conds))] <- paste(namesData[3:(3+sum(nr)-1)],"_Norm",sep="")
  return(list(countsEvents, conds, dim, n, nr, sortedconditions,ASSBinfo=ASSBinfo))
}

.eventtable <- function(df,startPosColumn4Counts, endPosCol4Counts){
  eventTab = data.frame(ID=rep(as.factor(df['ID']), endPosCol4Counts-startPosColumn4Counts+1),
  cond=as.factor(unlist(lapply(strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts],"_"),FUN=function(d){d[2]}))),
  counts=as.numeric(df[startPosColumn4Counts:endPosCol4Counts]),
  path=as.factor(unlist(lapply(strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts],"|"),FUN=function(d){d[1]}))),
  row.names = NULL)
  return(eventTab)
}

.fitNBglmModelsDSSPhi <- function(eventdata, phiDSS, phiDSScond, phiGlobal, nbAll){
    
  nbglmA0 <- negbin(counts~cond + path, data=eventdata, random=~1, fixpar=list(4,0))
  nbglmI0 <- negbin(counts~cond * path, data=eventdata, random=~1, fixpar=list(5,0))  
    # S: simple, A: additive, I : interaction models
    # Poisson model  
  nbAnov0 <- anova(nbglmA0, nbglmI0)
  nbAIC0 <- c(AIC(nbglmA0,k = log(nbAll))@istats$AIC, AIC(nbglmI0,k = log(nbAll))@istats$AIC)
    # singular.hessian:  true when fitting provided a singular hessian, indicating an overparamaterized model.
  nbSingHes0 <- c(nbglmA0@singular.hessian, nbglmI0@singular.hessian)
    # code: ‘code’ An integer (returned by ‘optim’) indicating why the optimization process terminated.
  nbCode0 <- c(nbglmA0@code, nbglmI0@code)  
  
    # binomial negative model, with global phi
  nbglmAgb <- negbin(counts~cond + path, data=eventdata, random=~1, fixpar=list(4,phiGlobal))
  nbglmIgb <- negbin(counts~cond * path, data=eventdata, random=~1, fixpar=list(5,phiGlobal)) 
  
  nbAnovgb <- anova( nbglmAgb, nbglmIgb)
  nbAICgb <- c(AIC(nbglmAgb,k = log(nbAll))@istats$AIC, AIC(nbglmIgb,k = log(nbAll))@istats$AIC)
    # the BIC in fact, since we use k = log(nobs)
  nbSingHesgb <- c( nbglmAgb@singular.hessian, nbglmIgb@singular.hessian)
  nbCodegb <- c(nbglmAgb@code, nbglmIgb@code)
    
    # binomial negative model, with phi DSS
  nbglmA <- negbin(counts~cond + path, data=eventdata, random=~1, fixpar=list(4, phiDSS))
  nbglmI <- negbin(counts~cond * path, data=eventdata, random=~1, fixpar=list(5, phiDSS))
  
  nbAnov <- anova( nbglmA, nbglmI)
  nbAIC <- c( AIC(nbglmA,k = log(nbAll))@istats$AIC, AIC(nbglmI,k = log(nbAll))@istats$AIC)
  nbSingHes <- c( nbglmA@singular.hessian, nbglmI@singular.hessian)
  nbCode <- c( nbglmA@code, nbglmI@code)
  
    # binomial negative model, with phi DSS, conditionally to the expression mean  
  nbglmAcond <- negbin(counts~cond + path, data=eventdata, random=~1, fixpar=list(4, phiDSScond))
  nbglmIcond <- negbin(counts~cond * path, data=eventdata, random=~1, fixpar=list(5, phiDSScond))

  nbAnovcond <- anova( nbglmAcond, nbglmIcond)
  nbAICcond <- c( AIC(nbglmAcond,k = log(nbAll))@istats$AIC, AIC(nbglmIcond,k = log(nbAll))@istats$AIC)
  nbSingHescond <- c(nbglmAcond@singular.hessian, nbglmIcond@singular.hessian) 
  nbCodecond <- c(nbglmAcond@code, nbglmIcond@code) 
  
  rslts <- c(nbAnov0@anova.table$'P(> Chi2)'[2],
  nbAnovgb@anova.table$'P(> Chi2)'[2],
  nbAnov@anova.table$'P(> Chi2)'[2],
  nbAnovcond@anova.table$'P(> Chi2)'[2],
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
  return(rslts)  
}

kissplice2counts <- function(fileName, counts=0, pairedEnd=FALSE, order=NULL, exonicReads =TRUE) {
  toConvert <- file(fileName, open = "r")
  lines <- readLines(toConvert)
  line <- lines[1]
  isQuality <- grepl("Q",line[1])
  resultLine1 <- .getInfoLine(line, counts, pairedEnd, order, exonicReads, isQuality)#get all the informations for the 1st line
  eventName <- resultLine1$eventName
  variantLength <- resultLine1$variantLength
  variantCounts <- resultLine1$variantCounts
  events.mat <- matrix(NA,length(lines)/2,length(variantCounts)+1)
  events.names <- rep(NA,length(lines)/2)
  events.mat[1,1] <- as.numeric(variantLength)
  events.mat[1,2:dim(events.mat)[2]] <- variantCounts
  events.names[1] <- eventName
  index <- 3
  indexNames <- 2
  firstLineChar <- substr(lines[index], start = 0, stop = 1)
  psiInfo <- matrix(NA,length(lines)/2,length(resultLine1$psiInfo))
  psiInfo[1,] <- resultLine1$psiInfo
  if (firstLineChar == '>') { #same for all other lines, ignore lines with sequences
    while (index <= length(lines)) {
      line <- lines[index]
      resultLine <- .getInfoLine(line, counts, pairedEnd, order, exonicReads, isQuality)
      eventName <- resultLine$eventName
      variantLength <- resultLine$variantLength
      variantCounts <- resultLine$variantCounts
      events.mat[indexNames,1] <- as.numeric(variantLength)
      events.mat[indexNames,2:dim(events.mat)[2]] <- variantCounts
      events.names[indexNames] <- eventName
      psiInfo[indexNames,] <- resultLine$psiInfo
      index <- index + 2
      indexNames <- indexNames + 1
    }
  }
  class(events.mat) <- "numeric"
  events.df <- as.data.frame(events.mat)
  events.df <- data.frame(events.names,events.df)
  close(toConvert)
  psidf <- as.data.frame(psiInfo)
  psiInfo.df <-  data.frame(events.names,psidf)
  return (list(countsEvents=events.df,psiInfo=psiInfo.df))
}

qualityControl <- function(countsData,conditions,storeFigs=FALSE, pathFigs="None") {

  options(warn=-1) # suppress the warning for the users

  if (storeFigs == TRUE){
    if (pathFigs == "None") {
      pathToFigs = "kissDEFigures"
    } else {
      pathToFigs = paste(pathFigs,"/kissDEFigures",sep="")
    }
    find = paste("find",pathToFigs)
    d<-system(find,TRUE,ignore.stderr=TRUE)
    if (length(d) == 0) { 
      command = paste("mkdir", pathToFigs)
      system(command,ignore.stderr=TRUE)
    }
  }
  ###################################################
  ### code chunk number 1: Read and prepare data
  ###################################################
  listData <-.readAndPrepareData(countsData,conditions)
  countsData <- listData[[1]]
  conds <- listData[[2]]
  dim <- listData[[3]]
  n <- listData[[4]]
  nr <- listData[[5]]

  ###################################################
  ### code chunk number 2: dendrogram
  ###################################################
  if (storeFigs == FALSE) {
    plot(hclust(as.dist(1-cor(countsData[ ,(dim+1):(dim+length(conds))])),"ward"))
    par(ask=TRUE)
    } else {
        filename = paste(pathToFigs,"/dendrogram.png",sep="")
        png(filename)
        plot(hclust(as.dist(1-cor(countsData[ ,(dim+1):(dim+length(conds))])),"ward"))
        void <- dev.off()
    }

  ###################################################
  ### code chunk number 3: replicates
  ###################################################
  if (storeFigs == FALSE) {
    heatmap(as.matrix(as.dist(1-cor(countsData[ ,(dim+1):(dim+length(conds))]))), margins = c(10,10))
    } else {
        filename = paste(pathToFigs,"/heatmap.png",sep="")
        png(filename)
        heatmap(as.matrix(as.dist(1-cor(countsData[ ,(dim+1):(dim+length(conds))]))), margins = c(10,10))
        void <- dev.off()
    }

  ###################################################
  ### code chunk number 4: intra-group and inter-group-variance
  ###################################################
  ## Mean and variance over all conditions and replicates (normalized counts!) 
  countsData$mn <- apply(countsData[ ,(dim+1):(dim+length(conds))], 1, mean)
  countsData$var <- apply(countsData[ ,(dim+1):(dim+length(conds))], 1, var)
  ## correction term
  nbAll <- sum(nr) # number of all observations in all groups
  countsData$ct  <- apply(countsData[ ,(dim+1):(dim+length(conds))], 1, sum)^2/nbAll
  ## sum of squares between groups
  countsData$ss <- apply(countsData[ ,(dim+1):(dim+length(conds)/n)],1,sum)^2/nr[1] + apply(countsData[ ,((dim+1)+length(conds)/n):(dim+length(conds))],1,sum)^2/nr[2] 
  ## substract the correction term from the SS and divide by the degrees of 
  df <- 1 # freedom(groups); here: df=2-1=1
  countsData$varInter <- (countsData$ss - countsData$ct)/df
  # intra-variability 
  countsData$varC1 <- apply( countsData[ ,(dim+1):(dim+nr[1])], 1, var )
  countsData$varC2 <- apply( countsData[ ,((dim+1)+nr[1]):(dim+nr[2]+nr[1])], 1, var )
  countsData$varIntra <- apply( data.frame(countsData$varC1, countsData$varC2), 1, mean )

  ###################################################
  ### code chunk number 5: intra-vs-inter
  ###################################################
  if (storeFigs == FALSE) {
      plot( x = countsData$varIntra, y = countsData$varInter, xlab = "Intra-variability", ylab = "Inter-variability", las = 1, log = "xy")
      abline( a = 0, b = 1, col = 2, lty = 2, lwd = 2 )
    } else {
        filename = paste(pathToFigs,"/InterIntraVariability.png",sep="")
        png(filename)
        plot( x = countsData$varIntra, y = countsData$varInter, xlab = "Intra-variability", ylab = "Inter-variability", las = 1, log = "xy")
        abline( a = 0, b = 1, col = 2, lty = 2, lwd = 2 )
        void <- dev.off()
    }
}

diffExpressedVariants <- function(countsData, conditions, storeFigs=FALSE, pathFigs="None", pvalue=0.05, filterLowCountsVariants=10, flagLowCountsConditions=10, readLength=75, overlap=42) {

  options(warn=-1) # suppress the warning for the users

  if (storeFigs == TRUE){
    if (pathFigs == "None") {
      pathToFigs = "kissDEFigures"
    } else {
      pathToFigs = paste(pathFigs,"/kissDEFigures",sep="")
    }
    find = paste("find",pathToFigs)
    d<-system(find,TRUE,ignore.stderr=TRUE)
    if (length(d) == 0) { 
      command = paste("mkdir", pathToFigs)
      system(command,ignore.stderr=TRUE)
    }
  }
  ###################################################
  ### code chunk number 1: Read and prepare data
  ###################################################
  listData <-.readAndPrepareData(countsData,conditions)#makes proper names for colums, normalizes data
  countsData <- listData[[1]]
  n <- listData[[4]]
  nr <- listData[[5]]
  sortedconditions <- listData[[6]]
  ASSBinfo <-  listData$ASSBinfo #in case counts option in kissplice2counts is at 1 or 2, we have info about junction counts (ASSB), that will be useful to correct the computation of delta psi in the end. They are stored here.
  if (! is.null(ASSBinfo)) {
    li<-c()
    for (i in (1:dim(ASSBinfo)[1])){
      if (i%%2!=0) {
        li<-c(li,i)
      }
    }
    ASSBinfo <- ASSBinfo[li,]
  }

  ##################################################
  ## code chunk number 2: event-list
  ##################################################
  # reduce data frame to the interesting columns
  nbAll <- sum(nr)
  dataPart <- countsData[ ,c(1:2,which(grepl("_Norm",names(countsData)))) ] 
  dataPart$Path <- gl( 2,1,dim(countsData)[1],labels=c("UP","LP") ) 
  dataPart2 <- cbind( dataPart[ seq( 1, dim( dataPart )[1], 2 ), ], dataPart[ seq( 2, dim(dataPart)[1] , 2 ), grepl("Norm", names( dataPart) ) ] )
  names(dataPart2)[3:(3+nbAll-1)] <- paste( "UP",names(dataPart2)[3:(3+nbAll-1)], sep="_" )
  names(dataPart2)[(3+nbAll+1):(3+2*nbAll+1-1)] <- paste( "LP", names(dataPart2)[(3+nbAll+1):(3+2*nbAll+1-1)], sep="_" )
  lengths <- data.frame(dataPart[seq(1,dim(dataPart)[1],2),2],dataPart[seq(2,dim(dataPart)[1],2),2])
  colnames(lengths) <- c("upper","lower")
  dataPart2[2] <- lengths$upper -  lengths$lower # computes the difference of length between the lower and upper paths 
  names(dataPart2)[2] <- "Length_diff"
  dataPart2 <- dataPart2[ ,c(-(3+nbAll))]
  if (anyDuplicated(dataPart2[ ,1]) > 0) {
    dataPart2 <- dataPart2[!duplicated(as.character(dataPart2[ ,1])), ]
  }
  rownames(dataPart2)=as.character(dataPart2[ ,1])
  if (! is.null(ASSBinfo)) {
    rownames(ASSBinfo)<-dataPart2[,1]
  }
  rownames(lengths) <-rownames(dataPart2) 
  # create list for the complete data set
  allEventtables  <- apply(dataPart2,1,.eventtable, startPosColumn4Counts = which(grepl("UP",names(dataPart2)))[1],endPosCol4Counts = ncol(dataPart2))

  ###################################################
  ### code chunk number 3: DSS dispersion estimation
  ###################################################
  dataNormCountsEvent <- as.matrix(dataPart2[ ,3:ncol(dataPart2)]) # the counts matrix
  colnames(dataNormCountsEvent) <- 1:ncol(dataNormCountsEvent)
  designs <- rep(c(1:(n*2)),c(nr,nr)) # the design matrix
  dispData <- newSeqCountSet(dataNormCountsEvent, as.data.frame(designs))
  dispData <-estDispersion(dispData)
  dispDataMeanCond <- newSeqCountSet(dataNormCountsEvent, as.data.frame(designs))
  # if (n > 2) {
  #     dispDataMeanCond <- estDispersion(dispData)
  #     } else {
  #       dispDataMeanCond <- estDispersion(dispData,trend=T)
  #     }
  dispDataMeanCond <- estDispersion(dispData)

  ###################################################
  ### code chunk number 4: variance - mean - Event level1
  ###################################################
  # compute mean and variance per Event (instead of per allele)
  event.mean.variance.df <- as.data.frame(cbind(apply(dataPart2[ ,which(grepl("_Norm",names(dataPart2)))],1,mean),apply(dataPart2[ ,which(grepl("_Norm",names(dataPart2)))],1,var)))
  names(event.mean.variance.df) <- c("Mean", "Variance")
  rownames(event.mean.variance.df) <- as.character(dataPart2[ ,1])
  # estimate the dispersion parameter D of the Quasi-Poisson distribution
  lm.D <- lm(event.mean.variance.df$Variance ~ event.mean.variance.df$Mean-1)
  # estimate the overdispersion parameter theta of the NB distritution
  modelNB <- Variance ~ Mean + 1/theta * Mean^2
  nls.modelNB <- nls(modelNB, data=event.mean.variance.df, start=list(theta=100))
  phi <- 1/coef(nls.modelNB) # to be used as fixed parameter later on

  ###################################################
  ### code chunk number 5: plot models
  ###################################################

  #compute model fit
  # log plot, so exclude 0 from the x values
  x <- c(seq(0.1,1,0.1),seq(2,5000,1))
  yQP <- x*coef(lm.D)
  yNB <- x + 1/coef(nls.modelNB) * x^2

  if (storeFigs == FALSE) {
    plot(event.mean.variance.df$Mean, event.mean.variance.df$Variance, 
    xlab="Mean Event count", 
    ylab="Variance Event count",
    log="xy", las=1)
    abline(a=0, b=1, col=2, lwd=2)
    lines(x,yQP,col=3, lwd=2)
    lines(x,yNB,col=6, lwd=2)
    legend("topleft", c("Poisson","Quasi-Poisson", "Negative Binomial"), 
    text.col=c(2,3,6), box.lty=0);
    } else {
        filename = paste(pathToFigs,"/models.png",sep="")
        png(filename)
        plot(event.mean.variance.df$Mean, event.mean.variance.df$Variance, 
        xlab="Mean Event count", 
        ylab="Variance Event count",
        log="xy", las=1)
        abline(a=0, b=1, col=2, lwd=2)
        lines(x,yQP,col=3, lwd=2)
        lines(x,yNB,col=6, lwd=2)
        legend("topleft", c("Poisson","Quasi-Poisson", "Negative Binomial"), 
        text.col=c(2,3,6), box.lty=0);
        void <- dev.off()
    }

  ###################################################
  ### code chunk number 6: exclude low counts
  ###################################################
  totLOW <- as.vector(apply(dataPart2[ ,(3 + sum(nr)):(3 + 2 * sum(nr) - 1)],1,sum)) #global counts for each variant (low/up) by event
  totUP <- as.vector(apply(dataPart2[ ,3:(3 + sum(nr) - 1)],1,sum))

  newindex <- dataPart2[-which(totUP <filterLowCountsVariants & totLOW<filterLowCountsVariants),1]#we filter out variants which counts do not reach the fixed limit
  dataPart3 <- dataPart2[newindex,]
  exprs(dispData) <- exprs(dispData)[-which(totUP<filterLowCountsVariants & totLOW<filterLowCountsVariants),]
  exprs(dispDataMeanCond) <- exprs(dispDataMeanCond)[-which(totUP<filterLowCountsVariants & totLOW<filterLowCountsVariants),]
  allEventtables  <- apply(dataPart3,1,.eventtable, startPosColumn4Counts = which(grepl("UP",names(dataPart3)))[1],endPosCol4Counts = ncol(dataPart3))

  ###################################################
  ### code chunk number 7: pALLGlobalPhi.glm.nb
  ###################################################
  pALLGlobalPhi.glm.nb=data.frame(t(rep(NA,28)))
  for (i in 1:length(allEventtables)) {
    pALLGlobalPhi.glm.nb[i, ] = try(.fitNBglmModelsDSSPhi(allEventtables[[i]],dispersion(dispData)[i],dispersion(dispDataMeanCond)[i], phi, nbAll) ,silent=T)
  }
  ###################################################
  ### code chunk number 8: excl_errors
  ###################################################
  sing.events <- which(grepl("Error",pALLGlobalPhi.glm.nb[ , 1]))
  if (length(sing.events) != 0) {
      pALLGlobalPhi.glm.nb <- pALLGlobalPhi.glm.nb[ - sing.events, ]
  }
  colnames(pALLGlobalPhi.glm.nb) <- c("(0)I vs A",
                                    "(gb)I vs A",
                                    "I vs A",
                                    "(c)I vs A",
                                    
                                    "(0)bicA","(0)bicI",
                                    "(gb)bicA","(gb)bicI",
                                    "bicA","bicI",
                                    "(c)bicA","(c)bicI",
                                    
                                    "(0)codeA","(0)codeI",
                                    "(gb)codeA","(gb)codeI", 
                                    "codeA","codeI",
                                    "(c)codeA","(c)codeI",
                                    
                                   "(0)shA","(0)shI",
                                    "(gb)shA","(gb)shI", 
                                    "shA","shI",
                                    "(c)shA","(c)shI")
  if (length(sing.events) != 0) {
    rownames(pALLGlobalPhi.glm.nb) <- dataPart3[ - sing.events, 1]
  } else {
    rownames(pALLGlobalPhi.glm.nb) <- dataPart3[ , 1]
  }
  pALLGlobalPhi.glm.nb = pALLGlobalPhi.glm.nb[!is.na(pALLGlobalPhi.glm.nb[ , 1]), ]
  matrixpALLGlobalPhi <- as.matrix(pALLGlobalPhi.glm.nb)
  storage.mode(matrixpALLGlobalPhi) <- 'numeric'

  ###################################################
  ### code chunk number 9: best model
  ###################################################
  bestmodel.table.n = apply(matrixpALLGlobalPhi[ ,c(6,8,10,12)],1,which.min)#######"%%%%%"
  bestmodel.table = bestmodel.table.n
  bestmodel.table[bestmodel.table == 1] = "Poisson"
  bestmodel.table[bestmodel.table == 2] = "NB, global phi"
  bestmodel.table[bestmodel.table == 3] = "NB, DSS phi"
  bestmodel.table[bestmodel.table == 4] = "NB, cond DSS phi"
  bestmodel.singhes = c()
  for (i in 1:length(bestmodel.table.n)) {
    bestmodel.singhes[i] = c(matrixpALLGlobalPhi[i,c(22,24,26,28)])[bestmodel.table.n[i]]##########%%%%%%%%%%
  }
  bestmodel.singhes = unlist(bestmodel.singhes)
  bestmodel = table(bestmodel.table)
  bestmodel2 = table(bestmodel.table,bestmodel.singhes)
  colnames(bestmodel2) = c("F","T")
  bestmodel3 = table(bestmodel.table,apply(matrixpALLGlobalPhi[ ,c(22,24,26,28)],1,sum),dnn=c("best model","singular hessian"))
  colnames(bestmodel3) = paste(colnames(bestmodel3), "Models")

  ###################################################
  ### code chunk number 10: glmnet
  ###################################################
  pALLGlobalPhi.glm.nb.glmnet = as.data.frame(matrixpALLGlobalPhi)
  pALLGlobalPhi.glm.nb.glmnet$glmnet.pval = 1
  pALLGlobalPhi.glm.nb.glmnet$glmnet.code = 0
  singhes0 = which(apply(pALLGlobalPhi.glm.nb[ ,c(6,8,10,12)],1,which.min) == 1)# Variants for which the Poisson model is better
  for (i in singhes0) {
    Xinter   = model.matrix(~cond*path,data= allEventtables[[which(rownames(dataPart3) == rownames(pALLGlobalPhi.glm.nb)[i])]]); 
    outinter = glmnet(Xinter,allEventtables[[which(rownames(dataPart3) == rownames(pALLGlobalPhi.glm.nb)[i])]]$counts,family="poisson",lambda=1e-4,alpha=0)
    Xprinc   = model.matrix(~path+cond,data= allEventtables[[which(rownames(dataPart3) == rownames(pALLGlobalPhi.glm.nb)[i])]]); 
    outprinc = glmnet(Xprinc,allEventtables[[which(rownames(dataPart3) == rownames(pALLGlobalPhi.glm.nb)[i])]]$counts,family="poisson",lambda=1e-4,alpha=0)
    Pv       = 1-pchisq(deviance(outprinc) - deviance(outinter),df=1)
    pALLGlobalPhi.glm.nb.glmnet$glmnet.pval[i] = Pv
    pALLGlobalPhi.glm.nb.glmnet$glmnet.code[i] = outinter$jerr
  }
  matrixpALLGlobalPhi.glmnet <- as.matrix(pALLGlobalPhi.glm.nb.glmnet)
  storage.mode(matrixpALLGlobalPhi.glmnet) <- 'numeric'

  #############################################################################################################
  ###  code chunk number 10.1 : Pseudo-counts  and glmnet                                                   ###
  #############################################################################################################
  singhes = which(apply(matrixpALLGlobalPhi[ ,c(6,8,10,12)],1,which.min) > 1 & apply(matrixpALLGlobalPhi[ ,c(22,24,26,28)],1,sum) != 0) 
  singhes_n = names(singhes) 

  pALLGlobalPhi.glm.nb.pen = as.data.frame(matrixpALLGlobalPhi)
  for(i in singhes){#we add pseudo-count (+1 in counts) for event with singular hessian for which the best model is not the Poisson model 
    pALLGlobalPhi.glm.nb.pen[i, ] = try(.fitNBglmModelsDSSPhi(.addOneCount(allEventtables[[i]]),
                                                            dispersion(dispData)[i],
                                                            dispersion(dispDataMeanCond)[i], phi, nbAll) ,silent=T)
  }

  pALLGlobalPhi.glm.nb <- as.data.frame(matrixpALLGlobalPhi)
  pALLGlobalPhi.glm.nb$final.pval.a.ia = 1
  i <- 1 #Poisson model
  li.singhes <- which(apply(pALLGlobalPhi.glm.nb[ ,c(6,8,10,12)],1,which.min) == i)
  pALLGlobalPhi.glm.nb$final.pval.a.ia[li.singhes] = pALLGlobalPhi.glm.nb.glmnet$glmnet.pval[li.singhes]
  
  i  <- 2 # negative binomial model, with global phi
  li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(6,8,10,12)],1,which.min) == i & apply(pALLGlobalPhi.glm.nb[ ,c(22,24,26,28)],1,sum) == 0)
  pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li,2]
  li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(6,8,10,12)],1,which.min) == i & apply(pALLGlobalPhi.glm.nb[ ,c(22,24,26,28)],1,sum) != 0)
  pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb.pen[li,2]

  i  <- 3 # negative binomial model, phi estimated with DSS
  li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(6,8,10,12)],1,which.min) == i & apply(pALLGlobalPhi.glm.nb[ ,c(22,24,26,28)],1,sum) == 0)
  pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li,3]
  li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(6,8,10,12)],1,which.min) == i & apply(pALLGlobalPhi.glm.nb[ ,c(22,24,26,28)],1,sum) != 0)
  pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb.pen[li,3]

  i  <- 4 # negative binomial model, phi estimated with DSS, conditionally to the expression mean
  li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(6,8,10,12)],1,which.min) == i & apply(pALLGlobalPhi.glm.nb[ ,c(22,24,26,28)],1,sum) == 0)
  pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li,4]
  li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(6,8,10,12)],1,which.min) == i & apply(pALLGlobalPhi.glm.nb[ ,c(22,24,26,28)],1,sum) != 0)
  pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb.pen[li,4]
  
  sing.events.final <- which(grepl("Error",pALLGlobalPhi.glm.nb[ , 29])) 
  if (length(sing.events.final) != 0) {
       pALLGlobalPhi.glm.nb <- pALLGlobalPhi.glm.nb[ - sing.events.final, ]
  }
  pALLGlobalPhi.glm.nb$final.padj.a.ia <- p.adjust(pALLGlobalPhi.glm.nb$final.pval.a.ia, method="fdr")

  #keeping only lines below pvalue in signifVariants:
  if (length(sing.events) != 0) {
    tmpdataPart3 <- dataPart3[ - sing.events, ]
    signifVariants <- cbind(tmpdataPart3[ - sing.events.final, ],pALLGlobalPhi.glm.nb$final.padj.a.ia )[ pALLGlobalPhi.glm.nb$final.padj.a.ia <= pvalue, ]
  } else {
    signifVariants <- cbind(dataPart3, pALLGlobalPhi.glm.nb$final.padj.a.ia )[ pALLGlobalPhi.glm.nb$final.padj.a.ia <= pvalue, ]
  }

  #############################################################################################################
  ###  code chunk number 11 :deltaPSI / deltaF computation
  #############################################################################################################
  if (! is.null(ASSBinfo)) {
    ASSBinfo <-subset(ASSBinfo,ASSBinfo$events.names %in% as.vector(signifVariants[,1])) #select only the lines corresponding to the remaining lines of signifVariants
  }
  sumLowCond <- matrix(data=rep(0,n*dim(signifVariants)[1]),nrow=dim(signifVariants)[1],ncol=n)#to  check later for low counts to flag
  pairsCond <-list()
  namesCond <- unique(sortedconditions)
  for ( i in 1:n) {
    j <- i+1
    while ( j <= n ) {
      pairsCond[[length(pairsCond)+1]] <- list(c(i,j))#creation of permutation of size 2 in c
      j <- j+1
    }
  } 
  namesPsiPairCond <- c()
  if (! is.null(ASSBinfo)) {
    newindex <- unlist(sapply(rownames(signifVariants),function(x)res <- which(ASSBinfo[,1]==x))) #to put the lines of the 2 data frames in the same order
    ASSBinfo <- ASSBinfo[newindex,]
  } else {
    rown <-row.names(signifVariants)
    lengths2 <- lengths[rown,]
  }
  deltapsi <- matrix(nrow=dim(signifVariants)[1], ncol=length(pairsCond))
  rownames(deltapsi) <- rownames(signifVariants)
  namesDeltaPsi <- c()
  for (pair in pairsCond) { #delta psi calculated for pairs of conditions, psi are calcuted for each replicateXcondition
      #for one pair
      index <- pair[[1]]
      condi <- namesCond[index]
      replicates <- nr[index]
      psiPairCond <- matrix(nrow=dim(signifVariants)[1],ncol=sum(replicates))
      colsPsiPairCond <- c()
      indexMatrixPsiPairCond <- 1
      indexdeltapsi <- 1
      for (nbRepli in 1:length(replicates)) {# for a given condition in the pair
        for (i in 1:replicates[nbRepli]) {# for each replicate (i) of the condition
          colsPsiPairCond <- c(colsPsiPairCond, paste(condi[nbRepli],'_r',i, sep=''))
          namesUp <- c(paste("UP_",condi[nbRepli],'_r',i,"_Norm", sep=''))
          namesLow <- c(paste("LP_",condi[nbRepli],'_r',i,"_Norm", sep=''))
          subsetUp <- signifVariants[namesUp]#the subsets are the counts we are going to use to compute all psis
          subsetLow <- signifVariants[namesLow]
          if (! is.null(ASSBinfo)) {# counts correction
            nameASSBinfo <- c(paste(condi[nbRepli],'_r',i, sep=''))
            subsetUp <- subsetUp/2-(ASSBinfo[,nameASSBinfo]/subsetUp)
          } else {#counts correction if there is no info about the junction counts
            subsetUp <- subsetUp/(lengths2$upper + readLength - 2*overlap + 1)
            subsetLow <- subsetLow/(lengths2$lower + readLength - 2*overlap + 1)
          }
          sumLowCond[,nbRepli] <- sumLowCond[,nbRepli] + as.matrix(subsetUp) + as.matrix(subsetLow)#sumLowCond sums up the counts for each condition
          subsetUp[which(subsetUp[,1] < 10),] <- NaN #NaN for counts <10
          subsetLow[which(subsetLow[,1] < 10),] <- NaN
          psiPairCond[,indexMatrixPsiPairCond] <- as.matrix(subsetUp/(subsetUp+subsetLow)) #psi is #incl/(#incl+#exclu) after all corrections
          indexMatrixPsiPairCond <- indexMatrixPsiPairCond + 1
          namesPsiPairCond <- c(namesPsiPairCond, as.character(condi[nbRepli]))
        }
      } 
  colnames(psiPairCond) <- namesPsiPairCond
  rownames(psiPairCond) <- rownames(signifVariants)
  rownames(sumLowCond) <- rownames(signifVariants)
  psiPairCond <- replace(psiPairCond, is.na(psiPairCond),0)
  NaNSums <- rowSums(( psiPairCond==0 )+0) #1 if NaN, 0 else
  #if there are 2 NaN and 3 values for a bcc, nanSums is at 2
  listNaN <- names(NaNSums[which(NaNSums>dim(psiPairCond)[2]/2)]) #when there are more NaN than nb of column/2, we do not calculate the psi
  psiPairCond[listNaN,]=NaN
  deltaPsiCond <- rowMeans(psiPairCond[,(replicates[1]+1):sum(replicates)]) - rowMeans(psiPairCond[,1:replicates[1]]) #delta psi is the mean of the psis of the 2nd condition (in terms of sorted condition) - the mean of the psis of the 1st condition 
  deltapsi[,indexdeltapsi] <- deltaPsiCond 
  indexdeltapsi <- indexdeltapsi + 1
  }
  #### when there are more than 2 conditions, we want to simplify the output :
  dPvector1 <-c(rep(0,dim(signifVariants)[1]))
  dPvector2 <-c(rep(0,dim(signifVariants)[1]))
  if (length(pairsCond) >1 ){
    for (l in 1:dim(deltapsi)[1]){ #if there are more than 2 conditions, we take the maximum of the deltaPSI of all pairs
      mindex <- which.max(abs(deltapsi[l,]))
      condA <- as.character(pairsCond[[mindex]][[1]][1])
      condB <- as.character(pairsCond[[mindex]][[1]][2])
      dP <- round(deltapsi[l,mindex],4)
      dPvector1[l] <- dP
      dPvector2[l] <-paste(as.character(dP),"(Cond",condB,",",condA,")",sep="")#we also return for which pair of conditions we found this max deltaPSI
    }
  } else {
    dPvector1 <- round(deltapsi,4)
    dPvector2 <- dPvector1
  }

  ###################################################
  ### code chunk 12 : final table
  ###################################################
  signifVariants <- cbind(signifVariants, dPvector1)
  signifVariants.sorted <- signifVariants[ order(-abs(dPvector1),signifVariants[dim(signifVariants)[2]-1]), ]#sorting by delta psi then by pvalue
  dPvector2.sorted <- dPvector2[order(-abs(dPvector1),signifVariants.sorted[dim(signifVariants.sorted)[2]-1])]
  signifVariants.sorted[dim(signifVariants.sorted)[2]] <- dPvector2.sorted

  colnames(signifVariants.sorted)[length(colnames(signifVariants.sorted))] <- 'Deltaf/DeltaPSI'# renaming last columns
  colnames(signifVariants.sorted)[length(colnames(signifVariants.sorted))-1] <- 'Adjusted_pvalue'
  rownames(signifVariants.sorted) <- signifVariants.sorted[,1]

  ###################################################
  ### code chunk 12.1 : flagging low counts
  ###################################################
  #Condition to flag a low count for an event :
  lowcounts <-apply(sumLowCond,1,function(x) length(which(x<flagLowCountsConditions))) >=n-1
  lowcounts <- lowcounts[rownames(signifVariants.sorted)] #to fit the order with the sorted order

  #### final tab ####
  signifVariants.sorted <- cbind(signifVariants.sorted, lowcounts)
  colnames(signifVariants.sorted[dim(signifVariants.sorted)[2]]) <- 'Low_counts'
  
  return(signifVariants.sorted)

  }