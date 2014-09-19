.lineParse <- function(lines, indexStart) {
  beginningLineToWrite <- ""
  splitElements <- strsplit(lines, "|", fixed=TRUE)[[1]] # splits the line
  if (indexStart == 6) {
    for (k in 1:(indexStart-2)) {
    beginningLineToWrite <- paste (beginningLineToWrite, splitElements[k], sep="|") #writes the firsts elements of the line : bcc, cycle... but NOT branching_nodes
    }
  } else {
    for (k in 1:(indexStart-1)) {
    beginningLineToWrite <- paste (beginningLineToWrite, splitElements[k], sep="|") #writes the firsts elements of the line : bcc, cycle... 
    }
  }
  endingLineToWrite <- paste(splitElements[length(splitElements)], sep="") #writes the rank
  ElementsNb <- length(splitElements) #number of elements in the line
  splitCounts <- splitElements[indexStart : (length(splitElements) - 1)] # avoids the name of the bcc ... and the rank (last one) to get only the counts
  s <- sapply(splitCounts, function(splitCounts) regmatches(splitCounts[[1]], gregexpr(pattern = "[0-9]+",splitCounts[[1]]))) #gets the junctions id (ex 1 in AS1) and the count (ex 6 in AS1_6)
  return(list(beginningLineToWrite,s,endingLineToWrite))
}

.countsSet <- function(lines, indexStart, counts=0, pairedEnd=FALSE, order=NULL, exonicReads=TRUE) {
  beginningLineToWrite <- .lineParse(lines, indexStart)[[1]]
  s <- .lineParse(lines, indexStart)[[2]] 
  endingLineToWrite <- .lineParse(lines, indexStart)[[3]]

  nbVec <- rep(0, length(s))
  countsVec <- rep(0, length(s))
  for (i in 1:length(s)) {
    nbVec[i] <- as.numeric(s[[i]][1])
    countsVec[i] <- as.numeric(s[[i]][2])
    if (counts > 1) {
      if (grepl("ASSB", names(s)[i]) == TRUE) { #so that counts on ASSB junction are not counted twice.
        countsVec[i] <- - countsVec[i]
      }
      
      if ((counts == 2) & (exonicReads == FALSE)) {
        if (grepl("^S[0-9]+", names(s)[i]) == TRUE) { #when exonic reads are not wanted we must discard reads counted in S_X
          countsVec[i] <- 0
        }
      }
    }
  }
  if (counts > 1) {
    d <- data.frame(nbVec,countsVec)
    names(d) <- c("NB", "COUNTS")
    sums <- aggregate(d$COUNTS, by=list(d$NB), sum) #sums the counts for each junction that belongs to the same event

    if (pairedEnd == TRUE) {
      if (is.null(order)) {
        order <- rep(1:((dim(sums)[1])/2), rep(2,((dim(sums)[1])/2)))
      } else {
        if (! is.vector(order)) {
          print("Error, order vector seems to be in a wrong format.")
        }
      }
      d2 <- data.frame(order, sums)
      dim(d2)
      names(d2)[3] <- 'sums'
      sums2 <- aggregate(d2$sums, by=list(d2$order), sum) # in case data is paired-end, there is one more sum to do, for each part of the pair
      sums <- sums2
    } 
  } else { ### counts == 0 
    if (pairedEnd == TRUE) {
      if (is.null(order)) {
        order <- rep(1:(length(s)/2), rep(2,(length(s)/2))) # for length(s)=8, will create a vector c(1,1,2,2,3,3,4,4) (assuming data is ordered)
      } else {
        if (! is.vector(order)) {
          print("Error, order vector seems to be in a wrong format.")
        }
      }
      d <- data.frame(order,countsVec)
      names(d) <- c("ORDER", "COUNTS")
      sums <- aggregate(d$COUNTS, by=list(d$ORDER), sum)
    }
  }
  lineToWrite <- ""
  for (j in 1:(dim(sums)[1]-1)) {
    lineToWrite <- paste(lineToWrite, sums[[1]][j], "_", sums[[2]][j], "|", sep="")
  }
  j <- j+1
  lineToWrite <- paste(lineToWrite, sums[[1]][j], "_", sums[[2]][j], sep="")
  lineToWrite <- paste(beginningLineToWrite, lineToWrite, endingLineToWrite, sep="|")
  lineToWrite <- substr(lineToWrite, start=2, stop=nchar(lineToWrite)) # lineToWrite looks like ">bcc_XXXX|Cycle_Y|Type_1|upper_path_length_ZZ|1_0|2_47|3_1|4_13|rank_1.00000"
  return(lineToWrite)
}

.addOneCount <- function(df)
{
  df$counts <- unlist(lapply(df[,'counts'],function(x){x+1}))
  return (df)
}


.replaceCounts <- function(lines, counts=0, pairedEnd=FALSE, order=NULL, exonicReads=TRUE) {
  index <- 1
  firstLineChar <- substr(lines[index], start = 0, stop = 1)
  if (grepl("branching_nodes", lines[index])) {
    indexStart <- 6
  } else {
    indexStart <- 5
  }
  lineToWrite <- c()
  if (firstLineChar == '>') {
    while (index <= length(lines)) {
      if (index%%2 == 1) {
        if (counts > 0) {
          lineToWrite <- c(lineToWrite, .countsSet(lines[index], indexStart, counts, pairedEnd, order, exonicReads))
        } else {
          if (pairedEnd == TRUE) {
            lineToWrite <- c(lineToWrite, .countsSet(lines[index], indexStart, counts, pairedEnd, order))
          }
        }
      
    } else {
      lineToWrite <- c(lineToWrite, lines[index])
    }
    index <- index + 1
    }
  }
  return(lineToWrite)
}

.readAndPrepareData <- function(countsData,conditions) {
  ###################################################
  ### code chunk number 1: Read data
  ###################################################
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
  }
  countsData[,-(1:2)] = countsData[,sortedindex]
  colnames(countsData) <- namesData
  countsData$Path <- gl( 2, 1, dim(countsData)[1], labels = c("UP", "LP"))

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
  cds <- newCountDataSet( countsData[ ,3:(3+length(conds)-1)], conds ) # create object
  cdsSF <- estimateSizeFactors(cds)
  sizeFactors( cdsSF )
  shouldWeNormalise=sum(is.na(sizeFactors(cdsSF))) < 1
  dim <- dim(countsData)[2]
  countsData[ ,(dim+1):(dim+length(conds)) ] <- round(counts(cdsSF, normalized=shouldWeNormalise))
  colnames(countsData)[(dim+1):(dim+length(conds))] <- paste(namesData[3:(3+sum(nr)-1)],"_Norm",sep="")
  return(list(countsData, conds, dim, n, nr, sortedconditions))
}

.eventtable <- function(df,startPosColumn4Counts, endPosCol4Counts){
  eventTab = data.frame(ID=rep(as.factor(df['ID']), endPosCol4Counts-startPosColumn4Counts+1),
  cond=as.factor(unlist(lapply(strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts],"_"),FUN=function(d){d[2]}))),
  counts=as.numeric(df[startPosColumn4Counts:endPosCol4Counts]),
  path=as.factor(unlist(lapply(strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts],"_"),FUN=function(d){d[1]}))),
  row.names = NULL)
  return(eventTab)
}

.fitNBglmModelsDSSPhi <- function(eventdata, phiDSS, phiDSScond, phiGlobal, nbAll){
    # S: simple, A: additive, I : interaction models
    # Poisson model  
  nbglmA0 <- negbin(counts~cond + path, data=eventdata, random=~1, fixpar=list(4,0))
  nbglmI0 <- negbin(counts~cond * path, data=eventdata, random=~1, fixpar=list(5,0))  

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

  isQuality <- grepl(toConvert[1],)
  line <- readLines(toConvert)
  if ( pairedEnd == TRUE || counts > 1 || grepl("Q",line[1]) ) {#if the data is paired-end or there is a count option or quality information (SNPS) we need a special treatment of data
    line <- .replaceCounts(line, counts, pairedEnd, order, exonicReads)
    } 

  index <- 1
  indexNames <- 1
  if (substr(line[index], start = 0, stop = 1) == '>') {
    len <- length(strsplit(line[index], "|", fixed = TRUE)[[1]])-6
  }
  events.mat <- matrix(NA,length(line)/2,len+2)
  events.names <- rep(NA,length(line)/2)
  firstLineChar <- substr(line[index], start = 0, stop = 1)
  if (firstLineChar == '>') {#checks if the line contains the header beginning with ">", not the sequence
    while (index <= length(line)) {
      nbCharLine <- length(strsplit(line[index],"|")[[1]]) # number of characters in the line
      lineInfos <- substr(line[index],start = 2,stop =nbCharLine) #the line without ">" that we want to avoid
      lineSplit <- strsplit(lineInfos, "|", fixed=TRUE)[[1]] # gets pieces of information separated by "|" in KisSplice format
      #example of lineSplit :  
      # lineSplit[1] : "bcc_3929"
      #          [2] : "Cycle_0"
      #          [3] : "Type_1"     
      #          [4] : "upper_path_length_90" 
      #          [5] : "C1_1" (first condition)
      #          [5+len] : "Cn_5" (last condition)              
      #          [5+len+1] : rank_0.90267" 
      eventName <- paste(lineSplit[1],lineSplit[2],sep="_") # concatenates the two first elements of lineSplit, ie the event name
      events.names[indexNames] <- eventName
      lengthInfo <- lineSplit[4]
      events.mat[indexNames,1] <- as.numeric(strsplit(lengthInfo,"_")[[1]][4]) #fills the first column of the matrix with length info
      lineCounts <- lineSplit[5:(5+len)] # gets every condition of the line and its associated count
      events.mat[indexNames,2:dim(events.mat)[2]] <- as.numeric(lapply(lineCounts, function(x) strsplit(x, "_")[[1]][2])) #fills the matrix others columns with the counts of the conditions of the line 
    index <- index + 2
    indexNames <- indexNames + 1
    }
  }
  class(events.mat) <- "numeric"
  events.df <- as.data.frame(events.mat)
  events.df <- data.frame(events.names,events.df)
  close(toConvert)
  return (events.df)
}

qualityControl <- function(countsData,conditions,storeFigs=FALSE) {
  options(warn=-1) # suppress the warning for the users

  if (storeFigs == TRUE){
    d<-system("find -type d -name Figures", TRUE)
    if (length(d) == 0) {
      system("mkdir Figures")
        } else if (d != "./Figures") {
          system("mkdir Figures") 
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
        png(filename="Figures/dendrogram.png")
        plot(hclust(as.dist(1-cor(countsData[ ,(dim+1):(dim+length(conds))])),"ward"))
        dev.off()
    }

  ###################################################
  ### code chunk number 3: replicates
  ###################################################
  if (storeFigs == FALSE) {
    heatmap(as.matrix(as.dist(1-cor(countsData[ ,(dim+1):(dim+length(conds))]))), margins = c(10,10))
    } else {
        png(filename="Figures/heatmap.png")
        heatmap(as.matrix(as.dist(1-cor(countsData[ ,(dim+1):(dim+length(conds))]))), margins = c(10,10))
        dev.off()
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
        png(filename="Figures/varInterIntra.png")
        plot( x = countsData$varIntra, y = countsData$varInter, xlab = "Intra-variability", ylab = "Inter-variability", las = 1, log = "xy")
        abline( a = 0, b = 1, col = 2, lty = 2, lwd = 2 )
        dev.off()
    }
}

diffExpressedVariants <- function(countsData, conditions, storeFigs=FALSE, pvalue=0.05) {
  if (storeFigs == TRUE){
    d<-system("find -type d -name Figures", TRUE)
    if (length(d) == 0) {
      system("mkdir Figures")
        } else if (d != "./Figures") {
          system("mkdir Figures") 
        }
  }

  options(warn=-1) # suppress the warning for the users
  ###################################################
  ### code chunk number 1: Read and prepare data
  ###################################################
  listData <-.readAndPrepareData(countsData,conditions)
  countsData <- listData[[1]]
  n <- listData[[4]]
  nr <- listData[[5]]
  sortedconditions <- listData[[6]]

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
  dataPart2[2] <- dataPart[seq(1,dim(dataPart)[1],2),2]-dataPart[seq(2,dim(dataPart)[1],2),2] # computes the difference of length between the lower and upper paths
  names(dataPart2)[2] <- "Length_diff"
  dataPart2 <- dataPart2[ ,c(-(3+nbAll))]
  if (anyDuplicated(dataPart2[ ,1]) > 0) {
    dataPart2 <- dataPart2[!duplicated(as.character(dataPart2[ ,1])), ]
  }
  rownames(dataPart2)=as.character(dataPart2[ ,1])

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
  ## estimate the overdispersion parameter theta of the NB distritution
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
        png(filename="Figures/models.png")
        plot(event.mean.variance.df$Mean, event.mean.variance.df$Variance, 
        xlab="Mean Event count", 
        ylab="Variance Event count",
        log="xy", las=1)
        abline(a=0, b=1, col=2, lwd=2)
        lines(x,yQP,col=3, lwd=2)
        lines(x,yNB,col=6, lwd=2)
        legend("topleft", c("Poisson","Quasi-Poisson", "Negative Binomial"), 
        text.col=c(2,3,6), box.lty=0);
        dev.off()
    }

  totLOW <- as.vector(apply(dataPart2[ ,(3 + sum(nr)):(3 + 2 * sum(nr) - 1)],1,sum)) #global counts for each variant (low/up) by event
  totUP <- as.vector(apply(dataPart2[ ,3:(3 + sum(nr) - 1)],1,sum))

  dataPart3 <- dataPart2[-which(totUP <10| totLOW<10),]#after the dispersion estimation, discard the events that have at least one variant with global count <10
  allEventtables  <- apply(dataPart3,1,.eventtable, startPosColumn4Counts = which(grepl("UP",names(dataPart3)))[1],endPosCol4Counts = ncol(dataPart3))
  ###################################################
  ### code chunk number 6: pALLGlobalPhi.glm.nb
  ###################################################
  pALLGlobalPhi.glm.nb=data.frame(t(rep(NA,28)))#mettre 0 a la place des NA
  for (i in 1:length(allEventtables)) {
    pALLGlobalPhi.glm.nb[i, ] = try(.fitNBglmModelsDSSPhi(allEventtables[[i]],dispersion(dispData)[i],dispersion(dispDataMeanCond)[i], phi, nbAll) ,silent=T)
  }
  ###################################################
  ### code chunk number 7: excl_errors
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


  #####
  matrixpALLGlobalPhi <- as.matrix(pALLGlobalPhi.glm.nb)
  storage.mode(matrixpALLGlobalPhi) <- 'numeric'
  #####
  #####

  ###################################################
  ### code chunk number 10: best model
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
  bestmodel3 = table(bestmodel.table,apply(matrixpALLGlobalPhi[ ,c(22,24,26,28)],1,sum),dnn=c("best model","singular hessian"))# ici pbl sum car character %%%%%%%%%%
  colnames(bestmodel3) = paste(colnames(bestmodel3), "Models")

  ###################################################
  ### code chunk number 11: glmnet
  ###################################################
  #pALLGlobalPhi.glm.nb.glmnet = pALLGlobalPhi.glm.nb ####
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
  ### Pseudo-count for event with singular hessian for which the best model is not the Poisson model ###
  #############################################################################################################
  singhes = which(apply(matrixpALLGlobalPhi[ ,c(6,8,10,12)],1,which.min) > 1 & apply(matrixpALLGlobalPhi[ ,c(22,24,26,28)],1,sum) != 0) ######%%%%%
  singhes_n = names(singhes) 

  pALLGlobalPhi.glm.nb.pen = as.data.frame(matrixpALLGlobalPhi)#########%%%%%%%%%%%%%
  for(i in singhes){
    pALLGlobalPhi.glm.nb.pen[i, ] = try(.fitNBglmModelsDSSPhi(.addOneCount(allEventtables[[i]]),
                                                            dispersion(dispData)[i],
                                                            dispersion(dispDataMeanCond)[i], phi, nbAll) ,silent=T)
  }

  # sing.events.pseudocounts <- which(grepl("Error",pALLGlobalPhi.glm.nb.pen[ , 1])) ########%%%%%%%%
  # if (length(sing.events.pseudocounts) != 0) {
  #     pALLGlobalPhi.glm.nb.pen <- pALLGlobalPhi.glm.nb.pen[ - sing.events.pseudocounts, ]
  # }

  # matrixpALLGlobalPhi.pen <- as.matrix(pALLGlobalPhi.glm.nb.pen)
  # storage.mode(matrixpALLGlobalPhi.pen) <- 'numeric'
  #singhes2 = which(apply(matrixpALLGlobalPhi.pen[ ,c(6,8,10,12)],1,which.min) > 1 & apply(matrixpALLGlobalPhi.pen[ ,c(22,24,26,28)],1,sum) != 0)######%%%%

  ###################################################
  ###################################################

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
  
  #############
  sing.events.final <- which(grepl("Error",pALLGlobalPhi.glm.nb[ , 29])) ########%%%%%%%%
  if (length(sing.events.final) != 0) {
       pALLGlobalPhi.glm.nb <- pALLGlobalPhi.glm.nb[ - sing.events.final, ]
  }
  #######

  pALLGlobalPhi.glm.nb$final.padj.a.ia <- p.adjust(pALLGlobalPhi.glm.nb$final.pval.a.ia, method="fdr")
  if (length(sing.events) != 0) {
    tmpdataPart3 <- dataPart3[ - sing.events, ] ######## %%%%%%%%%%
    signifVariants <- cbind(tmpdataPart3[ - sing.events.final, ],pALLGlobalPhi.glm.nb$final.padj.a.ia )[ pALLGlobalPhi.glm.nb$final.padj.a.ia <= pvalue, ]########### %%%%%%%%
  } else {
    signifVariants <- cbind(dataPart3, pALLGlobalPhi.glm.nb$final.padj.a.ia )[ pALLGlobalPhi.glm.nb$final.padj.a.ia <= pvalue, ]
  }
  # sorting by deltaPSI / deltaF
  finalDelta <- NA
  #1st condition
  PSI.replicat1 <- c()
  tmp.psi <-c()
  for (j1 in 1:nr[1]) {
    nameUp <- paste('UP_',sortedconditions[j1], '_r', j1, '_Norm', sep='')
    nameLow <- paste('LP_', sortedconditions[j1],'_r', j1, '_Norm', sep='') 
    for (indexSignif in 1:dim(signifVariants)[1]) {
      if (signifVariants[indexSignif,nameUp] + signifVariants[indexSignif,nameLow] != 0) {
        tmp.psi[indexSignif] <- signifVariants[indexSignif, nameUp] / (signifVariants[indexSignif,nameUp] + signifVariants[indexSignif,nameLow])
      } else {
        tmp.psi[indexSignif] <- 0
      }
    }
    PSI.replicat1 <- cbind( PSI.replicat1, tmp.psi)
    PSIcondI <- apply(PSI.replicat1, MARGIN=1, mean, na.rm=T)
    for (k1 in (2:n)) { # = 2nd condition
      PSI.replicat2 <- c()
      for (l in (1:nr[n])) {
        nameUp <- paste('UP_',sortedconditions[cumsum(nr)[n-1]+l], '_r', l, '_Norm', sep='')
        nameLow <- paste('LP_', sortedconditions[cumsum(nr)[n-1]+l],'_r', l, '_Norm', sep='') 
        for (indexSignif in 1:dim(signifVariants)[1]) {
          if (signifVariants[indexSignif,nameUp] + signifVariants[indexSignif,nameLow] != 0) {
            tmp.psi[indexSignif] <- signifVariants[indexSignif, nameUp] / (signifVariants[indexSignif,nameUp] + signifVariants[indexSignif,nameLow])
          } else {
            tmp.psi[indexSignif] <- 0
          }
        }
        PSI.replicat2 <- cbind( PSI.replicat2, tmp.psi)
      }
      PSIcondJ <- apply(PSI.replicat2, MARGIN=1, mean, na.rm=T)
      #deltaPSI for cond i,j
      deltaPSIij <- PSIcondJ- PSIcondI
      finalDelta <- apply( cbind( finalDelta, deltaPSIij) , MARGIN=1, max, na.rm=T)
    }
  }

  if (n > 2) {
    tmp.psi <-c()
    for (i in 2:(n-1)) {
    #1st condition
      PSI.replicat1 <- c()
      for (k2 in (1:nr[n-1])) {    
        nameUp <- paste('UP_',sortedconditions[cumsum(nr)[i-1]+k2], '_r', k2, '_Norm', sep='')
        nameLow <- paste('LP_', sortedconditions[cumsum(nr)[i-1]+k2],'_r', k2, '_Norm', sep='') 
        for (indexSignif in 1:dim(signifVariants)[1]) {
          if (signifVariants[indexSignif,nameUp] + signifVariants[indexSignif,nameLow] != 0) {
            tmp.psi[indexSignif] <- signifVariants[indexSignif, nameUp] / (signifVariants[indexSignif,nameUp] + signifVariants[indexSignif,nameLow])
          } else {
            tmp.psi[indexSignif] <- 0
          }
        }
        PSI.replicat1 <- cbind( PSI.replicat1, tmp.psi)
      }
      PSIcondI <- apply(PSI.replicat1, MARGIN=1, mean, na.rm=T)
      for (k3 in (i+1):n) {
        PSI.replicat2 <- c()
        for (l in 1:nr[n]) {
          nameUp <- paste('UP_',sortedconditions[cumsum(nr)[k3-1]+l], '_r', l, '_Norm', sep='')
          nameLow <- paste('LP_', sortedconditions[cumsum(nr)[k3-1]+l],'_r', l, '_Norm', sep='') 
          for (indexSignif in 1:dim(signifVariants)[1]) {
            if (signifVariants[indexSignif,nameUp] + signifVariants[indexSignif,nameLow] != 0) {
              tmp.psi[indexSignif] <- signifVariants[indexSignif, nameUp] / (signifVariants[indexSignif,nameUp] + signifVariants[indexSignif,nameLow])
            } else {
              tmp.psi[indexSignif] <- 0
            }
          }
          PSI.replicat2 <- cbind( PSI.replicat2, tmp.psi)
        }
        PSIcondJ <- apply(PSI.replicat2, MARGIN=1, mean, na.rm=T)
        # deltaPSI for cond i,j
        deltaPSIij <- PSIcondJ- PSIcondI
        finalDelta <- apply( cbind( finalDelta, deltaPSIij) , MARGIN=1, max, na.rm=T)
      }
    }
  }
  

  signifVariants <- cbind(signifVariants, finalDelta)# adding DeltaPsi/f to the final table
  colnames(signifVariants)[length(colnames(signifVariants))] <- 'Deltaf/DeltaPSI'# renaming last columns
  colnames(signifVariants)[length(colnames(signifVariants))-1] <- 'Adjusted_pvalue'# renaming last columns
  signifVariants.sorted <- signifVariants[order( abs(finalDelta), decreasing=T), ]



  lowcounts <- c()
 
###################################################
### Low counts
###################################################
#Condition to flag a low count for an event :
   # if at least n-1 conditions have counts <10 


#conditions
todo1 <- 3
todo2 <- 0
done <- 0
vectCond <- c()
for (i in 1:length(nr)) { #calculating the total count per condition (summing by variants and replicates) per event
  todo1 <- todo1 + done
  done <- nr[i]
  todo2 <- todo1 + done - 1

 sums <- apply(signifVariants.sorted[,todo1:todo2],1,sum) + apply(signifVariants.sorted[,(todo1 + sum(nr)):(todo2 + sum(nr))],1,sum) #up +low for 1 condition
 vectCond <- c(vectCond,sums)
}

m <- matrix(vectCond, ncol = n)
totCOND <- c()
for (i in 1:dim(m)[1]){
  totCOND <- c(totCOND,length(m[i, m[i, ]<10]) >= n-1) #at least n-1 conditions have counts below 10
} 

# lowcounts <- totUP <10 | totLOW <10 | totCOND
lowcounts <- totCOND
signifVariants.sorted <- cbind(signifVariants.sorted, lowcounts)
colnames(signifVariants.sorted[dim(signifVariants.sorted)[2]]) <- 'Low_counts'
return(signifVariants.sorted)

}