kissplice2counts <- function(fileName, counts = 0, pairedEnd = FALSE, order = NULL, exonicReads = TRUE, discoSNP = FALSE, k2rg = FALSE, keepInDel = TRUE) {
  toConvert <- file(fileName, open = "r")
  lines <- readLines(toConvert)
  cDupBcc <- c()
  if (k2rg == FALSE) {
    line <- lines[1]
    isQuality <- grepl("Q", line[1])
    resultLine1 <- .getInfoLine(line, counts, pairedEnd, order, exonicReads, isQuality, discoSNP)  # get all the informations for the 1st line
    eventName <- resultLine1$eventName
    variantLength <- resultLine1$variantLength
    variantCounts <- resultLine1$variantCounts
    events.mat <- matrix(NA, length(lines) / 2, length(variantCounts) + 1)
    events.names <- rep(NA, length(lines) / 2)
    events.mat[1, 1] <- as.numeric(variantLength)
    events.mat[1, 2:NCOL(events.mat)] <- variantCounts
    events.names[1] <- eventName
    index <- 3
    indexNames <- 2
    firstLineChar <- substr(lines[index], start = 0, stop = 1)
    psiInfo <- matrix(NA, length(lines) / 2, length(resultLine1$psiInfo))
    psiInfo[1, ] <- resultLine1$psiInfo
    #### 6/09 ####
    # nbcol <- length(variantCounts) + 2
    # events.df <- as.data.frame(matrix(nrow = length(lines) / 2, ncol = nbcol))
    ####
    if (firstLineChar == ">") {  # same for all other lines, ignore lines with sequences
      while (index <= length(lines)) {
        line <- lines[index]
        resultLine <- .getInfoLine(line, counts, pairedEnd, order, exonicReads)
        eventName <- resultLine$eventName
        variantLength <- resultLine$variantLength
        variantCounts <- resultLine$variantCounts
        events.mat[indexNames, 1] <- as.numeric(variantLength)
        events.mat[indexNames, 2:NCOL(events.mat)] <- variantCounts
        events.names[indexNames] <- eventName
        psiInfo[indexNames, ] <- resultLine$psiInfo
        index <- index + 2
        indexNames <- indexNames + 1
        class(events.mat) <- "numeric"
        #### 6/09 ####
        # events.df[, 1] <- events.names
        # events.df[, 2:nbcol] <- events.mat
        
        # events.df <- as.data.frame(events.mat)
        # events.df <- data.frame(events.names, events.df)
        ####
      }
    }
    # events.df <- as.data.frame(events.mat)
    events.df <- data.frame(events.names, events.mat)
    
  } else {
    GENEID <- 1
    GENENAME <- 2
    POS <- 3
    STRAND <- 4
    EVENT <- 5
    VARPARTLENGTH <- 6
    FRAMESHIFT <- 7
    CDS <- 8
    GENEBIOTYPE <- 9
    SPLICESITE <- 10
    BLOCSIZEUP <- 11
    SPLICESITEPOSUP <- 12
    PARALOGS <- 13
    COMPLEX <- 14
    SNPVARREGION <- 15
    EVENTNAME <- 16
    BLOCSIZELOW <- 17
    SPLICESITEPOSLOW <- 18
    PSIS <- 19
    COVERAGEUP <- 20
    COVERAGELOW <- 21
    CANONICAL <- 22
    
    if(keepInDel)
      FILESUFIX <- "_no_duplicate"
    else
      FILESUFIX <- "_no_indel_no_duplicate"
    fileSplit <- strsplit(fileName, split = "\\.")

    if(length(fileSplit[[1]]) > 1)
      FILE <- paste(fileSplit[[1]][1], FILESUFIX, ".", fileSplit[[1]][2], sep = "")
    else
      FILE <- paste(fileSplit[[1]][1], FILESUFIX, sep = "")
    
    i <- 1
    line <- lines[i]
    write(line, file = FILE)
    i <- 2
    line <- lines[i]
    filterOut <- c("deletion", "insertion")
    resultLine <- .getInfoLineK2rg(line, counts, pairedEnd, order, exonicReads)
    variantCountsUp <- resultLine$variantCountsUp
    iBcc <- 1
    lBcc <- list()
    while (i <= length(lines)) {
      if ((!strsplit(line, split = "\t")[[1]][16] %in% lBcc) && (keepInDel || !strsplit(line, split = "\t")[[1]][5] %in% filterOut)) {
        lBcc[iBcc] <- strsplit(line, split = "\t")[[1]][16]
        iBcc <- iBcc + 1
        write(line, file = FILE, append = TRUE)
      }
      i <- i + 1
      line <- lines[i]
    }
    events.mat <- matrix(NA, iBcc * 2 - 2, length(variantCountsUp) + 1)
    events.names <- rep(NA, iBcc * 2 - 2)
    psiInfo <- matrix(NA, iBcc * 2 - 2, length(resultLine$psiInfoUp))
    indexNames <- 1
    matBccApp <- matrix(0,nrow = length(lBcc)) # nombre d'apparition pour chaque BCC
    rownames(matBccApp) <- lBcc
    iDupBcc <- 1
    i <- 2
    while (i <= length(lines)) {
      line <- lines[i]
      if(keepInDel || !(strsplit(line, split = "\t")[[1]][5] %in% filterOut)){
        bcc <- strsplit(line, split = "\t")[[1]][16]
        matBccApp[bcc, 1] <- matBccApp[bcc, 1] + 1
        
        if (matBccApp[bcc, 1] > 1) {
          if(!bcc %in% cDupBcc) {
            cDupBcc[iDupBcc] <- bcc
            iDupBcc <- iDupBcc + 1
          }
          
        } else {
          resultLine <- .getInfoLineK2rg(line, counts, pairedEnd, order, exonicReads)
          eventName <- resultLine$eventName
          variantLength <- resultLine$variantLength
          variantCountsUp <- resultLine$variantCountsUp
          variantCountsLow <- resultLine$variantCountsLow
          events.mat[indexNames, 1] <- as.numeric(variantLength)
          events.mat[indexNames, 2:NCOL(events.mat)] <- variantCountsUp
          events.names[indexNames] <- eventName
          psiInfo[indexNames, ] <- resultLine$psiInfoUp
          events.mat[indexNames + 1, 1] <- 0
          events.mat[indexNames + 1, 2:NCOL(events.mat)] <- variantCountsLow
          events.names[indexNames + 1] <- eventName
          psiInfo[indexNames + 1, ] <- resultLine$psiInfoLow
          class(events.mat) <- "numeric"
          indexNames <- indexNames + 2
        }
      }
      i <- i + 1
      # events.df <- as.data.frame(events.mat)
      # events.df <- data.frame(events.names, events.df)
    }
    events.df <- data.frame(events.names, events.mat)
  }
  
  close(toConvert)
  psiInfo <- data.frame(events.names, as.data.frame(psiInfo))
  if (length(cDupBcc) > 0)
    dupBcc.df <- data.frame(matBccApp[cDupBcc, 1])
  else
    dupBcc.df <- data.frame()

  output <- list(countsEvents = events.df, psiInfo = psiInfo, discoInfo = discoSNP, dupBcc = dupBcc.df)
  class(output) <- c("list", "countsData")
  return(output)
}


print.countsData <- function(x, ...){
  print(x$countsEvents)
}


qualityControl <- function(countsData, conditions, storeFigs = FALSE) {
  
  options(warn = -1)  # suppress the warning for the users
  
  if(storeFigs == FALSE)
    pathToFigs <- NA
  else if (isTRUE(storeFigs))
    pathToFigs <- "kissDEFigures"
  else
    pathToFigs <- storeFigs
  
  # create a new folder if it doesn't exist
  if (!is.na(pathToFigs)){
    find <- paste("find", pathToFigs)
    d <- system(find, TRUE, ignore.stderr = TRUE)
    if (length(d) == 0) { 
      command <- paste("mkdir", pathToFigs)
      system(command, ignore.stderr = TRUE)
    }
  }
  
  ###################################################
  ### code chunk number 1: Read and prepare data
  ###################################################
  listData <- .readAndPrepareData(countsData, conditions)
  countsData <- listData$countsData
  conds <- listData$conditions
  dimns <- listData$dim
  n <- listData$n
  nr <- listData$nr
  
  ###################################################
  ### code chunk number 2: dendrogram
  ###################################################
  if (storeFigs == FALSE) {
    plot(hclust(as.dist(1 - cor(countsData[, (dimns + 1):(dimns + length(conds))])), "ward.D"))
    par(ask = TRUE)
  } else {
    filename <- paste(storeFigs, "/dendrogram.png", sep = "")
    png(filename)
    plot(hclust(as.dist(1 - cor(countsData[, (dimns + 1):(dimns + length(conds))])), "ward.D"))
    void <- dev.off()
  }
  
  ###################################################
  ### code chunk number 3: replicates
  ###################################################
  if (storeFigs == FALSE) {
    heatmap(as.matrix(as.dist(1 - cor(countsData[, (dimns + 1):(dimns + length(conds))]))), margins = c(10, 10))
  } else {
    filename <- paste(storeFigs, "/heatmap.png", sep = "")
    png(filename)
    heatmap(as.matrix(as.dist(1 - cor(countsData[, (dimns + 1):(dimns + length(conds))]))), margins = c(10, 10))
    void <- dev.off()
  }
  
  ###################################################
  ### code chunk number 4: intra-group and inter-group-variance
  ###################################################
  # Mean and variance over all conditions and replicates (normalized counts!)
  countsData$mn <- rowMeans(countsData[, (dimns + 1):(dimns + length(conds))])
  countsData$var <- apply(countsData[, (dimns + 1):(dimns + length(conds))], 1, var)
  # correction term
  nbAll <- sum(nr)  # number of all observations in all groups
  countsData$ct <- rowSums(countsData[, (dimns + 1):(dimns + length(conds))])^2 / nbAll
  # sum of squares between groups
  countsData$ss <- rowSums(countsData[, (dimns + 1):(dimns + length(conds) / n)])^2 / nr[1] + rowSums(countsData[, ((dimns + 1) + length(conds) / n):(dimns + length(conds))])^2 / nr[2]
  # substract the correction term from the SS and divide by the degrees of 
  df <- 1 # freedom(groups); here: df=2-1=1
  countsData$varInter <- (countsData$ss - countsData$ct) / df
  # intra-variability 
  countsData$varC1 <- apply(countsData[, (dimns + 1):(dimns + nr[1])], 1, var)
  countsData$varC2 <- apply(countsData[, ((dimns + 1) + nr[1]):(dimns + nr[2] + nr[1])], 1, var)
  countsData$varIntra <- rowMeans(data.frame(countsData$varC1, countsData$varC2))
  
  ###################################################
  ### code chunk number 5: intra-vs-inter
  ###################################################
  if (storeFigs == FALSE) {
    plot(x = countsData$varIntra, y = countsData$varInter, xlab = "Intra-variability", ylab = "Inter-variability", las = 1, log = "xy")
    abline(a = 0, b = 1, col = 2, lty = 2, lwd = 2)
  } else {
    filename <- paste(storeFigs, "/InterIntraVariability.png", sep = "")
    png(filename)
    plot(x = countsData$varIntra, y = countsData$varInter, xlab = "Intra-variability", ylab = "Inter-variability", las = 1, log = "xy")
    abline(a = 0, b = 1, col = 2, lty = 2, lwd = 2)
    void <- dev.off()
  }
}



diffExpressedVariants <- function(countsData, conditions, storeFigs = FALSE, pvalue = 0.05, filterLowCountsVariants = 10, flagLowCountsConditions = 10, discoSNP = FALSE) {
  
  options(warn = -1)  # suppress the warning for the users
  if(storeFigs == FALSE)
    pathToFigs <- NA
  else if (isTRUE(storeFigs))
    pathToFigs <- "kissDEFigures"
  else
    pathToFigs <- storeFigs
  
  print("Pre-processing the data...")
  chunk0 <- tryCatch({.readAndPrepareData(countsData, conditions)
    #### chunk 0 var ####
    # chunk0$countsData
    # chunk0$conditions
    # chunk0$dim
    # chunk0$n
    # chunk0$nr
    # chunk0$sortedconditions
    # chunk0$ASSBinfo
  }, error = function(err) {
    print(err)
    return(NA)
  })
  
  if (!is.na(chunk0)){  # no error in chunk 0
    ASSBinfo <- chunk0$ASSBinfo  # in case counts option in kissplice2counts is at 1 or 2, we have info about junction counts (ASSB), that will be useful to correct the computation of delta psi in the end. They are stored here.
    if (!is.null(ASSBinfo)) {
      li <- c()
      for (i in (1:NROW(ASSBinfo))){
        if (i%%2 != 0) {
          li <- c(li, i)
        }
      }
      ASSBinfo <- ASSBinfo[li, ]
    }
    print("Trying to fit models on data...")
    chunk1 <- tryCatch({.modelFit(chunk0$countsData, chunk0$n, chunk0$nr, ASSBinfo, storeFigs, pathToFigs, filterLowCountsVariants)
      #### chunk 1 var ####
      # chunk1$pALLGlobalPhi.glm.nb 
      # chunk1$sing.events
      # chunk1$dataPart3
      # chunk1$ASSBinfo
      # chunk1$allEventtables
      # chunk1$length
      # chunk1$phi
      # chunk1$dispData
    }, error = function(err) {
      print(paste(err, "An error occured, unable to fit models on data." ))
      return(NA)
    }) 
  } else {  # error in chunk 0
    chunk1 <- NA
  }
  
  if (!is.na(chunk1)) {  # no error in chunk 1 nor in chunk 0
    print("Searching for best model and computing pvalues...")
    chunk2 <- tryCatch({.bestModelandSingular(chunk1$pALLGlobalPhi.glm.nb, chunk1$sing.events, chunk1$dataPart3, chunk1$allEventtables, pvalue, chunk1$phi, chunk0$nr, chunk1$dispData)
      #### chunk 2 var ####  
      # chunk2$noCorrectPVal
      # chunk2$correctedPVal
      # chunk2$signifVariants
    }, error = function(err) {
      print(paste(err, "Returning only resultFitNBglmModel and sing.events")) 
      return(list(resultFitNBglmModel = chunk1$pALLGlobalPhi.glm.nb, sing.events = chunk1$sing.events))
    })
  } else {
    chunk2 <- NA
  }
  
  if (!is.na(chunk2)) {  # no error during chunk1
    if (length(chunk2) > 2) {  # no error during chunk2
      print("Computing size of the effect and last cutoffs...")
      chunk3 <- tryCatch({
        if (discoSNP != FALSE){
          if (is.na(countsData$discoInfo)){
            discoSNP <- FALSE
          } else {
            discoSNP <- countsData$discoInfo
          }
        } 
        
        signifVariants.sorted <- .sizeOfEffectCalc(chunk2$signifVariants, chunk1$ASSBinfo, chunk0$n, chunk0$nr, chunk0$sortedconditions, 
                                                   flagLowCountsConditions, chunk1$lengths, discoSNP)
        return(list(finalTable = signifVariants.sorted, 
                    correctedPVal = chunk2$correctedPVal, 
                    uncorrectedPVal = chunk2$noCorrectPVal, 
                    resultFitNBglmModel = chunk1$pALLGlobalPhi.glm.nb))
      }, error = function(err) {
        print(paste(err, "Returning only resultFitNBglmModel and pvalues tab"))
        return(list(correctedPVal = chunk2$correctedPVal,
                    uncorrectedPVal = chunk2$noCorrectPVal,
                    resultFitNBglmModel = chunk1$pALLGlobalPhi.glm.nb))
      })
    } else {  # error in chunk 2 does not allow to compute chunk 3
      return(chunk2)
    }
  } else {
    return(NA)
  }
}
