kissplice2counts <- function(fileName, counts = 0, pairedEnd = FALSE, order = NULL, exonicReads = TRUE, k2rg = FALSE, keep = c("All"), remove = NULL) {
  toConvert <- file(fileName, open = "r")
  nbLines <- countLines(fileName)
  if (k2rg == FALSE) {
    fileNameK2RG <- NULL
    index <- 1
    while (TRUE) {
      line = readLines(toConvert, n = 1)
      if (length(line) == 0) {
        break
      }
      if (substr(line, start = 0, stop = 1) != ">"){
        next
      }
      if (index == 1){
        isQuality <- grepl("Q", line[1])
      }
      resultLine <- .getInfoLine(line, counts, pairedEnd, order, exonicReads, isQuality)  # get all the informations for the line
      eventName <- resultLine$eventName
      variantLength <- resultLine$variantLength
      variantCounts <- resultLine$variantCounts
      if (index == 1){
        events.mat <- matrix(NA, nbLines[1] / 2, length(variantCounts) + 1)
        events.names <- rep(NA, nbLines[1] / 2)
        psiInfo <- matrix(NA, nbLines[1] / 2, length(resultLine$psiInfo))
      }
      events.mat[index, 1] <- as.numeric(variantLength)
      events.mat[index, 2:NCOL(events.mat)] <- variantCounts
      events.names[index] <- eventName
      psiInfo[index, ] <- resultLine$psiInfo
      index <- index + 1
      class(events.mat) <- "numeric"
    }
    events.df <- data.frame(events.names, events.mat)
    
  } else {
    fileNameK2RG <- fileName
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
    
    keepEvents <- wantedEvents(keep,remove)
    
    index <- 1
    iEvents <- 0  # nombre de bcc unique + duplique = nombre d'evenements 
    lEvents <- list()
    while (TRUE) {
      line = readLines(toConvert, n = 1)
      if (length(line) == 0) {
        break
      }
      if(substr(line[1], 0, 1)=="#"){
        index <- index + 1
        next
      }
      bcc <- strsplit(line, split = "\t")[[1]][EVENTNAME]
      if (strsplit(line, split = "\t")[[1]][EVENT] %in% keepEvents){
        lEvents[iEvents + 1] <- bcc
        iEvents <- iEvents + 1
      }
      index <- index + 1
    }
    lBcc <- unique(lEvents)
    iBcc <- length(lBcc)  # nombre de bcc unique
    matBccApp <- matrix(0, nrow = iBcc) # nombre d'apparition pour chaque BCC
    rownames(matBccApp) <- lBcc
    iDupBcc <- 1
    index <- 1
    indexNames <- 1
    seek(toConvert, 0) # reinitialize the cursor at the beginning of the file
    while (TRUE) {
      line = readLines(toConvert, n = 1)
      if (length(line) == 0) {
        break
      }
      if(substr(line[1], 0, 1)=="#"){
        index <- index + 1
        indexNames <- 1
        next
      }
      lLine <- strsplit(line, split = "\t")[[1]]
      if(lLine[EVENT] %in% keepEvents){
        bcc <- lLine[EVENTNAME]
        matBccApp[bcc, 1] <- matBccApp[bcc, 1] + 1
        if (matBccApp[bcc, 1] == 1) {
          resultLine <- .getInfoLineK2rg(line, counts, pairedEnd, order, exonicReads)
          if (indexNames == 1){
            events.mat <- matrix(NA, iBcc * 2, length(resultLine$variantCountsUp) + 1)
            events.names <- rep(NA, iBcc * 2)
            psiInfo <- matrix(NA, iBcc * 2, length(resultLine$variantCountsUp))
          }
          resultLine$variantLengthUp <- as.numeric(resultLine$variantLengthUp) + as.numeric(resultLine$variantLengthLow)
          events.mat[indexNames, 1] <- as.numeric(resultLine$variantLengthUp)
          events.mat[indexNames, 2:NCOL(events.mat)] <- resultLine$variantCountsUp
          events.names[indexNames] <- resultLine$eventName
          psiInfo[indexNames, ] <- resultLine$psiInfoUp
          events.mat[indexNames + 1, 1] <- as.numeric(resultLine$variantLengthLow)
          events.mat[indexNames + 1, 2:NCOL(events.mat)] <- resultLine$variantCountsLow
          events.names[indexNames + 1] <- resultLine$eventName
          psiInfo[indexNames + 1, ] <- resultLine$psiInfoLow
          class(events.mat) <- "numeric"
          indexNames <- indexNames + 2
        }
      }
      index <- index + 1
    }
    events.df <- data.frame(events.names, events.mat)
  }
  
  colnames(events.df) <- c("events.names", "events.length", paste("counts", 1:(length(colnames(events.df)) - 2), sep="")) # change col names
  
  close(toConvert)
  psiInfo <- data.frame(events.names, as.data.frame(psiInfo))
  
  output <- list(countsEvents = events.df, psiInfo = psiInfo, exonicReadsInfo = exonicReads, k2rgFile = fileNameK2RG)
  class(output) <- c("list", "countsData")
  return(output)
}

wantedEvents <- function(keep = c("All"), remove = NULL){
  EVENTS <- c("deletion", "insertion", "IR", "ES", "altA", "altD", "altAD", "alt", "unclassified", "-", " ", "", "unclassifiedSNP")
  ES_EVENTS <- c("MULTI", "alt", "altA", "altD", "altAD")
  wEvents <- c()
  if (keep == c("All") && is.null(remove)) {
    wEvents <- EVENTS
    for (i in 1:length(ES_EVENTS)) {
      wEvents <- append(wEvents, paste("ES_", ES_EVENTS[i], sep = ""))
    }
    return(wEvents)
  }
  
  if (!is.null(remove)) {
    for (i in 1:length(remove)) {
      if (!remove[i] %in% append(EVENTS, "MULTI")) {
        print(paste("In remove : couldn't find", remove[i]))
        stop("One of the element(s) of the remove vector is not part of : deletion, insertion, IR, ES, altA, altD, altAD, alt, unclassified, -, MULTI, , unclassifiedSNP")
      }
    }
  }
  ES <- FALSE
  if (keep[1] == "All") {
    for (i in 1:length(EVENTS)) {
      if (!EVENTS[i] %in% remove) {
        wEvents <- append(wEvents, EVENTS[i])
      }
    }
    if ("ES" %in% remove) {
      ES <- TRUE
    }
    if (ES == FALSE) {
      for (i in 1:length(ES_EVENTS)) {
        wEvents <- append(wEvents, paste("ES_", ES_EVENTS[i], sep = ""))
      }
    }
    return(wEvents)
  }
  for (i in 1:length(keep)) {
    if (!keep[i] %in% EVENTS) {
      print(paste("In keep : couldn't find", keep[i]))
      stop("One of the element(s) of the keep vector is not part of : deletion, insertion, IR, ES, altA, altD, altAD, alt, unclassified, -, , unclassifiedSNP")
    }
    if (ES == FALSE && keep[i] == "ES") {
      ES <- TRUE
    }
    wEvents <- append(wEvents, keep[i])
  }
  if (ES == FALSE && !is.null(remove)) {
    stop("Keep and remove can not be set together, unless keep contain ES (in that case, remove will act on ES events)")
  }
  if (ES == FALSE) {
    return(wEvents)
  }
  if (is.null(remove)) {
    for (i in 1:length(ES_EVENTS)) {
      wEvents <- append(wEvents, paste("ES_", ES_EVENTS[i], sep = ""))
    }
    return(wEvents)
  }
  for (i in 1:length(remove)){
    if (!remove[i] %in% ES_EVENTS) {
      print(paste("In remove : couldn't find",remove[i]))
      stop("One of the element(s) of the remove vector is not part of : altA, altD, altAD, alt,MULTI")
    }
  }
  for (i in 1:length(ES_EVENTS)) {
    if (!ES_EVENTS[i] %in% remove) {
      wEvents <- append(wEvents, paste("ES_", ES_EVENTS[i], sep = ""))
    }
  }
  return(wEvents)
}


print.countsData <- function(x, ...) {
  print(x$countsEvents)
}


qualityControl <- function(countsData, conditions, storeFigs = FALSE) {
  
  options(warn = -1)  # suppress the warning for the users
  
  if (storeFigs == FALSE) {
    pathToFigs <- NA
  } else {
    if (isTRUE(storeFigs)) {
      pathToFigs <- "kissDEFigures"
    } else {
      pathToFigs <- storeFigs
    }
  }
  
  # create a new folder if it doesn't exist
  if (!is.na(pathToFigs)) {
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
  ### select events with highest variance (on PSI)
  ###################################################
  countsData2 <- reshape(countsData[, c(1,(dimns):(dimns + length(conds)))], timevar = "Path", idvar = "ID", direction = "wide")
  for (i in 1:n){
    for (j in 1:nr[i]){
      countsData2$PSI <- countsData2[,(1+j+sum(nr[0:(i-1)]))]/
        (countsData2[,(1+j+sum(nr[0:(i-1)]))] + countsData2[,(1+sum(nr)+j+sum(nr[0:(i-1)]))])
      # replace PSI by NA if count less or equal to 10 reads for the 2 isoforms
      indexNA <- intersect(which(countsData2[,(1+j+sum(nr[0:(i-1)]))] < 10), which(countsData2[,(1+sum(nr)+j+sum(nr[0:(i-1)]))] < 10))
      countsData2$PSI[indexNA] <- NA
      colnames(countsData2)[9+j+(i-1)*2] <- paste(paste0("Cond",i),paste0("repl",j),sep="_")
    }
  }
  countsData2$vars <- apply(as.matrix(countsData2[, 10:(9 + length(conds))]), 1, var, na.rm = TRUE)
  ntop <- min(500, dim(countsData2)[1])
  selectntop <- order(countsData2$vars, decreasing=TRUE)[seq_len(ntop)]
  countsData2Selected <- countsData2[selectntop,]
  # remove all NAs
  countsData2Selected <- countsData2Selected[complete.cases(countsData2Selected[, 10:(9 + length(conds))]),]
  
  ###################################################
  ### code chunk number 3: heatmap
  ###################################################
  if (storeFigs == FALSE) {
    heatmap.2(as.matrix(as.dist(1 - cor(countsData2Selected[, 10:(9 + length(conds))]))), margins = c(10, 10), 
              cexRow = 1, cexCol = 1, density.info = "none", trace = "none")
    par(ask = TRUE)
  } else {
    filename <- paste(storeFigs, "/heatmap.png", sep = "")
    png(filename)
    heatmap.2(as.matrix(as.dist(1 - cor(countsData2Selected[, 10:(9 + length(conds))]))), margins = c(10, 10), 
              cexRow = 1, cexCol = 1, density.info = "none", trace = "none")
    void <- dev.off()
  }
  
  ###################################################
  ### PCA plot
  ###################################################
  pca <- prcomp(t(countsData2Selected[, 10:(9 + length(conds))]))
  fac <- factor(conds)
  colorpalette <- c("#192823", "#DD1E2F", "#EBB035", "#06A2CB", "#218559", "#D0C6B1")
  colors <- colorpalette[1:n]
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  if (storeFigs == FALSE) {
    par(oma = c(2, 1, 1, 1))
    plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main="PCA plot")
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend=levels(fac), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", pch = 20, col = colors)
  } else {
    filename <- paste(storeFigs, "/pca.png", sep = "")
    png(filename)
    par(oma = c(2, 1, 1, 1))
    plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main="PCA plot")
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend=levels(fac), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", pch = 20, col = colors)
    void <- dev.off()
  }
  
  # ###################################################
  # ### code chunk number 4: intra-group and inter-group-variance
  # ###################################################
  # # Mean and variance over all conditions and replicates (normalized counts!)
  # countsData$mn <- rowMeans(countsData[, (dimns + 1):(dimns + length(conds))])
  # countsData$var <- apply(countsData[, (dimns + 1):(dimns + length(conds))], 1, var)
  # # correction term
  # nbAll <- sum(nr)  # number of all observations in all groups
  # countsData$ct <- rowSums(countsData[, (dimns + 1):(dimns + length(conds))])^2 / nbAll
  # # sum of squares between groups
  # countsData$ss <- rowSums(countsData[, (dimns + 1):(dimns + length(conds) / n)])^2 / nr[1] + rowSums(countsData[, ((dimns + 1) + length(conds) / n):(dimns + length(conds))])^2 / nr[2]
  # # substract the correction term from the SS and divide by the degrees of 
  # df <- 1 # freedom(groups); here: df=2-1=1
  # countsData$varInter <- (countsData$ss - countsData$ct) / df
  # # intra-variability 
  # countsData$varC1 <- apply(countsData[, (dimns + 1):(dimns + nr[1])], 1, var)
  # countsData$varC2 <- apply(countsData[, ((dimns + 1) + nr[1]):(dimns + nr[2] + nr[1])], 1, var)
  # countsData$varIntra <- rowMeans(data.frame(countsData$varC1, countsData$varC2))
  # 
  # ###################################################
  # ### code chunk number 5: intra-vs-inter
  # ###################################################
  # if (storeFigs == FALSE) {
  #   plot(x = countsData$varIntra, y = countsData$varInter, xlab = "Intra-variability", ylab = "Inter-variability", las = 1, log = "xy")
  #   abline(a = 0, b = 1, col = 2, lty = 2, lwd = 2)
  # } else {
  #   filename <- paste(storeFigs, "/InterIntraVariability.png", sep = "")
  #   png(filename)
  #   plot(x = countsData$varIntra, y = countsData$varInter, xlab = "Intra-variability", ylab = "Inter-variability", las = 1, log = "xy")
  #   abline(a = 0, b = 1, col = 2, lty = 2, lwd = 2)
  #   void <- dev.off()
  # }
}



diffExpressedVariants <- function(countsData, conditions, pvalue = 1, filterLowCountsVariants = 10, flagLowCountsConditions = 10) {
  options(warn = -1)  # suppress the warning for the users
  
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
  
  if (!is.na(chunk0)) {  # no error in chunk 0
    ASSBinfo <- chunk0$ASSBinfo  # in case counts option in kissplice2counts is at 1 or 2, we have info about junction counts (ASSB), that will be useful to correct the computation of delta psi in the end. They are stored here.
    if (!is.null(ASSBinfo)) {
      li <- c()
      for (i in (1:NROW(ASSBinfo))) {
        if (i%%2 != 0) {
          li <- c(li, i)
        }
      }
      ASSBinfo <- ASSBinfo[li, ]
    }
    print("Trying to fit models on data...")
    chunk1 <- tryCatch({.modelFit(chunk0$countsData, chunk0$n, chunk0$nr, ASSBinfo, filterLowCountsVariants)
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
        sizeOfEffect <- .sizeOfEffectCalc(chunk2$signifVariants, chunk1$ASSBinfo, chunk0$n, chunk0$nr, chunk0$sortedconditions, 
                                          flagLowCountsConditions, chunk1$lengths, countsData$exonicReadsInfo)
        return(list(finalTable = sizeOfEffect$signifVariants.sorted, 
                    correctedPVal = chunk2$correctedPVal, 
                    uncorrectedPVal = chunk2$noCorrectPVal, 
                    resultFitNBglmModel = chunk1$pALLGlobalPhi.glm.nb,
                    `f/psiTable` = sizeOfEffect$psiTable,
                    k2rgFile = countsData$k2rgFile))
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


# writeMergeOutput(sizeOfEffect$signifVariants.sorted,sizeOfEffect$psiTable,output,countsData$k2rgFile)
writeOutputKissDE <- function(resDiffExprVariant, adjPvalMax = 1, dPSImin = 0, output, writePSI = FALSE) {
  if (adjPvalMax > 1 || adjPvalMax < 0) {
    print("ERROR: Invalid pvalMax (0 <= pvalMax <= 1).")
    return
  }
  if (dPSImin > 1 || dPSImin < 0) {
    print("ERROR: Invalid dPSImin (0 <= dPSImin <= 1). A dPSImin = 0.1 will catch all dPSI from -1 to -0.1 and all dPSI from 0.1 to 1.")
    return
  }
  k2rgFile <- resDiffExprVariant$k2rgFile
  if (writePSI){
    .writePSITable(resDiffExprVariant, adjPvalMax, dPSImin, output)
  } else{
    if (is.null(k2rgFile)) {
      .writeTableOutput(resDiffExprVariant$finalTable, adjPvalMax, dPSImin, output)
    }
    else {
      .writeMergeOutput(resDiffExprVariant, k2rgFile, adjPvalMax, dPSImin, output)
    }
  }
}