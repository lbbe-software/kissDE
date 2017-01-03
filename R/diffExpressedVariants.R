kissplice2counts <- function(fileName, counts = 0, pairedEnd = FALSE, order = NULL, exonicReads = TRUE, discoSNP = FALSE, k2rg = FALSE, keep = c("All"), remove = NULL) {
  toConvert <- file(fileName, open = "r")
  lines <- readLines(toConvert)
  if (k2rg == FALSE) {
    fileNameK2RG = NULL
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
    
    i <- 1
    nextI <- 1
    line <- lines[i]
    if(substr(line[1], 0, 1)=="#"){
      i <- 2
      nextI <- 2
      line <- lines[i]
    }
    resultLine <- .getInfoLineK2rg(line, counts, pairedEnd, order, exonicReads)
    variantCountsUp <- resultLine$variantCountsUp
    iEvents <- 0  # nombre de bcc unique + duplique = nombre d'evenements 
    lEvents <- list()
    while (i <= length(lines)) {
      bcc <- strsplit(line, split = "\t")[[1]][EVENTNAME]
      if (strsplit(line, split = "\t")[[1]][EVENT] %in% keepEvents){
        lEvents[iEvents + 1] <- bcc
        iEvents <- iEvents + 1
      }
      i <- i + 1
      line <- lines[i]
    }
    lBcc <- unique(lEvents)
    iBcc <- length(lBcc)  # nombre de bcc unique
    events.mat <- matrix(NA, iBcc * 2, length(variantCountsUp) + 1)
    events.names <- rep(NA, iBcc * 2)
    psiInfo <- matrix(NA, iBcc * 2, length(resultLine$psiInfoUp))
    indexNames <- 1
    
    matBccApp <- matrix(0,nrow = iBcc) # nombre d'apparition pour chaque BCC
    rownames(matBccApp) <- lBcc
    iDupBcc <- 1
    i <- nextI
    while (i <= length(lines)) {
      line <- lines[i]
      lLine <- strsplit(line, split = "\t")[[1]]
      if(lLine[EVENT] %in% keepEvents){
        bcc <- lLine[EVENTNAME]
        matBccApp[bcc, 1] <- matBccApp[bcc, 1] + 1
        if (matBccApp[bcc, 1] == 1) {
          resultLine <- .getInfoLineK2rg(line, counts, pairedEnd, order, exonicReads)
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
      i <- i + 1
      # events.df <- as.data.frame(events.mat)
      # events.df <- data.frame(events.names, events.df)
    }
    events.df <- data.frame(events.names, events.mat)
  }
  
  close(toConvert)
  psiInfo <- data.frame(events.names, as.data.frame(psiInfo))
  
  output <- list(countsEvents = events.df, psiInfo = psiInfo, discoInfo = discoSNP, exonicReadsInfo = exonicReads, k2rgFile = fileNameK2RG)
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



diffExpressedVariants <- function(countsData, conditions, storeFigs = FALSE, pvalue = 0.05, filterLowCountsVariants = 10, flagLowCountsConditions = 10, discoSNP = FALSE, output = "./kissDE_result") {
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
        if (discoSNP != FALSE) {
          if (is.na(countsData$discoInfo)) {
            discoSNP <- FALSE
          } else {
            discoSNP <- countsData$discoInfo
          }
        } 
        
        sizeOfEffect <- .sizeOfEffectCalc(chunk2$signifVariants, chunk1$ASSBinfo, chunk0$n, chunk0$nr, chunk0$sortedconditions, 
                                          flagLowCountsConditions, chunk1$lengths, discoSNP, countsData$exonicReadsInfo)
        if(!is.null(countsData$k2rgFile)) {
          print("Writing output...")
          writeMergeOutput(sizeOfEffect$signifVariants.sorted,sizeOfEffect$psiTable,output,countsData$k2rgFile)
        }
        return(list(finalTable = sizeOfEffect$signifVariants.sorted, 
                    correctedPVal = chunk2$correctedPVal, 
                    uncorrectedPVal = chunk2$noCorrectPVal, 
                    resultFitNBglmModel = chunk1$pALLGlobalPhi.glm.nb,
                    `f/psiTable` = sizeOfEffect$psiTable))
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



plotPSI <- function(diffVariants, conditions, thresholdPvalue = 0.05, thresholdDeltaPSI = 0.10, show = "PSI", storeFigs = FALSE) {
  
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
  
  ######################
  #### Prepare data ####
  ######################
  psi <- diffVariants$psiTable
  psi <- psi[order(psi$ID), ]
  signifVariants <- diffVariants$finalTable
  signifVariants <- signifVariants[order(signifVariants$ID), ]
  nbAll <- length(conditions)
  sortedconditions <- sort(conditions)
  uniqconditions <- unique(sortedconditions)
  n <- length(uniqconditions)
  nr <- rle(sortedconditions)$lengths
  psi$totCount <- apply(signifVariants[ ,(2 + 1):(2 + nbAll + nbAll)], 1, sum)
  psi$deltaPSI <- signifVariants$`Deltaf/DeltaPSI`
  psi$Adjusted_pvalue <- signifVariants$Adjusted_pvalue
  # calculate the mean of PSI by condition if show = "both" or show = "mean"
  if (show == "both" || show == "mean") {
    namePsi <- names(psi)
    meanPsi <- rowMeans(psi[, 2:(1+nr[1])], na.rm = TRUE)
    psi <- cbind(psi, meanPsi)
    nameColumn <- paste("Mean_PSI_", uniqconditions[1], sep = "")
    namePsi <- c(namePsi, nameColumn)
    for (i in 2:n) {
      meanPsi <- rowMeans(psi[, (2 + cumsum(nr)[i - 1]):(1 + cumsum(nr)[i])], na.rm = TRUE)
      psi <- cbind(psi, meanPsi)
      nameColumn <- paste("Mean_PSI_", uniqconditions[i], sep = "")
      namePsi <- c(namePsi, nameColumn)
    }
    names(psi) <- namePsi
  }
  new_data <- psi[rowSums(is.na(psi[, (2:(1 + cumsum(nr)[n]))])) == 0, ]  # remove events with only NA in PSI value
  new_data_ordered <- new_data[order(new_data$totCount), ] # re-order the data by total count
  # add a column containing the info to draw the correct shape and color in function of the group (significant or not) : significant events = 2 and not significant events = 1
  new_data_ordered$signif <- ifelse((!is.na(new_data_ordered$Adjusted_pvalue) & !is.na(new_data_ordered$deltaPSI) & new_data_ordered$Adjusted_pvalue <= thresholdPvalue & abs(new_data_ordered$deltaPSI) >= thresholdDeltaPSI), 2, 1)
  
  ##############
  #### Plot ####
  ##############
  # define the colors for the plot (max 8 conditions)
  colors <- list()
  current_colors <- palette()
  colors_transparent <- adjustcolor(palette(), alpha.f = 0.5)
  for (i in 1:n) {
    condColors <- c(colors_transparent[i], current_colors[i])
    colors[[i]] <- condColors
  }
  point_type <- c(1, 16)
  point_type_mean <- c(4, 8)
  par(xpd = TRUE)
  if (storeFigs == FALSE) {
    columnName <- paste(uniqconditions[1], "_repl", 1, sep = "")
    legendLabels <- c()
    legendCol <- c()
    legendPch <- c()
    legendNCol <- 0
    if (show == "PSI" || show == "both") {
      plot(x = c(1:length(new_data_ordered[, columnName])), y = new_data_ordered[, columnName], 
           ylab="Percent Spliced In (PSI)", col = colors[[1]][new_data_ordered$signif], pch = point_type[new_data_ordered$signif], ylim = c(0, 1), xaxt = "n", xlab = "")
      for (i in 1:n) {
        for (j in 1:nr[i]) {
          if (i != 1 || j != 1) {
            columnName <- paste(uniqconditions[i], "_repl", j, sep = "")
            points(x = c(1:length(new_data_ordered[, columnName])), y = new_data_ordered[, columnName], 
                   col = colors[[i]][new_data_ordered$signif], pch = point_type[new_data_ordered$signif])
          }
        }
        legendLabels <- c(legendLabels, uniqconditions[i], paste(uniqconditions[i], "signif", sep = " "))
        legendCol <- c(legendCol, colors[[i]])
        legendPch <- c(legendPch, point_type)
        legendNCol <- legendNCol + 1
        if (show == "both") {
          columnNameMean <- paste("Mean_PSI_", uniqconditions[i], sep = "")
          points(x = c(1:length(new_data_ordered[, columnNameMean])), y = new_data_ordered[, columnNameMean], 
                 col = colors[[i]][new_data_ordered$signif], pch = point_type_mean[new_data_ordered$signif])
          legendLabels <- c(legendLabels, paste(uniqconditions[i], "mean", sep = " "), paste(uniqconditions[i], "mean\nsignif", sep = " "))
          legendCol <- c(legendCol, colors[[i]])
          legendPch <- c(legendPch, point_type_mean)
          legendNCol <- legendNCol + 1
        }
      }
      legend(10, -0.07, legendLabels, col = legendCol, pch = legendPch, cex = 0.75, ncol = legendNCol)
    } else if (show == "mean") {
      for (i in 1:n) {
        columnNameMean <- paste("Mean_PSI_", uniqconditions[i], sep = "")
        if (i == 1) {
          plot(x = c(1:length(new_data_ordered[, columnNameMean])), y = new_data_ordered[,columnNameMean], ylab = "Percent Spliced In (PSI)", 
               col = colors[[i]][new_data_ordered$signif], pch = point_type_mean[new_data_ordered$signif], ylim = c(0, 1), xaxt = "n", xlab = "")
        } else {
          points(x = c(1:length(new_data_ordered[, columnNameMean])), y = new_data_ordered[, columnNameMean], 
                 col = colors[[i]][new_data_ordered$signif], pch = point_type_mean[new_data_ordered$signif])
        }
        legendLabels <- c(legendLabels, paste(uniqconditions[i], "mean", sep = " "), paste(uniqconditions[i], "mean\nsignif", sep = " "))
        legendCol <- c(legendCol, colors[[i]])
        legendPch <- c(legendPch, point_type_mean)
        legendNCol <- legendNCol + 1
      }
      legend(10, -0.07, legendLabels, col = legendCol, pch = legendPch, cex = 0.75, ncol = legendNCol)
    }
  } else {
    filename <- paste(storeFigs, "/psi.png", sep = "")
    widthPNG <- as.integer(length(new_data_ordered$ID) * 4.8)
    png(filename, width = widthPNG)
    columnName <- paste(uniqconditions[1], "_repl", 1, sep = "")
    legendLabels <- c()
    legendCol <- c()
    legendPch <- c()
    legendNCol <- 0
    if (show == "PSI" || show == "both") {
      plot(x = c(1:length(new_data_ordered[, columnName])), y = new_data_ordered[, columnName], ylab = "Percent Spliced In (PSI)", 
           col = colors[[1]][new_data_ordered$signif], pch = point_type[new_data_ordered$signif], ylim = c(0, 1), xaxt = "n", xlab = "")
      for (i in 1:n) {
        for (j in 1:nr[i]) {
          if (i != 1 || j != 1) {
            columnName <- paste(uniqconditions[i], "_repl", j, sep = "")
            points(x = c(1:length(new_data_ordered[, columnName])), y = new_data_ordered[, columnName], col = colors[[i]][new_data_ordered$signif], 
                   pch = point_type[new_data_ordered$signif])
          }
        }
        legendLabels <- c(legendLabels, uniqconditions[i], paste(uniqconditions[i], "signif", sep = " "))
        legendCol <- c(legendCol, colors[[i]])
        legendPch <- c(legendPch, point_type)
        legendNCol <- legendNCol + 1
        if (show == "both") {
          columnNameMean <- paste("Mean_PSI_", uniqconditions[i], sep = "")
          points(x = c(1:length(new_data_ordered[, columnNameMean])), y = new_data_ordered[, columnNameMean], 
                 col = colors[[i]][new_data_ordered$signif], pch = point_type_mean[new_data_ordered$signif])
          legendLabels <- c(legendLabels, paste(uniqconditions[i], "mean", sep = " "), paste(uniqconditions[i], "mean\nsignif", sep = " "))
          legendCol <- c(legendCol, colors[[i]])
          legendPch <- c(legendPch, point_type_mean)
          legendNCol <- legendNCol + 1
        }
      }
      legend(10, -0.07, legendLabels, col = legendCol, pch = legendPch, cex = 0.75, ncol = legendNCol)
    } else if (show == "mean") {
      for (i in 1:n) {
        columnNameMean <- paste("Mean_PSI_", uniqconditions[i], sep = "")
        if (i == 1) {
          plot(x = c(1:length(new_data_ordered[, columnNameMean])), y = new_data_ordered[, columnNameMean], ylab = "Percent Spliced In (PSI)", 
               col = colors[[i]][new_data_ordered$signif], pch = point_type_mean[new_data_ordered$signif], ylim = c(0, 1), xaxt = "n", xlab = "")
        } else {
          points(x = c(1:length(new_data_ordered[, columnNameMean])), y = new_data_ordered[, columnNameMean], 
                 col = colors[[i]][new_data_ordered$signif], pch = point_type_mean[new_data_ordered$signif])
        }
        legendLabels <- c(legendLabels, paste(uniqconditions[i], "mean", sep = " "), paste(uniqconditions[i], "mean\nsignif", sep = " "))
        legendCol <- c(legendCol, colors[[i]])
        legendPch <- c(legendPch, point_type_mean)
        legendNCol <- legendNCol + 1
      }
      legend(10, -0.07, legendLabels, col = legendCol, pch = legendPch, cex = 0.75, ncol = legendNCol)
    }
    void <- dev.off()
  }
}

writeMergeOutput <- function(finalTable,psiTable,output,k2rgFile) {
  EVENTNAME <- 16
  COUNTSSTART <- 3
  COUNTSENDBEFOREEND <- 3
  PSISTART <- 2
  
  COUNTSEND <- ncol(finalTable)-COUNTSENDBEFOREEND
  PSIEND <- ncol(psiTable)
  
  lBcc <- rownames(finalTable)
  fK2RG <- file(k2rgFile, open = "r")
  lines <- readLines(fK2RG)
  fOut <- file(output, open = "w")
  i <- 1
  line <- lines[i]
  if(substr(line[1],0,1)=="#"){
    nCol <- length(strsplit(line, split = "\t")[[1]])
    countsLabel <- paste(colnames(finalTable)[COUNTSSTART:COUNTSEND],collapse=",")
    countsName <- paste("CountsNorm(",countsLabel,")",sep="")
    psiLabel <- paste(colnames(psiTable)[PSISTART:PSIEND],collapse=",")
    psiName <- paste("psiNorm(",psiLabel,")",sep="")
    countsHead <- paste(nCol+1,countsName,sep=":")
    psiHead <- paste(nCol+2,psiName,sep=":")
    pvHead <- paste(nCol+3,"adjusted_pvalue",sep=":")
    dPSIHead <- paste(nCol+4,"dPSI",sep=":")
    warnHead <- paste(nCol+5,"warnings",sep=":")
    toWrite <- paste(line,countsHead,psiHead,pvHead,dPSIHead,warnHead,sep="\t")
    writeLines(toWrite,fOut)
    i <- 2
  }
  while(i <= length(lines)) {
    line <- lines[i]
    bcc <- strsplit(line, split = "\t")[[1]][EVENTNAME]
    bcc <- strsplit(bcc, split = "Type_")[[1]][1]
    bcc <- substr(bcc,0,nchar(bcc)-1)
    if(bcc %in% lBcc) {
      countsValue <- paste(as.character(finalTable[rownames(finalTable)==bcc,][COUNTSSTART:COUNTSEND]),collapse=",")
      psiValue <- as.numeric(psiTable[psiTable$ID==bcc,][PSISTART:PSIEND])
      psiValue <- paste(as.character(psiValue),collapse=",")
      toWrite <- paste(line,countsValue,psiValue,finalTable[rownames(finalTable)==bcc,]$Adjusted_pvalue,finalTable[rownames(finalTable)==bcc,]$`Deltaf/DeltaPSI`,finalTable[rownames(finalTable)==bcc,]$lowcounts,sep="\t")
      writeLines(toWrite,fOut)
    }
    i <- i + 1
  }
  close(fK2RG)
  close(fOut)
}