.lineParse <- function(line, counts=0) {
    if (counts == 0) {
        splittedCounts <- grep(pattern="C[[:digit:]]+_[[:digit:]]+", x=line, 
                                perl=TRUE, value=TRUE)
    }
    else {  # no C in the header => ASSB counts
        splittedCounts <- grep(pattern="[ASB]{1,4}[[:digit:]]+_[[:digit:]]+", 
                            x=line, perl=TRUE, value=TRUE)
    }
    ## gets the junctions id (ex 1 in AS1) and the count (ex 6 in AS1_6)
    split1 <- sub("(?=[0-9])(?<=[A-Z])", "_", splittedCounts, perl=TRUE)
    split2 <- strsplit(split1, "_", fixed=TRUE)
    countsperCond <- matrix(unlist(split2), nrow=3)
    
    return(countsperCond)
}



.getJunctionCounts <- function(X){
    return(unlist(regmatches(x=X[1], 
        m=gregexpr(pattern="[0-9]+", text=X[1]))))
}



.testAggregate <- function(nbVec, countsVec){
    nbcond <- max(nbVec)
    agg <- rep_len(0,nbcond)
    for (i in seq_along(countsVec))
        agg[nbVec[i]] = agg[nbVec[i]] + countsVec[i]
    
    return(agg)
}



.countsSetk2rg <- function(splittedCounts, counts=0, pairedEnd=FALSE, 
        order=NULL, exonicReads=TRUE) {
    split1 <- sub("(?=[0-9])(?<=[A-Z])", "_", splittedCounts[[1]], perl=TRUE)
    split2 <- strsplit(split1, "_", fixed=TRUE)
    countsperCond <- matrix(unlist(split2), nrow=3)
    ## Create the vectors of counts and samples
    namesVec <- countsperCond[1,]
    nbVec <- as.numeric(countsperCond[2,])
    countsVec <- as.numeric(countsperCond[3,])
    psiVec <- rep.int(0, length(countsVec))
    if (counts >= 1) {
        
        ## Replace the maximum of AS or SB by the minimum of AS or SB
        asIndex <- namesVec == "AS"
        sbIndex <- namesVec == "SB"
        if(countsVec[asIndex]>=countsVec[sbIndex]) {
            countsVec[asIndex]=countsVec[sbIndex]
            
        } else {
            countsVec[sbIndex]=countsVec[asIndex]
        }
        ## Replace all ASSB counts by their opposite and fill psiVec
        assbIndex <- namesVec == "ASSB"
        countsVec[assbIndex] <- -countsVec[assbIndex]
        psiVec[assbIndex] <- countsVec[assbIndex]
        
        if ((counts == 2) & (exonicReads == FALSE)) {
            ## Replace all S counts by 0
            sIndex <- namesVec == "S"
            countsVec[sIndex] <- 0
            
        }
    }
    if (counts >= 1) {
        ## sums the counts for each junction that belongs to the same event
        sums <- .testAggregate(nbVec, countsVec)
        ## dpsi will store counts of ASSB counts 
        assbPsi <- .testAggregate(nbVec, psiVec)
        
        if (pairedEnd == TRUE) {
            if (is.null(order)) {
                nr <- length(sums) / 2
                order <- rep(x = seq_len(nr), times = rep(2, nr))
            } else {
                if (!is.vector(order)) {
                    stop("Order vector seems to be in a wrong format.")
                }
            }
            ## in case data is paired-end, there is one more sum to do, 
            ## for each part of the pair
            sums <- .testAggregate(order, sums)
            assbPsi <- .testAggregate(order, assbPsi)
        } 
    } else { ## counts == 0 
        if (pairedEnd == TRUE) {
            if (is.null(order)) {
                ## for length(s)=8, will create a vector c(1,1,2,2,3,3,4,4) 
                ## (assuming data is ordered)
                nr <- length(countsVec) / 2
                order <- rep(x = seq_len(nr), times = rep(2, nr))
            } else {
                if (!is.vector(order)) {
                    stop("'order' vector seems to be in a wrong format.")
                }
            }
        } else {
            order <- c(seq_along(countsVec))
        }
        myby <- list(order)
        sums <- .testAggregate(order, countsVec)
        assbPsi <- NULL
    }
    listCounts <- sums
    
    return(list(vCounts=listCounts, 
                psiCounts=assbPsi))
}



.countsSet <- function(line, counts=0, pairedEnd=FALSE, 
        order=NULL, exonicReads=TRUE) {
    
    countsperCond <- .lineParse(line, counts)
    ## Create the vectors of counts and samples
    namesVec <- countsperCond[1,]
    nbVec <- as.numeric(countsperCond[2,])
    countsVec <- as.numeric(countsperCond[3,])
    psiVec <- rep.int(0, length(countsVec))
    if (counts >= 1) {
        ## Replace the maximum of AS or SB by the minimum of AS or SB
        asIndex <- namesVec == "AS"
        sbIndex <- namesVec == "SB"
        if(countsVec[asIndex]>=countsVec[sbIndex]) {
            countsVec[asIndex]=countsVec[sbIndex]
            
        } else {
            countsVec[sbIndex]=countsVec[asIndex]
        }
        ## Replace all ASSB counts by their opposite and fill psiVec
        assbIndex <- namesVec == "ASSB"
        countsVec[assbIndex] <- -countsVec[assbIndex]
        psiVec[assbIndex] <- countsVec[assbIndex]
        
        if ((counts == 2) & (exonicReads == FALSE)) {
            ## Replace all S counts by 0
            sIndex <- namesVec == "S"
            countsVec[sIndex] <- 0
            
        }
    }
    if (counts >= 1) {
        ## sums the counts for each junction that belongs to the same event
        sums <- .testAggregate(nbVec, countsVec)
        ## dpsi will store counts of ASSB counts 
        assbPsi <- .testAggregate(nbVec, psiVec)
        
        if (pairedEnd == TRUE) {
            if (is.null(order)) {
                nr <- length(sums) / 2
                order <- rep(x = seq_len(nr), times = rep(2, nr))
            } else {
                if (!is.vector(order)) {
                    stop("Order vector seems to be in a wrong format.")
                }
            }
            ## in case data is paired-end, there is one more sum to do, 
            ## for each part of the pair
            sums <- .testAggregate(order, sums)
            assbPsi <- .testAggregate(order, assbPsi)
        } 
    } else { ## counts == 0 
        if (pairedEnd == TRUE) {
            if (is.null(order)) {
                ## for length(s)=8, will create a vector c(1,1,2,2,3,3,4,4) 
                ## (assuming data is ordered)
                nr <- length(countsVec) / 2
                order <- rep(x = seq_len(nr), times = rep(2, nr))
            } else {
                if (!is.vector(order)) {
                    stop("'order' vector seems to be in a wrong format.")
                }
            }
        } else {
            order <- c(seq_along(countsVec))
        }
        myby <- list(order)
        sums <- .testAggregate(order, countsVec)
        assbPsi <- NULL
    }
    listCounts <- sums
    
    return(list(vCounts=listCounts, 
                psiCounts=assbPsi))
}



.getInfoLine <- function(line, counts=0, pairedEnd=FALSE, order=NULL, 
        exonicReads=TRUE, indexStart=5) {
    beginningLine <- line[seq_len(indexStart)]
    eventName <- paste(beginningLine[1], beginningLine[2], sep="|")
    variantLength <- as.numeric(strsplit(x = beginningLine[4], 
        split = "_", fixed=TRUE)[[1]][4])
    endLine <- line[indexStart:length(line)]
    resultCountsSet <- .countsSet(endLine, counts, pairedEnd, order, 
        exonicReads)
    
    return(list(eventName=eventName, 
                variantLength=variantLength, 
                variantCounts=resultCountsSet$vCounts, 
                psiInfo=resultCountsSet$psiCounts))
}



.getInfoLineK2rg <- function(line, counts=0, pairedEnd=FALSE, order=NULL, 
        exonicReads=TRUE) {
    firstPart <- line[16]
    firstPartSplit <- strsplit(x = firstPart, split = "|", fixed = TRUE)[[1]]
    eventName <- paste(firstPartSplit[1], firstPartSplit[2], sep="|")
    if(grepl(",", line[20], fixed = TRUE)) { # old k2rg file version
      countsUp <- strsplit(x = line[20], split = ",", fixed = TRUE)
      countsLow <- strsplit(x = line[21], split = ",", fixed = TRUE)
    } else { # new k2rg file version
      countsUp <- strsplit(x = line[20], split = "|", fixed = TRUE)
      countsLow <- strsplit(x = line[21], split = "|", fixed = TRUE)
    }
    resultCountsSetUp <- .countsSetk2rg(countsUp, counts, pairedEnd, 
                                        order, exonicReads)
    resultCountsSetLow <- .countsSetk2rg(countsLow, counts, pairedEnd, 
        order, exonicReads)
    variantLengthUp <- sum(as.numeric(strsplit(x = line[11], split = ",", 
        fixed = TRUE)[[1]]))
    variantLengthLow <- sum(as.numeric(strsplit(x = line[17], split = ",", 
        fixed = TRUE)[[1]]))
    
    return(list(eventName=eventName, 
                variantLengthUp=variantLengthUp, 
                variantLengthLow=variantLengthLow, 
                variantCountsUp=resultCountsSetUp$vCounts, 
                variantCountsLow=resultCountsSetLow$vCounts, 
                psiInfoUp=resultCountsSetUp$psiCounts, 
                psiInfoLow=resultCountsSetLow$psiCounts))
}



.readAndPrepareData <- function(countsData, conditions) {
    ###################################################
    ### code chunk number 1: Read data
    ###################################################
    if (is.null(countsData$psiInfo)){
        countsEvents <- countsData  ## count table provided by the user
        psiInfo <- NULL
    } else {
        countsEvents <- countsData$countsEvents  ## provided by kissplice2counts
        if (NCOL(countsData$psiInfo) > 1) {
            psiInfo <- countsData$psiInfo  ## info about ASSB counts
        } else {
            psiInfo <- NULL
        }
    }
    
    # Check options
    countsID <- countsEvents[[1]]
    i <- 1
    for(thisID in countsID) {
        if(i%%2) {
            savedID <- thisID
        }
        else{
            if(thisID!=savedID) {
                stop("Input error: 'countsData' must contains two following 
    rows for the same event name (first column). The row with the ID ", savedID,
    " is alone. See the vignette for more informations.")
            }
        }
        i <- i + 1
    }
    if(!i%%2) {
        stop("Input error: 'countsData' must contains two following rows for 
    the same event name (first column). The last row (ID: ", savedID,
    ") is alone. See the vignette for more informations.")
    }
    
    tableID <- table(countsEvents[[1]])
    if(max(tableID) != 2){
        stop("Input error: in 'countsData', the event(s) ",
    names(tableID[tableID==max(tableID)]), " has(have) more than two lines.")
    }
    
    countsLength <- countsEvents[[2]]
    if(!is.vector(countsLength,mode="numeric")){
        stop("Input error: in 'countsData', the second column ('length') is 
    not composed of numerical values.")
    }
    if(min(countsLength)<0) {
        stop("Input error: in 'countsData', the smallest value for the 
    second column (length) should be 0 (the minimal value in the current 
    data is ", min(countsLength), ").")
    }
    
    if(!is.vector(conditions)){
        stop("Input error: 'condition' must be a vector.")
    }
    if(ncol(countsEvents)-2!=length(conditions)){
        stop("Input error: not the same amount of conditions in 'countsData' 
    (starting at column 3) and 'condition'.")
    }
    
    if("*"%in%conditions){
        toRm <- which(conditions=="*")
        countsEvents <- countsEvents[, -(toRm+2)]
        conditions <- conditions[-toRm]
        if(!is.null(psiInfo)) {
            psiInfo <- psiInfo[, -(toRm+1)]
        }
    }
    
    sortedconditions <- sort(conditions)
    n <- length(unique(sortedconditions))
    nr <- rle(sortedconditions)$lengths
    ## if at least 1 condition does not contain replicates, stop the analysis
    if (length(nr[nr == 1]) > 0){
        stop("The data does not contain replicates for all conditions. 
    Please change the input data and relaunch.")
    }
    sortedindex <- order(conditions) + 2
    namesData <- c("ID", "Length", rep(NA, length(conditions)))
    for (k in seq_len(nr[1])) {
        namesData[2 + k] <- paste(sortedconditions[k], "_repl", 
                                    k, sep="", collapse="")
    }
    
    for (i in 2:n) {
        for (j in seq_len(nr[i])) {
            namesData[2 + cumsum(nr)[i - 1] + j] <- 
                paste(sortedconditions[cumsum(nr)[i - 1] + j], "_repl", j, 
                    sep="", collapse="")
        }
    }  ## proper names for conditionsXreplicates
    countsEvents[, -(seq_len(2))] <- countsEvents[, sortedindex]
    colnames(countsEvents) <- namesData
    
    if (!is.null(psiInfo)){
        psiInfo[, -1] <- psiInfo[, sortedindex - 1]
        colnames(psiInfo) <- c("events.names", namesData[c(-1, -2)])
        ASSBinfo <- data.frame(psiInfo)
        colnames(ASSBinfo) <- c("events.names", namesData[c(-1, -2)])
        colnames(psiInfo) <- namesData[c(-1, -2)]
    } else {
        ASSBinfo <- NULL
    }
    countsEvents$Path <- gl(2, 1, NROW(countsEvents), labels=c("UP", "LP"))
    
    ###################################################
    ### code chunk number 2: Normalization
    ###################################################
    # Normalization with DESeq2
    conds <- c()
    for (i in seq_len(n)) {
        conds <- c(conds, rep(paste("Cond", i, sep="", collapse=""), nr[i]))
    } 
    # Creating every 2 rows sum table :
    ## This was the object previously used
    countsEventsCounts <- countsEvents[, !(names(countsEvents) %in% 
        c("ID", "Length", "Path"))]
    
    ## This is the object that we use now
    countsEventsSum <- rowsum(countsEventsCounts, 
        as.integer(gl(nrow(countsEventsCounts), 2, nrow(countsEventsCounts))))
    
    ## create a DESeqDataSet object
    suppressMessages(dds <- DESeqDataSetFromMatrix(
        countData=countsEventsSum, 
        colData=data.frame(condition=as.factor(conds)),
        design=~ condition))
    
    ddsSF <- estimateSizeFactors(dds)
    sizeFactorsSum <- sizeFactors(ddsSF)
    suppressMessages(ddsSF <- DESeqDataSetFromMatrix(
        countData=countsEventsCounts, 
        colData=data.frame(condition=as.factor(conds)),
        design=~ condition))
    sizeFactors(ddsSF) <- sizeFactorsSum
    shouldWeNormalize <- sum(is.na(sizeFactors(ddsSF))) < 1
    dimns <- NCOL(countsEvents)
    
    # add columns containing normalized data
    
    countsEvents[, (dimns + 1):(dimns + length(conds))] <- 
        round(counts(ddsSF, normalized=shouldWeNormalize))
    colnames(countsEvents)[(dimns + 1):(dimns + length(conds))] <- 
        paste(namesData[3:(3 + sum(nr) - 1)], "_Norm", sep="")
    
    return(list(countsData=countsEvents, 
                conditions=conds,
                dim=dimns, 
                n=n, 
                nr=nr, 
                sortedconditions=sortedconditions, 
                ASSBinfo=ASSBinfo))
}



.eventtable <- function(df, startPosColumn4Counts, endPosCol4Counts) {
    eventTab <- data.frame(
        ID=rep(as.factor(df["ID"]), endPosCol4Counts-startPosColumn4Counts+1),
        cond=as.factor(unlist(lapply(
            strsplit(x = names(df)[startPosColumn4Counts:endPosCol4Counts], 
                split = "_", fixed = TRUE), 
            FUN=function(d){paste(d[2:(length(d) - 2)], collapse="_")}), 
            use.names=FALSE)), 
        ## to manage the conditions names which contains "_"
        counts=as.numeric(df[startPosColumn4Counts:endPosCol4Counts]),
        path=as.factor(unlist(lapply(
            strsplit(x = names(df)[startPosColumn4Counts:endPosCol4Counts], 
                split = "|"), 
            FUN=function(d){d[1]}), 
            use.names=FALSE)),
        row.names=NULL)
    
    return(eventTab)
}



.addOneCount <- function(df){
    df$counts <- unlist(lapply(df[, "counts"], function(x){x + 1}), 
                        use.names=FALSE)
    return(df)
}



.fitNBglmModelsDSSPhi <- function(eventdata, phiDSS, nbAll){
    
    ## binomial negative model, with phi DSS
    nbglmA <- negbin(counts~cond + path, data=eventdata, random=~1, 
                fixpar=list(4, phiDSS))
    nbglmI <- negbin(counts~cond * path, data=eventdata, random=~1, 
                fixpar=list(5, phiDSS))
    
    nbAnov <- anova(nbglmA, nbglmI)
    nbAIC <- c(AIC(nbglmA, k=log(nbAll))@istats$AIC, 
                AIC(nbglmI, k=log(nbAll))@istats$AIC)
    nbSingHes <- c(nbglmA@singular.hessian, nbglmI@singular.hessian)
    nbCode <- c(nbglmA@code, nbglmI@code)
    
    rslts <- c(nbAnov@anova.table$'P(> Chi2)'[2],
                nbAIC,
                nbCode,
                nbSingHes)
    
    return(rslts)  
}



.modelFit <-function(countsData, n, nr, ASSBinfo, filterLowCountsVariants, 
                    techRep, nbCore){
    ##################################################
    ## code chunk number 1: event-list
    ##################################################
    # reduce data frame to the interesting columns
    nbAll <- sum(nr)
    dataPart <- countsData[, c(seq_len(2), 
                        which(grepl("_Norm", names(countsData), fixed = TRUE)))]
    dataPart$Path <- gl(2, 1, NROW(countsData), labels=c("UP","LP"))
    dataPart2 <- cbind(dataPart[seq(1, NROW(dataPart), 2), ], 
                    dataPart[seq(2, NROW(dataPart), 2), 
                    grepl("Norm", names(dataPart), fixed = TRUE)])
    names(dataPart2)[3:(3 + nbAll - 1)] <- 
        paste("UP", names(dataPart2)[3:(3 + nbAll - 1)], sep="_")
    names(dataPart2)[(3 + nbAll + 1):(3 + 2 * nbAll + 1 - 1)] <- 
        paste("LP", names(dataPart2)[(3 + nbAll + 1):(3 + 2 * nbAll + 1 - 1)], 
            sep="_")
    lengths <- data.frame(dataPart[seq(1, NROW(dataPart), 2), 2], 
                    dataPart[seq(2, NROW(dataPart), 2), 2])
    colnames(lengths) <- c("upper", "lower")
    
    ## computes the difference of length between the lower and upper paths 
    dataPart2[2] <- lengths$upper - lengths$lower
    names(dataPart2)[2] <- "Length_diff"
    dataPart2 <- dataPart2[, -which(colnames(dataPart2) == "Path")]
    if (anyDuplicated(dataPart2$ID) > 0) {
        dataPart2 <- dataPart2[!duplicated(as.character(dataPart2$ID)), ]
    }
    rownames(dataPart2) <- as.character(dataPart2$ID)
    if (!is.null(ASSBinfo)) {
        rownames(ASSBinfo) <- dataPart2$ID
    }
    rownames(lengths) <- rownames(dataPart2) 
    
    ###################################################
    ### code chunk number 2: DSS dispersion estimation
    ###################################################
    ## the counts matrix
    dataNormCountsEvent <- as.matrix(dataPart2[, 3:ncol(dataPart2)])
    colnames(dataNormCountsEvent) <- seq_len(ncol(dataNormCountsEvent))
    
    # We select the variant with the largest number of reads
    dataNormCountsEventSum <- matrix(nrow=nrow(dataNormCountsEvent), ncol=nbAll)
    rownames(dataNormCountsEventSum) <- rownames(dataNormCountsEvent)
    colnames(dataNormCountsEventSum) <- seq_len(nbAll)
    for(i in seq_len(nrow(dataNormCountsEvent))){
        current <- dataNormCountsEvent[i,]
        dataNormCountsEventSum[i,] <- current[seq_len(sum(nr))]+
            current[(sum(nr)+1):(2*sum(nr))]
    }
    
    designs <- rep(0:(n-1),nr)
    dispData <- newSeqCountSet(dataNormCountsEventSum, as.data.frame(designs), 
                        normalizationFactor = rep(1, nbAll))
    ## fix the seed to avoid the stochastic outputs of the 
    ## DSS:estDispersion function
    set.seed(40)
    dispData <- estDispersion(dispData)
    names(exprs(dispData)) <- rownames(dataPart2)
    
    ###################################################
    ### code chunk number 5: exclude low counts
    ###################################################
    ## global counts for each variant (low/up) by event
    totLOW <- as.vector(apply(dataPart2[, (3 + sum(nr)):(3 + 2 * sum(nr) - 1)], 
                        1, sum))
    totUP <- as.vector(apply(dataPart2[, 3:(3 + sum(nr) - 1)], 1, sum))
    names(totLOW) <- rownames(dataPart2)
    names(totUP) <- rownames(dataPart2)
    
    if (length(-which(totUP < filterLowCountsVariants & 
                totLOW < filterLowCountsVariants)) > 0){
        
        dataPart3 <- dataPart2[-which(totUP < filterLowCountsVariants & 
                                        totLOW < filterLowCountsVariants), ]
        exprsDataWithoutLowCounts <- 
            exprs(dispData)[-which(totUP < filterLowCountsVariants & 
                totLOW < filterLowCountsVariants), ]
        
        dispersionDataWithoutLowCounts <- 
            dispersion(dispData)[-which(totUP < filterLowCountsVariants & 
                totLOW < filterLowCountsVariants)]
        
        normFactorWithoutLowCounts <- dispData@normalizationFactor
        dispData <- newSeqCountSet(exprsDataWithoutLowCounts, 
                                    as.data.frame(designs))
        dispData@normalizationFactor <- normFactorWithoutLowCounts
        dispData@dispersion <- dispersionDataWithoutLowCounts
        
    } else {
        dataPart3 <- dataPart2
    }
    
    allEventtables <- apply(dataPart3, 1, 
                            .eventtable, 
                            startPosColumn4Counts=which(
                                grepl("UP", names(dataPart3), fixed = TRUE))[1],
                            endPosCol4Counts=ncol(dataPart3))
    
    ###################################################
    ### code chunk number 6: pALLGlobalPhi.glm.nb
    ###################################################
    pALLGlobalPhiGlmNb <- data.frame(t(rep(NA, 7)))
    ## create the cluster
    cl <- parallel::makeCluster(nbCore)
    doParallel::registerDoParallel(cl)
    if(techRep) {
        pALLGlobalPhiGlmNb_list <- foreach(i=seq_along(allEventtables)) %dopar% 
            .fitNBglmModelsDSSPhi(allEventtables[[i]], 0, nbAll)
        pALLGlobalPhiGlmNb <- do.call(rbind.data.frame, pALLGlobalPhiGlmNb_list)
    }
    else {
        pALLGlobalPhiGlmNb_list <- foreach(i=seq_along(allEventtables)) %dopar% 
            .fitNBglmModelsDSSPhi(allEventtables[[i]], dispersion(dispData)[i],
                nbAll)
        pALLGlobalPhiGlmNb <- do.call(rbind.data.frame, pALLGlobalPhiGlmNb_list)
    }
    
    ## close the cluster
    parallel::stopCluster(cl)
    
    ###################################################
    ### code chunk number 7: excl_errors
    ###################################################
    sing.events <- which(grepl("Error", pALLGlobalPhiGlmNb[, 1], fixed = TRUE))
    if (length(sing.events) != 0) {
        pALLGlobalPhiGlmNb <- pALLGlobalPhiGlmNb[-sing.events, ]
    }
    colnames(pALLGlobalPhiGlmNb) <- c("I vs A", "bicA", "bicI", "codeA", 
                                        "codeI", "shA", "shI")
    if (length(sing.events) != 0) {
        rownames(pALLGlobalPhiGlmNb) <- dataPart3[-sing.events, 1]
    } else {
        rownames(pALLGlobalPhiGlmNb) <- dataPart3[, 1]
    }
    return(list(pALLGlobalPhi.glm.nb=pALLGlobalPhiGlmNb, 
                sing.events=sing.events, 
                dataPart3=dataPart3, 
                ASSBinfo=ASSBinfo, 
                allEventtables=allEventtables, 
                lengths=lengths, 
                dispData=dispData))
}



.bestModelandSingular <- function(pALLGlobalPhi.glm.nb, sing.events, dataPart3, 
                            allEventtables, pvalue, nr, dispData) {
    nbAll <- sum(nr)
    pALLGlobalPhi.glm.nb <- 
        pALLGlobalPhi.glm.nb[!is.na(pALLGlobalPhi.glm.nb[, 1]), ]
    
    if (NROW(pALLGlobalPhi.glm.nb) > 0) {
        matrixpALLGlobalPhi <- as.matrix(pALLGlobalPhi.glm.nb)
        storage.mode(matrixpALLGlobalPhi) <- "numeric"
        
        ###################################################
        ###  code chunk number 4: Pseudo-counts
        ###################################################
        singhes <- which(rowSums(matrixpALLGlobalPhi[, c(6, 7)]) != 0)
        
        singhes_n <- names(singhes)
        pALLGlobalPhi.glm.nb.pen <- as.data.frame(matrixpALLGlobalPhi)
        for (i in singhes_n) {
            pALLGlobalPhi.glm.nb.pen[i, ] <- 
                try(.fitNBglmModelsDSSPhi(.addOneCount(allEventtables[[i]]),
                    dispersion(dispData)[i], nbAll), 
                    silent=TRUE)
        }
        
        pALLGlobalPhi.glm.nb <- as.data.frame(matrixpALLGlobalPhi)
        
        pALLGlobalPhi.glm.nb$final.pval.a.ia <- 1
        
        ## negative binomial model, phi estimated with DSS
        li <- which(rowSums(pALLGlobalPhi.glm.nb[, c(6, 7)]) == 0)
        pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li, 1]
        li <- which(rowSums(pALLGlobalPhi.glm.nb[, c(6, 7)]) != 0)
        pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- 
            pALLGlobalPhi.glm.nb.pen[li, 1]
        
        sing.events.final <- which(grepl("Error", pALLGlobalPhi.glm.nb[, 8]))
        if (length(sing.events.final) != 0) {
            pALLGlobalPhi.glm.nb <- pALLGlobalPhi.glm.nb[-sing.events.final, ]
        }
        noCorrectPVal <- pALLGlobalPhi.glm.nb$final.pval.a.ia
        names(noCorrectPVal) <- rownames(pALLGlobalPhi.glm.nb)
        pALLGlobalPhi.glm.nb$final.padj.a.ia <- 
            p.adjust(pALLGlobalPhi.glm.nb$final.pval.a.ia, method="fdr")
        correctedPVal <- pALLGlobalPhi.glm.nb$final.padj.a.ia
        names(correctedPVal) <- rownames(pALLGlobalPhi.glm.nb)
        
        if (length(sing.events) != 0) {
            tmpdataPart3_1 <- dataPart3[-sing.events, ]
        } else {
            tmpdataPart3_1 <- dataPart3
        }
        
        if (length(sing.events.final) != 0) {
            tmpdataPart3 <- tmpdataPart3_1[-sing.events.final, ]
            signifVariants <- 
                cbind(tmpdataPart3, pALLGlobalPhi.glm.nb$final.padj.a.ia)
            signifVariants <-
                signifVariants[pALLGlobalPhi.glm.nb$final.padj.a.ia <= pvalue, ]
        } else {
            signifVariants <- 
                cbind(dataPart3, pALLGlobalPhi.glm.nb$final.padj.a.ia)
            signifVariants <- 
                signifVariants[pALLGlobalPhi.glm.nb$final.padj.a.ia <= pvalue, ]
        }
        
        return(list(noCorrectPVal=noCorrectPVal, 
                    correctedPVal=correctedPVal, 
                    signifVariants=signifVariants))
    } else {
        return(NA)
    }
}



.sizeOfEffectCalc <- function(signifVariants, ASSBinfo, n, nr, sortedconditions,
                        flagLowCountsConditions, lengths, exonicReads) {
    ###################################################
    ### code chunk 1: compute delta PSI/f
    ###################################################
    if (!is.null(ASSBinfo)) {
        ## select only the lines corresponding to the remaining lines of 
        ## signifVariants
        ASSBinfo <- subset(ASSBinfo, ASSBinfo$events.names %in% 
                            as.vector(signifVariants[, 1]))
    }
    
    ## to check later for low counts to flag
    sumLowCond <- matrix(data=rep(0, n * NROW(signifVariants)), 
                        nrow=NROW(signifVariants), ncol=n)
    pairsCond <- list()
    namesCond <- unique(sortedconditions)
    for (i in seq_len(n)) {
        j <- i + 1
        while (j <= n) {
            ## creation of permutation of size 2 in c
            pairsCond[[length(pairsCond) + 1]] <- list(c(i, j))
            j <- j + 1
        }
    }
    
    namesPsiPairCond <- c()
    if (!is.null(ASSBinfo)) {
        if (exonicReads){
            rown <- row.names(signifVariants)
            lengths2 <- lengths[rown, ]
        }
        else{
            ## to put the lines of the 2 data frames in the same order
            newindex  <- unlist(vapply(rownames(signifVariants),
                                function(x) res <- which(ASSBinfo[, 1] == x),
                                FUN.VALUE = integer(1)),
                                use.names=FALSE)
            
            ASSBinfo <- ASSBinfo[newindex, ]
        }
        
    } else {
        rown <- row.names(signifVariants)
        lengths2 <- lengths[rown, ]
    }
    
    deltapsi <- matrix(nrow=NROW(signifVariants), ncol=length(pairsCond))
    rownames(deltapsi) <- rownames(signifVariants)
    namesDeltaPsi <- c()
    ## initialize empty data frame to save the PSI values
    psi <- data.frame(ID=rownames(signifVariants))
    
    ## delta psi calculated for pairs of conditions, psi are calcuted for 
    ## each replicateXcondition
    for (pair in pairsCond) {
        ## for one pair
        index <- pair[[1]]
        if (is.list(index)) {
            condi <- namesCond[index[[1]]]
            replicates <- nr[index[[1]]]
        } else {
            condi <- namesCond[index]
            replicates <- nr[index]
        }
        psiPairCond <- matrix(nrow=NROW(signifVariants), ncol=sum(replicates))
        colsPsiPairCond <- c()
        indexMatrixPsiPairCond <- 1
        indexdeltapsi <- 1
        namesPsiPairCond <- c()
        nbLoop <- 0
        
        ## for a given condition in the pair
        for (nbRepli in seq_along(replicates)) {
            ## for each replicate (i) of the condition
            for (i in seq_len(replicates[nbRepli])) {
                nbLoop <- nbLoop + 1
                
                colsPsiPairCond <- 
                    c(colsPsiPairCond, paste(condi[nbRepli], 
                            "_repl", i, sep=""))
                
                namesUp <- 
                    c(paste("UP_", condi[nbRepli], 
                            "_repl", i, "_Norm", sep=""))
                
                namesLow <- 
                    c(paste("LP_", condi[nbRepli], 
                            "_repl", i, "_Norm", sep=""))
                
                ## the subsets are the counts we are going to use to compute 
                ## all psis
                subsetUp <- signifVariants[namesUp]
                subsetLow <- signifVariants[namesLow]
                ## counts correction
                if (!is.null(ASSBinfo) & !is.null(exonicReads) 
                    && !exonicReads) {
                    
                    ## ExonicReads=FALSE and counts=2
                    nameASSBinfo <- c(paste(condi[nbRepli], "_repl", i, sep=""))
                    subsetUp[which(subsetUp > 0), ] <- 
                        subsetUp[which(subsetUp > 0), ] / 
                        (2 - ASSBinfo[which(subsetUp > 0), nameASSBinfo] / 
                            subsetUp[which(subsetUp > 0), ])
                } else {
                    ## counts correction if there is no info about the junction
                    ## counts apparent size of upper path other apparent size 
                    ## of lower path. Default is: the user don't know the 
                    ## length for each events/eR=FALSE
                    correctFactor <- rep(2, length(lengths2$lower))
                    ## We adjust this factor for every other cases 
                    ## (eR=F and counts=1, eR=T and counts=0 or 2)
                    ## Note that if eR=F and counts=1, the correctFactor will 
                    ## be 2 or very close to 2
                    correctFactor[lengths2$lower!=0] <- 
                        lengths2$upper[lengths2$lower!=0] / 
                        lengths2$lower[lengths2$lower!=0]
                    subsetUp <- subsetUp / correctFactor
                }
                
                ## sumLowCond sums up the counts         
                sumLowCond[, nbRepli] <- 
                    sumLowCond[, nbRepli] + 
                    as.matrix(subsetUp) + 
                    as.matrix(subsetLow)
                
                ## psi is #incl/(#incl+#exclu) after all corrections for each 
                ## replicate
                psiPairCond[, indexMatrixPsiPairCond] <- 
                    as.matrix(subsetUp / (subsetUp + subsetLow))
                
                ## if counts are too low we will put NaN        
                indexNan <- 
                    intersect(which(subsetUp[, 1] < 10), 
                        which(subsetLow[, 1] < 10))
                
                psiPairCond[indexNan, nbLoop] <- NaN
                indexMatrixPsiPairCond <- indexMatrixPsiPairCond + 1
                namesPsiPairCond <- 
                    c(namesPsiPairCond, paste(as.character(condi[nbRepli]), 
                        "_repl", i, sep=""))
            }
        } 
        colnames(psiPairCond) <- namesPsiPairCond
        rownames(psiPairCond) <- rownames(signifVariants)
        rownames(sumLowCond) <- rownames(signifVariants)
        NaNSums <- rowSums((is.na(psiPairCond)) + 0)  ## 1 if NaN, 0 else
        
        ## if there are 2 NaN and 3 values for a bcc, nanSums is at 2 when
        ## there are more NaN than nb of column/2, we don't calculate the psi
        listNaN <- names(NaNSums[which(NaNSums > NCOL(psiPairCond) / 2)])
        psiPairCond[listNaN, ] <- NaN
        
        ## delta psi is the mean of the psis of the 2nd condition 
        ## (in terms of sorted condition) - the mean of the psis of the 1st 
        ## condition 
        deltaPsiCond <- 
            rowMeans(psiPairCond[, (replicates[1] + 1):sum(replicates)], 
                na.rm=TRUE) - 
            rowMeans(psiPairCond[, seq_len(replicates[1])], na.rm=TRUE)
        
        deltapsi[, indexdeltapsi] <- deltaPsiCond 
        indexdeltapsi <- indexdeltapsi + 1
        
        ## add the PSI in the data frame
        psi <- merge(psi, psiPairCond, by.x="ID", by.y="row.names")
    }
    
    ## when there are more than 2 conditions, we want to simplify the output:
    dPvector1 <- c(rep(0, NROW(signifVariants)))
    dPvector2 <- c(rep(0, NROW(signifVariants)))
    if (length(pairsCond) > 1) {
        ## if there are more than 2 conditions, we take the maximum of the 
        ## deltaPSI of all pairs
        for (l in seq_len(NROW(deltapsi))) {
            mindex <- which.max(abs(deltapsi[l, ]))
            if (length(mindex) != 0) {
                condA <- as.character(pairsCond[[mindex]][[1]][1])
                condB <- as.character(pairsCond[[mindex]][[1]][2])
                dP <- round(deltapsi[l, mindex], 4)
                dPvector1[l] <- dP
                dPvector2[l] <- paste(as.character(dP), 
                                    "(Cond", condB, ",", condA, ")", sep="")
            }
        }
        
    } else {
        dPvector1 <- round(deltapsi, 4)
        dPvector2 <- dPvector1
    }
    colnames(signifVariants) <- gsub("UP", "Variant1", colnames(signifVariants))
    colnames(signifVariants) <- gsub("LP", "Variant2", colnames(signifVariants))
    
    
    ###################################################
    ### code chunk 2: final table
    ###################################################
    signifVariants <- cbind(signifVariants, dPvector1)
    sortOrder <- 
        order(-abs(dPvector1), as.matrix(signifVariants[NCOL(signifVariants) - 1]))
    
    ## sorting by delta psi then by pvalue
    signifVariants.sorted <- signifVariants[sortOrder, ]
    
    names(dPvector2) <- rownames(signifVariants)
    dPvector2.sorted <- dPvector2[rownames(signifVariants.sorted)]
    signifVariants.sorted[NCOL(signifVariants.sorted)] <- dPvector2.sorted
    
    ## renaming last columns
    colnames(signifVariants.sorted)[length(colnames(signifVariants.sorted))] <-
        "Deltaf/DeltaPSI"
    colnames(signifVariants.sorted)[length(colnames(signifVariants.sorted)) - 
        1] <- "Adjusted_pvalue"
    class(signifVariants.sorted$Adjusted_pvalue) <- c("pval", 
        class(signifVariants.sorted$Adjusted_pvalue))
    
    ###################################################
    ### code chunk 3: flagging low counts
    ###################################################
    ## Condition to flag a low count for an event:
    lowcounts <- apply(sumLowCond, 1, function(x) 
        length(which(x < flagLowCountsConditions))) >= n - 1
    
    ## to fit the order with the sorted order
    
    lowcounts <- lowcounts[rownames(signifVariants.sorted)]
    
    ## final tab
    signifVariants.sorted <- cbind(signifVariants.sorted, lowcounts)
    colnames(signifVariants.sorted[NCOL(signifVariants.sorted)]) <- "Low_counts"
    
    ## return the final table and the table of PSI
    return(list(signifVariants.sorted=signifVariants.sorted,
        psiTable=psi))
}

.writeTableOutput <- function(finalTable, pvalMax=1, dPSImin=0, output) {
    fOut <- file(output, open="w")
    colNames <- colnames(finalTable)
    rowNames <- rownames(finalTable)
    nbRow <- length(rowNames)
    idDPSI <- length(colNames) - 1
    idPV <- length(colNames) - 2
    lHead <- c()
    j <- 1
    for (colName in colNames) {
        lHead <- append(x = lHead, values = paste(j, colName, sep="."))
        j <- j + 1
    }
    head <- paste(lHead, collapse="\t")
    head <- paste("#", head, sep="")
    writeLines(head,fOut)
    
    #fT <- finalTable
    #rownames(fT) <- NULL
    #colnames(fT) <- NULL
    i <- 1
    while (i <= nbRow) {
        apv <- finalTable[i, idPV]
        if (is.na(apv)) {
            apv <- 1
        }
        dpsi <- abs(finalTable[i, idDPSI])
        if (is.na(dpsi)) {
            dpsi <- 0
        }
        if (dpsi >= dPSImin && apv <= pvalMax) {
            writeLines(paste(rowNames[i], 
                paste(as.character(finalTable[i, ])[-1], collapse="\t"), 
                sep="\t"), fOut)
        }
        i <- i + 1
    }
    close(fOut)
} 

.writeMergeOutput <- function(resDiffExpr, k2rgFile, pvalMax=1, 
                        dPSImin=0, output) {
    
    K2RG_GENEID <- "Gene_Id"
    K2RG_GENENAME <- "Gene_name"
    K2RG_POS <- "Chromosome_and_genomic_position"
    K2RG_STRAND <- "Strand"
    K2RG_TYPE <- "Event_type"
    K2RG_VARLENGTH <- "Variabble_part_length"
    K2RG_FS <- "Frameshift_?"
    K2RG_CDS <- "CDS_?"
    K2RG_BIO <- "Gene_biotype"
    K2RG_KNOWNSS <- "number_of_known_splice_sites/number_of_SNPs"
    K2RG_BLOCUP <- "genomic_blocs_size_(upper_path)"
    K2RG_POSUP <-"genomic_position_of_each_splice_site_(upper_path)/of_each_SNP"
    K2RG_PARA <- "paralogs_?"
    K2RG_COMPLEX <- "Complex_event_?"
    K2RG_SNP <- "snp_in_variable_region"
    K2RG_EVENT <- "Event_name"
    K2RG_BLOCLOW <- "genomic_blocs_size_(lower_path)"
    K2RG_POSLOW <- "genomic_position_of_each_splice_site_(lower_path)"
    K2RG_PSI <- "Psi_for_each_replicate"
    K2RG_COVUP <- "Read_coverage(upper_path)"
    K2RG_COVLOW <- "Read_coverage(lower_path)"
    K2RG_CANON <- "Canonical_sites?"
    LK2RG <- c(K2RG_GENEID, K2RG_GENENAME, K2RG_POS, K2RG_STRAND, K2RG_TYPE, 
        K2RG_VARLENGTH, K2RG_FS, K2RG_CDS, K2RG_BIO, K2RG_KNOWNSS, 
        K2RG_BLOCUP, K2RG_POSUP, K2RG_PARA, K2RG_COMPLEX, K2RG_SNP, 
        K2RG_EVENT, K2RG_BLOCLOW, K2RG_POSLOW, K2RG_PSI, K2RG_COVUP, 
        K2RG_COVLOW, K2RG_CANON)
    
    K2RGKDE_APV <- "adjusted_pvalue"
    K2RGKDE_DPSI <- "dPSI"
    K2RGKDE_WARN <- "warnings"
    
    EVENTNAME <- 16
    COUNTSSTART <- 3
    COUNTSENDBEFOREEND <- 3
    
    finalTable <- resDiffExpr$finalTable
    finalTable <- finalTable[(is.na(finalTable$`Deltaf/DeltaPSI`) | 
        abs(finalTable$`Deltaf/DeltaPSI`)>=dPSImin) & 
        (is.na(finalTable$Adjusted_pvalue) | 
        finalTable$Adjusted_pvalue<=pvalMax), ]
    if(pvalMax<1){
        finalTable <- finalTable[!is.na(finalTable$Adjusted_pvalue),]
    }
    if(dPSImin>0){
        finalTable <- finalTable[!is.na(finalTable$`Deltaf/DeltaPSI`),]
    }
    psiTable <- resDiffExpr$`f/psiTable`
    finalAndPsiTable <- merge.data.frame(x = finalTable,
        y = psiTable, by.x = 1, by.y = 1)
    rownames(finalAndPsiTable) <- finalAndPsiTable$ID
    
    COUNTSEND <- ncol(finalTable)-COUNTSENDBEFOREEND
    PSISTART <- ncol(finalTable)+1
    PSIEND <- ncol(finalAndPsiTable)
    
    countsLabel <- 
        paste(colnames(finalAndPsiTable)[COUNTSSTART:COUNTSEND], collapse=",")
    countsName <- paste("CountsNorm(", countsLabel,")", sep="")
    psiLabel <- paste(colnames(finalAndPsiTable)[PSISTART:PSIEND], collapse=",")
    psiName <- paste("psiNorm(", psiLabel, ")", sep="")
    
    finalAndPsiTable$counts <- apply(X = 
        finalAndPsiTable[,COUNTSSTART:COUNTSEND], 
        MARGIN = 1, FUN = "paste", collapse=",")
    finalAndPsiTable$psis <- apply(X = finalAndPsiTable[,PSISTART:PSIEND],
        MARGIN = 1, FUN = "paste", collapse=",")
    finalAndPsiTable <- finalAndPsiTable[,colnames(finalAndPsiTable)%in%
        c("counts","psis","Adjusted_pvalue", "Deltaf/DeltaPSI","lowcounts")]
    finalAndPsiTable <- finalAndPsiTable[, c(4,5,1,2,3)]
    
    fK2RG <- file(k2rgFile, open="r")
    lines <- readLines(fK2RG)
    fOut <- file(output, open="w")
    
    line <- lines[1]
    if(startsWith(line, "#")) {
        nCol <- length(strsplit(x = line, split = "\t", fixed = TRUE)[[1]])
        countsHead <- paste(nCol+1, countsName, sep=".")
        psiHead <- paste(nCol+2, psiName, sep=".")
        pvHead <- paste(nCol+3, K2RGKDE_APV, sep=".")
        dPSIHead <- paste(nCol+4, K2RGKDE_DPSI, sep=".")
        warnHead <- paste(nCol+5, K2RGKDE_WARN, sep=".")
        toWrite <- paste(line, countsHead, psiHead, pvHead, 
            dPSIHead, warnHead, sep="\t")
        writeLines(toWrite, fOut)
        lines <- tail(lines, -1)
    }
    
    else {
        ## a header is created from the k2rg format with 22 columns
        lHead <- c()
        j <- 1
        for (k2rg_field in LK2RG) {
            lHead <- append(x = lHead, values = paste(j, k2rg_field, sep="."))
            j <- j + 1
        }
        
        lHead <- append(x = lHead, values = paste(j, countsName, sep="."))
        lHead <- append(x = lHead, values = paste(j + 1, psiName, sep="."))
        lHead <- append(x = lHead, values = paste(j + 2, K2RGKDE_APV, sep="."))
        lHead <- append(x = lHead, values = paste(j + 3, K2RGKDE_DPSI, sep="."))
        lHead <- append(x = lHead, values = paste(j + 4, K2RGKDE_WARN, sep="."))
        toWrite <- paste(lHead, collapse="\t")
        toWrite <- paste("#", toWrite, sep="")
        writeLines(toWrite,fOut)
    }
    
    # Transform lines in a data.frame
    t <- lapply(lapply(lapply(strsplit(x = lines, split = "\t"), "[[", 16),
                        strsplit, "\\|"), "[[", 1)
    t1 <- lapply(t,"[[",1)
    t2 <- lapply(t,"[[",2)
    lines <- data.frame(paste(t1,t2,sep="|"),lines)
    
    # Merge
    finalAndPsiTable <- merge.data.frame(x = lines,
                                    y = finalAndPsiTable,
                                    by.x = 1,
                                    by.y = 0)
    writeLines(paste(finalAndPsiTable[, 2], finalAndPsiTable[, 3],
            finalAndPsiTable[, 4], finalAndPsiTable[, 5], 
            finalAndPsiTable[, 6],finalAndPsiTable[, 7], sep="\t"), fOut)
    
    close(fK2RG)
    close(fOut)
}



.writePSITable <- function(resDiffExprVariant, adjPvalMax=1, dPSImin=0, output){
    
    selected_id <- 
        resDiffExprVariant$finalTable$ID[
            which(resDiffExprVariant$finalTable$Adjusted_pvalue <= adjPvalMax &
                abs(resDiffExprVariant$finalTable$`Deltaf/DeltaPSI`) >= 
                    dPSImin)]
    selected_PSI <- resDiffExprVariant$`f/psiTable`[selected_id,]
    write.table(resDiffExprVariant$`f/psiTable`, file=output, 
        quote=FALSE, sep="\t", row.names=FALSE)
}


.wantedEvents <- function(keep=c("All"), remove=NULL){
  if("-"%in%keep) {
    keep <- append(x = keep, values = c("-", "SNP","del"))
  }
  if("-"%in%remove) {
    remove <- append(x = remove, 
                     values = c("-", "SNP", "del"))
  }
  EVENTS <- c("deletion", "insertion", "IR", "ES", "altA", "altD", "altAD", 
              "-", "SNP", "indel")
  ES_EVENTS <- c("MULTI", "altA", "altD", "altAD", "MULTI_altA",
                 "MULTI_altD", "MULTI_altAD")
  wEvents <- c()
  if (keep == c("All") && is.null(remove)) {
    wEvents <- EVENTS
    for (i in seq_along(ES_EVENTS)) {
      wEvents <- append(x = wEvents, 
                        values = paste("ES_", ES_EVENTS[i], sep=""))
    }
    return(wEvents)
  }
  
  ES <- FALSE
  if (keep[1] == "All") {
    for (i in seq_along(EVENTS)) {
      if (!EVENTS[i] %in% remove) {
        wEvents <- append(x = wEvents, values = EVENTS[i])
      }
    }
    if ("ES" %in% remove) {
      ES <- TRUE
    }
    if (ES == FALSE) {
      for (i in seq_along(ES_EVENTS)) {
        wEvents <- append(x = wEvents, 
                          values = paste("ES_", ES_EVENTS[i], sep=""))
      }
    }
    return(wEvents)
  }
  for (i in seq_along(keep)) {
    if (ES == FALSE && keep[i] == "ES") {
      ES <- TRUE
    }
    wEvents <- append(x = wEvents, values = keep[i])
  }
  if (ES == FALSE) {
    return(wEvents)
  }
  if (is.null(remove)) {
    for (i in seq_along(ES_EVENTS)) {
      wEvents <- append(x = wEvents, 
                        values = paste("ES_", ES_EVENTS[i], sep=""))
    }
    return(wEvents)
  }
  for (i in seq_along(ES_EVENTS)) {
    if (!ES_EVENTS[i] %in% remove) {
      wEvents <- append(x = wEvents, 
                        values = paste("ES_", ES_EVENTS[i], sep=""))
    }
  }
  return(wEvents)
}
