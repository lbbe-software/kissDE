.lineParse <- function(line, indexStart, isQuality) {
	beginningLineToWrite <- ""
	splitElements <- strsplit(line, "|", fixed=TRUE)[[1]]  ## splits the line
	if (indexStart == 6) {
		for (k in seq_len(4)) { ## indexStart - 2 = 4
			## writes the firsts elements of the line: bcc, cycle... 
			## but NOT branching_nodes
			beginningLineToWrite <- paste(beginningLineToWrite, 
																		splitElements[k], sep="|")
		}
	} else {
		for (k in seq_len(indexStart - 1)) {
			## writes the firsts elements of the line: bcc, cycle... 
			beginningLineToWrite <- paste(beginningLineToWrite, 
																		splitElements[k], sep="|")
		}
	}
	allcondi <- gregexpr(pattern="C[[:digit:]]+_[[:digit:]]+", text=line, 
											 perl=TRUE)
	if (allcondi[[1]][1] == -1) {  # no C in the header => ASSB counts
		allcondi <- gregexpr(pattern="[ASB]{1,4}[[:digit:]]+_[[:digit:]]+", 
												 text=line, perl=TRUE)
	}
	splittedCounts <- regmatches(x=line, m=allcondi)
	## gets the junctions id (ex 1 in AS1) and the count (ex 6 in AS1_6)
	countsperCond <- vapply(splittedCounts[[1]], .getJunctionCounts, 
													c("jct_id"="0", "count"="0"))
	
	return(list(beginning=beginningLineToWrite, 
							countsperCond=countsperCond))
}



.getJunctionCounts <- function(X){
	return(unlist(regmatches(x=X[[1]], 
													 m=gregexpr(pattern="[0-9]+", text=X[[1]]))))
}



.countsSetk2rg <- function(splittedCounts, counts=0, pairedEnd=FALSE, 
													 order=NULL, exonicReads=TRUE) {
	countsperCond <- vapply(splittedCounts[[1]], .getJunctionCounts, 
													c("jct_id"="0", "count"="0"))
	nbVec <- rep.int(0, dim(countsperCond)[2])
	countsVec <- rep.int(0, dim(countsperCond)[2])
	psiVec <- rep.int(0, dim(countsperCond)[2])
	for (i in seq_len(dim(countsperCond)[2])) {
		nbVec[i] <- as.numeric(countsperCond[1,i])
		countsVec[i] <- as.numeric(countsperCond[2,i])
		if (counts >= 1) {  ## specific issues linked with --counts option
			if (grepl("ASSB", colnames(countsperCond)[i]) == TRUE) {
				## so that counts on ASSB junction are not counted twice.
				psiVec[i] <- countsVec[i]
				countsVec[i] <- -countsVec[i]
			}
			if ((counts == 2) & (exonicReads == FALSE)) {
				if (grepl("^S[0-9]+", colnames(countsperCond)[i]) == TRUE) {
					## when exonic reads are not wanted we must discard 
					## reads counted in S_X
					countsVec[i] <- 0
				}
			}
		}
	}
	if (counts >= 1) {
		d <- data.frame(nbVec, countsVec)
		names(d) <- c("NB", "COUNTS")
		## sums the counts for each junction that belongs to the same event
		sums <- aggregate(d$COUNTS, by=list(d$NB), sum)
		## dpsi will store counts of ASSB counts 
		dpsi <- data.frame(nbVec, psiVec)
		names(dpsi) <- c("NB", "ASSB")
		assbPsi <- aggregate(dpsi$ASSB, by=list(dpsi$NB), sum)
		if (pairedEnd == TRUE) {
			if (is.null(order)) {
				order <- rep(seq_len((NROW(sums)) / 2), rep(2, ((NROW(sums)) / 2)))
			} else {
				if (!is.vector(order)) {
					stop("Order vector seems to be in a wrong format.")
				}
			}
			d2 <- data.frame(order, sums)
			names(d2)[3] <- "sums"
			## in case data is paired-end, there is one more sum to do, 
			## for each part of the pair
			sums2 <- aggregate(d2$sums, by=list(d2$order), sum)
			sums <- sums2
			dpsi2 <- data.frame(order, assbPsi)
			names(dpsi2) [3] <- "sums"
			assbPsi2 <- aggregate(dpsi2$sums, by=list(dpsi2$order), sum)
			assbPsi <- assbPsi2
		} 
		listASSB <- t(assbPsi)[2, ]
	} else {  # counts == 0 
		if (pairedEnd == TRUE) {
			if (is.null(order)) {
				## for length(s)=8, will create a vector c(1,1,2,2,3,3,4,4) 
				## (assuming data is ordered)
				order <- rep(seq_len(dim(countsperCond)[2] / 2), 
										 rep(2, dim(countsperCond)[2] / 2))
			} else {
				if (!is.vector(order)) {
					stop("Order vector seems to be in a wrong format.")
				}
			}
		} else {
			order <- c(seq_along(countsperCond[2,]))
		}
		d <- data.frame(order, countsVec)
		names(d) <- c("ORDER", "COUNTS")
		sums <- aggregate(d$COUNTS, by=list(d$ORDER), sum)
		listASSB <- NULL
	}
	listCounts <- t(sums)[2, ]
	
	return(list(vCounts=listCounts, 
							psiCounts=listASSB))
}



.countsSet <- function(line, indexStart, counts=0, pairedEnd=FALSE, 
											 order=NULL, exonicReads=TRUE, isQuality) {
	resultParsing <- .lineParse(line, indexStart, isQuality)
	beginningLineInfo <- resultParsing$beginning
	countsperCond <- resultParsing$countsperCond
	nbVec <- rep.int(0, dim(countsperCond)[2])
	countsVec <- rep.int(0, dim(countsperCond)[2])
	psiVec <- rep.int(0, dim(countsperCond)[2])
	for (i in seq_len(dim(countsperCond)[2])) {
		nbVec[i] <- as.numeric(countsperCond[1, i])
		countsVec[i] <- as.numeric(countsperCond[2, i])
		if (counts >= 1) {  ## specific issues linked with --counts option
			if (grepl("ASSB", colnames(countsperCond)[i]) == TRUE) {
				## so that counts on ASSB junction are not counted twice.
				psiVec[i] <- countsVec[i]
				countsVec[i] <- -countsVec[i]
			}
			if ((counts == 2) & (exonicReads == FALSE)) {
				if (grepl("^S[0-9]+", colnames(countsperCond)[i]) == TRUE) {
					## when exonic reads are not wanted we must discard 
					## reads counted in S_X
					countsVec[i] <- 0
				}
			}
		}
	}
	if (counts >= 1) {
		d <- data.frame(nbVec, countsVec)
		names(d) <- c("NB", "COUNTS")
		## sums the counts for each junction that belongs to the same event
		sums <- aggregate(d$COUNTS, by=list(d$NB), sum)
		## dpsi will store counts of ASSB counts 
		dpsi <- data.frame(nbVec, psiVec)
		names(dpsi) <- c("NB", "ASSB")
		assbPsi <- aggregate(dpsi$ASSB, by=list(dpsi$NB), sum)
		if (pairedEnd == TRUE) {
			if (is.null(order)) {
				order <- rep(seq_len((NROW(sums)) / 2), rep(2, ((NROW(sums)) / 2)))
			} else {
				if (!is.vector(order)) {
					stop("Order vector seems to be in a wrong format.")
				}
			}
			d2 <- data.frame(order, sums)
			names(d2)[3] <- "sums"
			## in case data is paired-end, there is one more sum to do, 
			## for each part of the pair
			sums2 <- aggregate(d2$sums, by=list(d2$order), sum)
			sums <- sums2
			dpsi2 <- data.frame(order, assbPsi)
			names(dpsi2)[3] <- "sums"
			assbPsi2 <- aggregate(dpsi2$sums, by=list(dpsi2$order), sum)
			assbPsi <- assbPsi2
		} 
		listASSB <- t(assbPsi)[2, ]
	} else {  ## counts == 0 
		if (pairedEnd == TRUE) {
			if (is.null(order)) {
				## for length(s)=8, will create a vector c(1,1,2,2,3,3,4,4) 
				## (assuming data is ordered)
				order <- rep(seq_len(dim(countsperCond)[2] / 2), 
										 rep(2, dim(countsperCond)[2] / 2))
			} else {
				if (!is.vector(order)) {
					stop("Order vector seems to be in a wrong format.")
				}
			}
		} else {
			order <- c(seq_len(dim(countsperCond)[2]))
		}
		d <- data.frame(order, countsVec)
		names(d) <- c("ORDER", "COUNTS")
		sums <- aggregate(d$COUNTS, by=list(d$ORDER), sum)
		listASSB <- NULL
	}
	listCounts <- t(sums)[2, ]
	
	return(list(firstPart=beginningLineInfo, 
							vCounts=listCounts, 
							psiCounts=listASSB))
}



.getInfoLine <- function(line, counts=0, pairedEnd=FALSE, order=NULL, 
												 exonicReads=TRUE, isQuality) {
	if (grepl("branching_nodes", line)) {
		indexStart <- 6 
	} else {
		indexStart <- 5
	}
	resultCountsSet <- .countsSet(line, indexStart, counts, pairedEnd, order, 
																exonicReads, isQuality)
	lineFirstPart <- resultCountsSet$firstPart
	lineFirstPartSplit <- strsplit(lineFirstPart, "|", fixed=TRUE)[[1]]
	eventName <- paste(lineFirstPartSplit[2], lineFirstPartSplit[3], sep="|")
	eventName <- substr(eventName, start=2, stop=nchar(eventName))
	variantLength <- strsplit(lineFirstPartSplit[5], "_")[[1]][4]
	
	return(list(eventName=eventName, 
							variantLength=variantLength, 
							variantCounts=resultCountsSet$vCounts, 
							psiInfo=resultCountsSet$psiCounts))
}



.getInfoLineK2rg <- function(line, counts=0, pairedEnd=FALSE, order=NULL, 
														 exonicReads=TRUE) {
	splitElements <- strsplit(line, "\t", fixed=TRUE)[[1]]
	firstPart <- splitElements[16]
	firstPartSplit <- strsplit(firstPart, "|", fixed=TRUE)[[1]]
	eventName <- paste(firstPartSplit[1], firstPartSplit[2], sep="|")
	countsUp <- strsplit(splitElements[20], ",")
	countsLow <- strsplit(splitElements[21], ",")
	resultCountsSetUp <- .countsSetk2rg(countsUp, counts, pairedEnd, 
																			order, exonicReads)
	resultCountsSetLow <- .countsSetk2rg(countsLow, counts, pairedEnd, 
																			 order, exonicReads)
	variantLengthUp <- as.integer(splitElements[6])
	variantLengthLow <- sum(as.integer(strsplit(splitElements[17],",")[[1]]))
	
	return (list(eventName=eventName, 
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
		for (j in seq_len(nr[n])) {
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
	
	## create a DESeqDataSet object
	suppressMessages(dds <- DESeqDataSetFromMatrix(
		countData=
			countsEvents[, !(names(countsEvents) %in% c("ID", "Length", "Path"))], 
		colData=data.frame(condition=conds),
		design=~ condition))
	
	ddsSF <- estimateSizeFactors(dds)
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
		ID=rep(as.factor(df["ID"]), endPosCol4Counts - startPosColumn4Counts + 1),
		cond=as.factor(unlist(lapply(
			strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts], "_"), 
			FUN=function(d){paste(d[2:(length(d) - 2)], collapse="_")}), 
			use.names=FALSE)), ## to manage the conditions names which contains "_"
		counts=as.numeric(df[startPosColumn4Counts:endPosCol4Counts]),
		path=as.factor(unlist(lapply(
			strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts], "|"), 
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



.fitNBglmModelsDSSPhi <- function(eventdata, phiDSS, phiDSScond, 
																	phiGlobal, nbAll){
	
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
  
  rslts <- c(1,       # [1]
             1,      # [2]
             nbAnov@anova.table$'P(> Chi2)'[2],        # [3]
             1,    # [4]	 
             c(9000,9000),       # [5:6] nbAIC0 
             c(9000,9000),      # [7:8] nbAICgb
             nbAIC,        # [9:10] nbAIC
             c(9000,9000),    # [11:12] nbAICcond
             
             c(0,0),      # [13:14]
             c(0,0),     # [15:16]
             nbCode,       # [17:18]
             c(0,0),   # [19:20]
             
             c(0,0),   # [21:22]
             c(0,0),  # [23:24]
             nbSingHes,    # [25:26]
             c(0,0))# [27:28]
  return(rslts)  
}



.modelFit <-function(countsData, n, nr, ASSBinfo, filterLowCountsVariants, techRep){
	##################################################
	## code chunk number 1: event-list
	##################################################
	# reduce data frame to the interesting columns
	nbAll <- sum(nr)
	dataPart <- countsData[, c(seq_len(2), 
														 which(grepl("_Norm", names(countsData))))]
	dataPart$Path <- gl(2, 1, NROW(countsData), labels=c("UP","LP"))
	dataPart2 <- cbind(dataPart[seq(1, NROW(dataPart), 2), ], 
										 dataPart[seq(2, NROW(dataPart), 2), 
										 				 grepl("Norm", names(dataPart))])
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
	# create list for the complete data set
	allEventtables <- apply(dataPart2, 1, 
													.eventtable, 
													startPosColumn4Counts=which(
														grepl("UP", names(dataPart2)))[1], 
													endPosCol4Counts=ncol(dataPart2))
	
	###################################################
	### code chunk number 2: DSS dispersion estimation
	###################################################
	## the counts matrix
	dataNormCountsEvent <- as.matrix(dataPart2[, 3:ncol(dataPart2)])
	colnames(dataNormCountsEvent) <- seq_len(ncol(dataNormCountsEvent))
	
	# We select the variant with the largest number of reads
	dataNormCountsEventSum <- matrix(nrow=nrow(dataNormCountsEvent), ncol=nbAll)
	rownames(dataNormCountsEventSum) <- rownames(dataNormCountsEvent)
	colnames(dataNormCountsEventSum) <- 1:nbAll
	for(i in c(1:nrow(dataNormCountsEvent))){
	  current <- dataNormCountsEvent[i,]
	  dataNormCountsEventSum[i,] <- current[1:(sum(nr))]+current[(sum(nr)+1):(2*sum(nr))]
	}
	
	designs <- rep(c(0,1),nr)
	dispData <- newSeqCountSet(dataNormCountsEventSum, as.data.frame(designs), 
	                           normalizationFactor = rep(1, nbAll))
	## fix the seed to avoid the stochastic outputs of the 
	## DSS:estDispersion function
	set.seed(40)
	dispData <- estDispersion(dispData)
	#hist(dispersion(dispData))
	names(exprs(dispData)) <- rownames(dataPart2)
	
	###################################################
	### code chunk number 3: variance - mean - Event level1
	###################################################
	## compute mean and variance per Event (instead of per allele)
	event.mean.variance.df <- as.data.frame(
		cbind(apply(dataPart2[, which(grepl("_Norm", names(dataPart2)))], 1, mean),
					apply(dataPart2[, which(grepl("_Norm", names(dataPart2)))], 1, var)))
	names(event.mean.variance.df) <- c("Mean", "Variance")
	rownames(event.mean.variance.df) <- as.character(dataPart2[, 1])
	## estimate the dispersion parameter D of the Quasi-Poisson distribution
	lm.D <- lm(event.mean.variance.df$Variance ~ event.mean.variance.df$Mean - 1)
	## estimate the overdispersion parameter theta of the NB distritution
	modelNB <- Variance ~ Mean + 1 / theta * Mean^2
	nls.modelNB <- nls(modelNB, data=event.mean.variance.df, 
										 start=list(theta=100))
	phi <- 1 / coef(nls.modelNB) ## to be used as fixed parameter later on
	
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
														grepl("UP", names(dataPart3)))[1], 
													endPosCol4Counts=ncol(dataPart3))
	
	###################################################
	### code chunk number 6: pALLGlobalPhi.glm.nb
	###################################################
	pALLGlobalPhiGlmNb <- data.frame(t(rep(NA, 28)))
	if(techRep) {
	  for (i in seq_along(allEventtables)) {
	    pALLGlobalPhiGlmNb[i, ] <- 
	      try(.fitNBglmModelsDSSPhi(allEventtables[[i]], 
	                                0, 
	                                0, 
	                                phi, 
	                                nbAll), 
	          silent=TRUE)
	  }
	}
	else {
	  for (i in seq_along(allEventtables)) {
	    pALLGlobalPhiGlmNb[i, ] <- 
	      try(.fitNBglmModelsDSSPhi(allEventtables[[i]], 
	                                dispersion(dispData)[i], 
	                                dispersion(dispData)[i], 
	                                phi, 
	                                nbAll), 
	          silent=TRUE)
	  }
	}
	
	
	###################################################
	### code chunk number 7: excl_errors
	###################################################
	sing.events <- which(grepl("Error", pALLGlobalPhiGlmNb[, 1]))
	if (length(sing.events) != 0) {
		pALLGlobalPhiGlmNb <- pALLGlobalPhiGlmNb[-sing.events, ]
	}
	colnames(pALLGlobalPhiGlmNb) <- c("(0)I vs A",
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
							phi=phi, 
							dispData=dispData))
}



.bestModelandSingular <- function(pALLGlobalPhi.glm.nb, sing.events, 
																	dataPart3, allEventtables, pvalue, phi, nr, 
																	dispData) {
	nbAll <- sum(nr)
	pALLGlobalPhi.glm.nb <- 
		pALLGlobalPhi.glm.nb[!is.na(pALLGlobalPhi.glm.nb[, 1]), ]
	
	if (NROW(pALLGlobalPhi.glm.nb) > 0) {
		matrixpALLGlobalPhi <- as.matrix(pALLGlobalPhi.glm.nb)
		storage.mode(matrixpALLGlobalPhi) <- "numeric"
		
		###################################################
		### code chunk number 1: best model
		###################################################
		bestmodel.table.n <- 
			apply(matrixpALLGlobalPhi[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, which.min)
		
		bestmodel.table <- bestmodel.table.n
		bestmodel.table[bestmodel.table == 1] <- "Poisson"
		bestmodel.table[bestmodel.table == 2] <- "Poisson"
		bestmodel.table[bestmodel.table == 3] <- "NB, global phi"
		bestmodel.table[bestmodel.table == 4] <- "NB, global phi"
		bestmodel.table[bestmodel.table == 5] <- "NB, DSS phi"
		bestmodel.table[bestmodel.table == 6] <- "NB, DSS phi"
		bestmodel.table[bestmodel.table == 7] <- "NB, cond DSS phi"
		bestmodel.table[bestmodel.table == 8] <- "NB, cond DSS phi"
		bestmodel.singhes <- c()
		for (i in seq_along(bestmodel.table.n)) {
			bestmodel.singhes[i] <- 
				c(matrixpALLGlobalPhi[i, c(21, 22, 23, 24, 25, 26, 27, 28)])[
					bestmodel.table.n[i]]
		}
		
		bestmodel.singhes <- unlist(bestmodel.singhes, use.names=FALSE)
		bestmodel <- table(bestmodel.table)
		
		###################################################
		### code chunk number 3: glmnet
		###################################################
		pALLGlobalPhi.glm.nb.glmnet <- as.data.frame(matrixpALLGlobalPhi)
		pALLGlobalPhi.glm.nb.glmnet$glmnet.pval <- 1
		pALLGlobalPhi.glm.nb.glmnet$glmnet.code <- 0
		
		## Variants for which the Poisson model is better
		singhes0 <- 
			which(apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, 
									which.min) == 1 | 
							apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, 
										which.min) == 2)
		
		singhes0_n <- names(singhes0)
		for (i in singhes0_n) {
			Xinter <- 
				model.matrix(~cond * path, 
										 data=allEventtables[[which(rownames(dataPart3) == i)]])
			
			outinter <- 
				glmnet(Xinter, 
							 allEventtables[[which(rownames(dataPart3) == i)]]$counts, 
							 family="poisson", lambda=1e-4, alpha=0)
			
			Xprinc <- 
				model.matrix(~path + cond, 
										 data=allEventtables[[which(rownames(dataPart3) == i)]])
			
			outprinc <- 
				glmnet(Xprinc, 
							 allEventtables[[which(rownames(dataPart3) == i)]]$counts, 
							 family="poisson", lambda=1e-4, alpha=0)
			
			Pv <- 1 - pchisq(deviance(outprinc) - deviance(outinter), df=1)
			pALLGlobalPhi.glm.nb.glmnet[i, "glmnet.pval"] <- Pv
			pALLGlobalPhi.glm.nb.glmnet[i, "glmnet.code"] <- outinter$jerr
		}
		matrixpALLGlobalPhi.glmnet <- as.matrix(pALLGlobalPhi.glm.nb.glmnet)
		storage.mode(matrixpALLGlobalPhi.glmnet) <- "numeric"
		
		###################################################
		###  code chunk number 4: Pseudo-counts and glmnet
		###################################################
		singhes <- 
			which(apply(matrixpALLGlobalPhi[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, 
									which.min) > 2 &
							rowSums(matrixpALLGlobalPhi[, c(21, 22, 23, 24,25, 26, 27, 28)])
						!= 0)
		
		singhes_n <- names(singhes)
		pALLGlobalPhi.glm.nb.pen <- as.data.frame(matrixpALLGlobalPhi)
		for (i in singhes_n) {
			pALLGlobalPhi.glm.nb.pen[i, ] <- 
				try(.fitNBglmModelsDSSPhi(.addOneCount(allEventtables[[i]]),
																	dispersion(dispData)[i],
																	dispersion(dispData)[i], phi, nbAll), 
						silent=TRUE)
		}
		
		pALLGlobalPhi.glm.nb <- as.data.frame(matrixpALLGlobalPhi)
		
		pALLGlobalPhi.glm.nb$final.pval.a.ia <- 1
		
		## Poisson model
		i <- 1
		li.singhes <- 
			which(apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, 
									which.min) == i | 
							apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, 
										which.min) == i+1)
		
		pALLGlobalPhi.glm.nb$final.pval.a.ia[li.singhes] <- 
			pALLGlobalPhi.glm.nb.glmnet$glmnet.pval[li.singhes]
		
		## negative binomial model, with global phi
		i <- 3
		li <- which((apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, 
											 which.min) == i | 
								 	apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1,
								 				which.min) == i+1 ) & 
									rowSums(pALLGlobalPhi.glm.nb[, c(21, 22, 23, 24, 25, 26, 27,
																									 28)]) == 0)
		pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li, 2]
		li <- which((apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, 
											 which.min) == i |
								 	apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, 
								 				which.min) == i+1 ) & 
									rowSums(pALLGlobalPhi.glm.nb[, c(21, 22, 23, 24, 25, 26, 27,
																									 28)]) != 0)
		pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb.pen[li, 2]
		
		## negative binomial model, phi estimated with DSS
		i <- 5
		li <- which((apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1, 
											 which.min) == i |
								 	apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1,
								 				which.min) == i+1 ) &
									rowSums(pALLGlobalPhi.glm.nb[, c(21, 22, 23, 24, 25, 26, 27,
																									 28)]) == 0)
		pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li, 3]
		li <- which((apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1,
											 which.min) == i |
								 	apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1,
								 				which.min) == i+1 ) &
									rowSums(pALLGlobalPhi.glm.nb[, c(21, 22, 23, 24, 25, 26, 27,
																									 28)]) != 0)
		pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb.pen[li, 3]
		
		## negative binomial model, phi estimated with DSS, conditionally to 
		## the expression mean    
		i <- 7
		li <- which((apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1,
											 which.min) == i |
								 	apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1,
								 				which.min) == i+1 ) &
									rowSums(pALLGlobalPhi.glm.nb[, c(21, 22, 23, 24, 25, 26, 27,
																									 28)]) == 0)
		pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li, 4]	 
		li <- which((apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1,
											 which.min) == i |
								 	apply(pALLGlobalPhi.glm.nb[, c(5, 6, 7, 8, 9, 10, 11, 12)], 1,
								 				which.min) == i+1 ) &
									rowSums(pALLGlobalPhi.glm.nb[, c(21, 22, 23, 24, 25, 26, 27,
																									 28)]) != 0)
		pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb.pen[li, 4]
		
		sing.events.final <- which(grepl("Error", pALLGlobalPhi.glm.nb[, 29]))
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



.sizeOfEffectCalc <- function(signifVariants, ASSBinfo, n, nr, 
															sortedconditions, flagLowCountsConditions, 
															lengths, exonicReads) {
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
			newindex <- unlist(sapply(rownames(signifVariants), 
																function(x) res <- which(ASSBinfo[, 1] == x)), 
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
					c(colsPsiPairCond, paste(condi[nbRepli], "_repl", i, sep=""))
				
				namesUp <- 
					c(paste("UP_", condi[nbRepli], "_repl", i, "_Norm", sep=""))
				
				namesLow <- 
					c(paste("LP_", condi[nbRepli], "_repl", i, "_Norm", sep=""))
				
				## the subsets are the counts we are going to use to compute all psis
				subsetUp <- signifVariants[namesUp]
				subsetLow <- signifVariants[namesLow]
				## counts correction
				if (!is.null(ASSBinfo)) {
					if (exonicReads){
						## apparent size of upper path other apparent size of lower path
						correctFactor <- lengths2$upper / lengths2$lower
						subsetUp <- subsetUp / correctFactor
					}
					else{
						nameASSBinfo <- c(paste(condi[nbRepli], "_repl", i, sep=""))
						subsetUp[which(subsetUp > 0), ] <- 
							subsetUp[which(subsetUp > 0), ] / 
							(2 - ASSBinfo[which(subsetUp > 0), nameASSBinfo] / 
							 	subsetUp[which(subsetUp > 0), ])
					}
					
				} else {
					## counts correction if there is no info about the junction counts
					
					## apparent size of upper path other apparent size of lower path
					correctFactor <- lengths2$upper / lengths2$lower
					subsetUp <- subsetUp / correctFactor
				}
				
				## sumLowCond sums up the counts         
				sumLowCond[, nbRepli] <- 
					sumLowCond[, nbRepli] + as.matrix(subsetUp) + as.matrix(subsetLow)
				
				## psi is #incl/(#incl+#exclu) after all corrections for each replicate
				psiPairCond[, indexMatrixPsiPairCond] <- 
					as.matrix(subsetUp / (subsetUp + subsetLow))
				
				## if counts are too low we will put NaN        
				indexNan <- 
					intersect(which(subsetUp[, 1] < 10), which(subsetLow[, 1] < 10))
				
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
		
		## if there are 2 NaN and 3 values for a bcc, nanSums is at 2
		## when there are more NaN than nb of column/2, we don't calculate the psi
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
		order(-abs(dPvector1), signifVariants[NCOL(signifVariants) - 1])
	
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
  K2RG_POSUP <- "genomic_position_of_each_splice_site_(upper_path)/of_each_SNP"
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
                                       y = psiTable,
                                       by.x = 1,
                                       by.y = 1)
  rownames(finalAndPsiTable) <- finalAndPsiTable$ID
  
  countsLabel <- 
    paste(colnames(finalAndPsiTable)[COUNTSSTART:COUNTSEND], collapse=",")
  countsName <- paste("CountsNorm(", countsLabel,")", sep="")
  psiLabel <- paste(colnames(finalAndPsiTable)[PSISTART:PSIEND], collapse=",")
  psiName <- paste("psiNorm(", psiLabel, ")", sep="")
  
  COUNTSEND <- ncol(finalTable)-COUNTSENDBEFOREEND
  PSISTART <- ncol(finalTable)+1
  PSIEND <- ncol(finalAndPsiTable)
  
  finalAndPsiTable$counts <- apply(X = finalAndPsiTable[,COUNTSSTART:COUNTSEND],
                                   MARGIN = 1,
                                   FUN = "paste",
                                   collapse=",")
  finalAndPsiTable$psis <- apply(X = finalAndPsiTable[,PSISTART:PSIEND],
                                 MARGIN = 1,
                                 FUN = "paste",
                                 collapse=",")
  finalAndPsiTable <- finalAndPsiTable[,colnames(finalAndPsiTable)%in%
                                         c("counts","psis","Adjusted_pvalue",
                                           "Deltaf/DeltaPSI","lowcounts")]
  finalAndPsiTable <- finalAndPsiTable[,c(4,5,1,2,3)]
  
  fK2RG <- file(k2rgFile, open="r")
  lines <- readLines(fK2RG)
  fOut <- file(output, open="w")
  
  line <- lines[1]
  if(substr(line[1], 0, 1) == "#"){
    nCol <- length(strsplit(line, split="\t")[[1]])
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
      lHead <- append(lHead, paste(j, k2rg_field, sep="."))
      j <- j + 1
    }
    
    lHead <- append(lHead, paste(j, countsName, sep="."))
    lHead <- append(lHead, paste(j + 1, psiName, sep="."))
    lHead <- append(lHead, paste(j + 2, K2RGKDE_APV, sep="."))
    lHead <- append(lHead, paste(j + 3, K2RGKDE_DPSI, sep="."))
    lHead <- append(lHead, paste(j + 4, K2RGKDE_WARN, sep="."))
    toWrite <- paste(lHead, collapse="\t")
    toWrite <- paste("#", toWrite, sep="")
    writeLines(toWrite,fOut)
  }
  
  # Transform lines in a data.frame
  t <- lapply(lapply(lapply(strsplit(lines,"\t"),"[[",16),strsplit,"\\|"),"[[",1)
  t1 <- lapply(t,"[[",1)
  t2 <- lapply(t,"[[",2)
  lines <- data.frame(paste(t1,t2,sep="|"),lines)
  
  # Merge
  finalAndPsiTable <- merge.data.frame(x = lines,
                                       y = finalAndPsiTable,
                                       by.x = 1,
                                       by.y = 0)
  writeLines(paste(finalAndPsiTable[,2],finalAndPsiTable[,3],finalAndPsiTable[,4],
                   finalAndPsiTable[,5],finalAndPsiTable[,6],finalAndPsiTable[,7],
                   sep="\t"),fOut)
  
  close(fK2RG)
  close(fOut)
}



.writePSITable <- function(resDiffExprVariant, adjPvalMax=1, 
													 dPSImin=0, output){
	
	selected_id <- 
		resDiffExprVariant$finalTable$ID[
			which(resDiffExprVariant$finalTable$Adjusted_pvalue <= adjPvalMax & 
							abs(resDiffExprVariant$finalTable$`Deltaf/DeltaPSI`) >= dPSImin)]
	selected_PSI <- resDiffExprVariant$`f/psiTable`[selected_id,]
	write.table(resDiffExprVariant$`f/psiTable`, file=output, 
							quote=FALSE, sep="\t", row.names=FALSE)
}


