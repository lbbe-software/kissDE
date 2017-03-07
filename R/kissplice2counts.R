kissplice2counts <- function(fileName, counts=0, pairedEnd=FALSE, order=NULL, 
    exonicReads=TRUE, k2rg=FALSE, keep=c("All"), remove=NULL) {
	## check options compatibility
	if (counts == 1 & exonicReads == TRUE) { 
		## when counts=1 set automatically exonicReads=TRUE
		exonicReads == FALSE
		warning("Changing 'exonicReads' value to FALSE for consistency with 
						counts=1.")
	}
	if (k2rg == FALSE & (keep != c("All") | !is.null(remove))) {
		## keep and remove should only be used when k2rg=TRUE
		keep <- c("All")
		remove <- NULL
		warning("Changing 'keep' and 'remove' options to default value for 
        consistency with k2rg=FALSE.")
	}
	
	toConvert <- file(fileName, open="r")
	nbLines <- countLines(fileName)
	if (k2rg == FALSE) {
		fileNameK2RG <- NULL
		index <- 1
		while (TRUE) {
			line <- readLines(toConvert, n=1)
			if (length(line) == 0) {
				break
			}
			if (substr(line, start=0, stop=1) != ">"){
				next
			}
			if (index == 1){
				isQuality <- grepl("Q", line[1])
			}
			## get all the informations for the line
			resultLine <- .getInfoLine(line, counts, pairedEnd, order, exonicReads, 
																 isQuality)
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
		
		keepEvents <- wantedEvents(keep, remove)
		
		index <- 1
		iEvents <- 0  ## number of unique and duplicated bcc = number of events
		lEvents <- list()
		while (TRUE) {
			line <- readLines(toConvert, n=1)
			if (length(line) == 0) {
				break
			}
			if(substr(line[1], 0, 1) == "#"){
				index <- index + 1
				next
			}
			bcc <- strsplit(line, split="\t")[[1]][EVENTNAME]
			if (strsplit(line, split="\t")[[1]][EVENT] %in% keepEvents){
				lEvents[iEvents + 1] <- bcc
				iEvents <- iEvents + 1
			}
			index <- index + 1
		}
		lBcc <- unique(lEvents)
		iBcc <- length(lBcc)  ## number of unique bcc
		matBccApp <- matrix(0, nrow=iBcc) ## number of occurrences of each bcc
		rownames(matBccApp) <- lBcc
		iDupBcc <- 1
		index <- 1
		indexNames <- 1
		seek(toConvert, 0) ## reinitialize the cursor at the beginning of the file
		while (TRUE) {
			line <- readLines(toConvert, n=1)
			if(length(line) == 0) {
				break
			}
			if(substr(line[1], 0, 1) == "#") {
				index <- index + 1
				indexNames <- 1
				next
			}
			lLine <- strsplit(line, split="\t")[[1]]
			if(lLine[EVENT] %in% keepEvents){
				bcc <- lLine[EVENTNAME]
				matBccApp[bcc, 1] <- matBccApp[bcc, 1] + 1
				if (matBccApp[bcc, 1] == 1) {
					resultLine <- .getInfoLineK2rg(line, counts, pairedEnd, order, 
																				 exonicReads)
					if (indexNames == 1){
						events.mat <- matrix(NA, iBcc * 2, 
																 length(resultLine$variantCountsUp) + 1)
						events.names <- rep(NA, iBcc * 2)
						psiInfo <- matrix(NA, iBcc * 2, length(resultLine$variantCountsUp))
					}
					resultLine$variantLengthUp <- 
						as.numeric(resultLine$variantLengthUp) + 
						as.numeric(resultLine$variantLengthLow)
					events.mat[indexNames, 1] <- as.numeric(resultLine$variantLengthUp)
					events.mat[indexNames, 2:NCOL(events.mat)] <- 
						resultLine$variantCountsUp
					events.names[indexNames] <- resultLine$eventName
					psiInfo[indexNames, ] <- resultLine$psiInfoUp
					events.mat[indexNames + 1, 1] <- as.numeric(resultLine$variantLengthLow)
					events.mat[indexNames + 1, 2:NCOL(events.mat)] <- 
						resultLine$variantCountsLow
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
	
	## update col names
	colnames(events.df) <- c("events.names", "events.length", paste("counts", 
													 			1:(length(colnames(events.df)) - 2), sep=""))
	
	close(toConvert)
	psiInfo <- data.frame(events.names, as.data.frame(psiInfo))
	
	output <- list(countsEvents=events.df, psiInfo=psiInfo, 
								 exonicReadsInfo=exonicReads, k2rgFile=fileNameK2RG)
	class(output) <- c("list", "countsData")
	return(output)
}

wantedEvents <- function(keep=c("All"), remove=NULL){
	EVENTS <- c("deletion", "insertion", "IR", "ES", "altA", "altD", "altAD", 
							"alt", "unclassified", "-", " ", "", "unclassifiedSNP")
	ES_EVENTS <- c("MULTI", "alt", "altA", "altD", "altAD")
	wEvents <- c()
	if (keep == c("All") && is.null(remove)) {
		wEvents <- EVENTS
		for (i in 1:length(ES_EVENTS)) {
			wEvents <- append(wEvents, paste("ES_", ES_EVENTS[i], sep=""))
		}
		return(wEvents)
	}
	
	if (!is.null(remove)) {
		for (i in 1:length(remove)) {
			if (!remove[i] %in% append(EVENTS, "MULTI")) {
				print(paste("In remove: couldn't find", remove[i]))
				stop("One of the element(s) of the remove vector is not part of: 
						 deletion, insertion, IR, ES, altA, altD, altAD, alt, unclassified, 
						 -, MULTI, unclassifiedSNP")
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
				wEvents <- append(wEvents, paste("ES_", ES_EVENTS[i], sep=""))
			}
		}
		return(wEvents)
	}
	for (i in 1:length(keep)) {
		if (!keep[i] %in% EVENTS) {
			print(paste("In keep: couldn't find", keep[i]))
			stop("One of the element(s) of the keep vector is not part of: 
					 deletion, insertion, IR, ES, altA, altD, altAD, alt, unclassified, 
					 -, unclassifiedSNP")
		}
		if (ES == FALSE && keep[i] == "ES") {
			ES <- TRUE
		}
		wEvents <- append(wEvents, keep[i])
	}
	if (ES == FALSE && !is.null(remove)) {
		stop("Keep and remove can not be set together, unless keep contain ES (in 
				 that case, remove will act on ES events)")
	}
	if (ES == FALSE) {
		return(wEvents)
	}
	if (is.null(remove)) {
		for (i in 1:length(ES_EVENTS)) {
			wEvents <- append(wEvents, paste("ES_", ES_EVENTS[i], sep=""))
		}
		return(wEvents)
	}
	for (i in 1:length(remove)){
		if (!remove[i] %in% ES_EVENTS) {
			print(paste("In remove: couldn't find",remove[i]))
			stop("One of the element(s) of the remove vector is not part of: 
					 altA, altD, altAD, alt, MULTI")
		}
	}
	for (i in 1:length(ES_EVENTS)) {
		if (!ES_EVENTS[i] %in% remove) {
			wEvents <- append(wEvents, paste("ES_", ES_EVENTS[i], sep=""))
		}
	}
	return(wEvents)
}


print.countsData <- function(x, ...) {
	print(x$countsEvents)
}
