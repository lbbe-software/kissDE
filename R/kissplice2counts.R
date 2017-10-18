kissplice2counts <- function(fileName, counts=0, pairedEnd=FALSE, order=NULL,
                             exonicReads=TRUE, k2rg=FALSE, keep=c("All"), 
                             remove=NULL) {
  ## check options 
  if(!file.exists(fileName)) {
    stop(paste("Input error : user's file ",fileName," does not exist. Is the path and/or file's name corect ?",sep=""))
  }
  if(!counts%in%[0,1,2]){
    stop("Input error : counts option is not equal to 0, 1 or 2.")
  }
  if(!logical(pairedEnd)){
    stop("Input error : pairedEnd option must be a boolean.")
  }
  if(!logical(exonicReads)){
    stop("Input error : exonicReads option must be a boolean.")
  }
  if(!logical(k2rg)){
    stop("Input error : k2rg option must be a boolean.")
  }
  if(!is.null(order) & !is.vector(order, mode="numeric")) {
    stop("Input error : order option must be a vector.")
  }
  if(k2rg & keep!=c("All")) {
    if(!is.vector(keep)) {
      stop("Input error : keep option must be a vector.")
    }
    keep <- unique(keep)
    for(element in keep) {
      if(!element%in%c("deletion", "insertion", "IR", "ES", "altA",
                       "altD", "altAD", "alt", "unclassified")) {
        stop(paste("Input error : element ",element," of option keep is not recognize. 
                   Each elements of the keep vector must be in the following list :
                   deletion, insertion, IR, ES, altA, altD, altAD, alt, unclassified.",
                   sep=""))
      }
    }
  }
  
	## check options compatibility
	if (counts == 1 & exonicReads == TRUE) { 
		## when counts=1 set automatically exonicReads=TRUE
		exonicReads <- FALSE
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
			resultLine <- .getInfoLine(line, counts, pairedEnd, order, 
			                           exonicReads, isQuality)
			eventName <- resultLine$eventName
			variantLength <- resultLine$variantLength
			variantCounts <- resultLine$variantCounts
			if (index == 1){
				events.mat <- matrix(NA, nbLines[1]/2, length(variantCounts)+1)
				events.names <- rep(NA, nbLines[1] / 2)
				psiInfo <- matrix(NA, nbLines[1]/2, length(resultLine$psiInfo))
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
		
		keepEvents <- .wantedEvents(keep, remove)
		
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
		## reinitialize the cursor at the beginning of the file
		seek(toConvert, 0) 
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
                    resultLine <- .getInfoLineK2rg(line, counts, pairedEnd, 
                                                   order, exonicReads)
					if (indexNames == 1){
						events.mat <- matrix(NA, iBcc * 2, 
                            length(resultLine$variantCountsUp) + 1)
						events.names <- rep(NA, iBcc * 2)
						psiInfo <- matrix(NA, iBcc * 2, 
						                  length(resultLine$psiInfoUp))
					}
					events.mat[indexNames, 1] <- 
					    as.numeric(resultLine$variantLengthUp)
					events.mat[indexNames, 2:NCOL(events.mat)] <- 
						resultLine$variantCountsUp
					events.names[indexNames] <- resultLine$eventName
					psiInfo[indexNames, ] <- resultLine$psiInfoUp
					events.mat[indexNames + 1, 1] <- 
						as.numeric(resultLine$variantLengthLow)
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
        seq_len(length(colnames(events.df)) - 2), sep=""))
	
	close(toConvert)
	psiInfo <- data.frame(events.names, as.data.frame(psiInfo))
	
	output <- list(countsEvents=events.df, psiInfo=psiInfo, 
	               exonicReadsInfo=exonicReads, k2rgFile=fileNameK2RG)
	class(output) <- c("list", "countsData")
	return(output)
}

print.countsData <- function(x, ...) {
	print(x$countsEvents)
}
