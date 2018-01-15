kissplice2counts <- function(fileName, counts=0, pairedEnd=FALSE, order=NULL,
                            exonicReads=TRUE, k2rg=FALSE, keep=c("All"),
                            remove=NULL) {
    
    ######## check function inputs
    
    if(!file.exists(fileName)) {
        stop("Input error: user's file does not exist.")
    }
    if(!counts%in%c(0,1,2)){
        stop("Input error: counts option is not equal to 0, 1 or 2.")
    }
    if(!is.logical(pairedEnd)){
        stop("Input error: pairedEnd option must be a boolean.")
    }
    if(!is.logical(exonicReads)){
        stop("Input error: exonicReads option must be a boolean.")
    }
    if(!is.logical(k2rg)){
        stop("Input error: k2rg option must be a boolean.")
    }
    if(!is.null(order) & !is.vector(order, mode="numeric")) {
        stop("Input error: order option must be a numeric vector.")
    }
    if(!pairedEnd & !is.null(order)) {
        stop("Input error: order option can be set only if pairedEnd is TRUE.")
    }
    
    if(k2rg & keep!=c("All")) {
        if(!is.vector(keep)) {
            stop("Input error: keep option must be a vector.")
        }

        keep <- unique(keep)
        values <- c("deletion", "insertion", "IR", "ES", "altA", "altD", 
            "altAD", "alt", "unclassified")
        
        # vector of unauthorized 'keep' values
        notinValues <- keep[!keep%in%values]
        if(length(notinValues) > 0) {
            notinValues <- paste(notinValues, collapse = ", ")
            stop("Input error: 'keep' must be a vector of values in 
    c('deletion', 'insertion', 'IR', 'ES', 'altA', 'altD', 'altAD', 
    'alt', 'unclassified'.")
        }
    }
    
    if(k2rg & !is.null(remove)) {
        if(!is.vector(remove)) {
            stop("Input error: 'remove' must be a vector.")
        }

        remove <- unique(remove)
        values <- c("deletion", "insertion", "IR", "ES", "altA", "altD", 
                    "altAD", "alt", "unclassified", "MULTI")
        
        # vector of unauthorized 'remove' values
        notinValues <- remove[!remove%in%values]
        if(length(notinValues) > 0) {
            notinValues <- paste(notinValues, collapse = ", ")
            stop("Input error: 'remove' must be a vector of values in 
    c('deletion', 'insertion', 'IR', 'ES', 'altA', 'altD', 'altAD', 
    'alt', 'unclassified', 'MULTI'.")
        }
        
        if(keep!=c("All") & "ES"%in%keep) {
            values <- c("altA", "altD", "altAD", "alt", "MULTI")
            
            # vector of unauthorized 'remove' values
            notinValues <- remove[!remove%in%values]
            if(length(notinValues) > 0) {
                notinValues <- paste(notinValues, collapse = ", ")
                stop("Input error: if 'keep' contains 'ES', 'remove' must be a 
    vector of values in: c('altA', 'altD', 'altAD', 'alt', 'MULTI'.")
            }
        } else if(keep!=c("All")) {
            stop("Input error: 'keep' and 'remove' can not be used together,
    unless 'keep' contains 'ES' (see the vignette for more informations).")
        }
    }
    
    ######## check options compatibility
    
    if (counts == 1 & exonicReads == TRUE) { 
        ## when counts=1 set automatically exonicReads=TRUE
        exonicReads <- FALSE
        warning("Changing 'exonicReads' value to FALSE for consistency with
    counts=1.")
    }
    if (counts == 0 & exonicReads == FALSE) { 
        ## when counts=1 set automatically exonicReads=TRUE
        exonicReads <- TRUE
        warning("Changing 'exonicReads' value to TRUE for consistency with
    counts=0.")
    }
    if (k2rg == FALSE & (keep != c("All") | !is.null(remove))) {
        ## keep and remove should only be used when k2rg=TRUE
        keep <- c("All")
        remove <- NULL
        warning("Changing 'keep' and 'remove' options to default value for 
    consistency with k2rg=FALSE.")
    }
    
    ########
    
    fpath <- file(fileName, open="r")
    nbLines <- countLines(fileName)
    if (k2rg == FALSE) {
        fileNameK2RG <- NULL
        lines <- readLines(fpath)
        lines <- lines[startsWith(lines, ">")]
        lines <- sub("^>", "", lines)
        lines <- strsplit(lines, "|", fixed=TRUE)
        ## get all the informations for all lines
        infoLines <- lapply(lines, .getInfoLine, counts, pairedEnd, order, 
                            exonicReads)
        lengthVariantCounts <- length(unlist(infoLines[[1]][3]))
        lengthPsiInfo <- length(unlist(infoLines[[1]][4]))
        events.mat <- matrix(NA, nbLines/2, lengthVariantCounts+1)
        events.names <- rep(NA, nbLines/ 2)
        psiInfo <- matrix(NA, nbLines/2, lengthPsiInfo)
        index <- 1
        while (index <= nbLines/2) {
            event <- infoLines[[index]]
            events.names[index] <- event[[1]]
            events.mat[index, 1] <- as.numeric(event[[2]])
            events.mat[index, 2:NCOL(events.mat)] <- event[[3]]
            psiInfo[index, ] <- event[[4]]
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
            line <- readLines(fpath, n=1)
            if (length(line) == 0) {
                break
            }
            if(startsWith(line, "#")) {
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
        seek(fpath, 0) 
        while (TRUE) {
            line <- readLines(fpath, n=1)
            if(length(line) == 0) {
                break
            }
            if(startsWith(line, "#")) {
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
    
    close(fpath)
    psiInfo <- data.frame(events.names, as.data.frame(psiInfo))
    
    output <- list(countsEvents=events.df, psiInfo=psiInfo, 
        exonicReadsInfo=exonicReads, k2rgFile=fileNameK2RG)
    class(output) <- c("list", "countsData")
    return(output)
}

print.countsData <- function(x, ...) {
    print(x$countsEvents)
}
