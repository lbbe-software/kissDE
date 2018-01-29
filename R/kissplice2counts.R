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
    
    if (k2rg == FALSE) {
        fileNameK2RG <- NULL
        lines <- readLines(fpath)
        lines <- lines[startsWith(lines, ">")]
        lines <- sub("^>", "", lines)
        lines <- strsplit(x = lines, split = "|", fixed = TRUE)
        
        ## get all the informations for all lines
        infoLines <- lapply(lines, .getInfoLine, counts, pairedEnd, order, 
                            exonicReads)
        
        eventsnames <- unlist(lapply(infoLines, function(X) X$eventName))
        
        nbLines <- length(infoLines)
        eventsmat <- cbind(
            unlist(lapply(infoLines, function(X) X$variantLength)),
            matrix(unlist(lapply(infoLines, function(X) X$variantCounts)), nbLines, byrow = TRUE))
            
        infoFirstLine <- infoLines[[1]]
        psiInfo <- matrix(NA, nbLines, length(infoFirstLine$psiInfo))
        for (index in seq_len(nbLines))
            psiInfo[index, ] <- infoLines[[index]]$psiInfo
        
        events.df <- data.frame(events.names = eventsnames,
                               events.mat = eventsmat)
        
    } else {
        fileNameK2RG <- fileName
        EVENT <- 5
        EVENTNAME <- 16
        keepEvents <- .wantedEvents(keep, remove)
        
        lines <- readLines(fpath)
        lines <- lines[!startsWith(lines, "#")]
        lines <- strsplit(x = lines, split = "\t")
        ## keep only events of the selected type (keepEvents)
        keptLines <- lapply(lines, function(X) {if(X[EVENT] %in% keepEvents) X})
        keptLines <- keptLines[!vapply(keptLines, is.null, isTRUE(1))]
        lEvents <- unlist(lapply(keptLines, function(X) X[EVENTNAME]))
        
        nbLines <- length(keptLines)
        
        ## get all the informations for all lines
        infoLines <- lapply(keptLines, .getInfoLineK2rg, counts, pairedEnd,
                            order, exonicReads)
        
        lBcc <- unique(lEvents)
        iBcc <- length(lBcc)  ## number of unique bcc
        matBccApp <- matrix(0, nrow=iBcc) ## number of occurrences of each bcc
        rownames(matBccApp) <- lBcc
        
        ## Initialization of matrix
        eventsmat <- matrix(NA, iBcc * 2, length(infoLines[[1]]$variantCountsUp) + 1)
        eventsnames <- rep(NA, iBcc * 2)
        psiInfo <- matrix(NA, iBcc * 2, length(infoLines[[1]]$psiInfoUp))
        index <- 1
        indexNames <- 1
        for (index in seq_len(nbLines)){
            bcc <- keptLines[[index]][EVENTNAME]
            matBccApp[bcc, 1] <- matBccApp[bcc, 1] + 1
            if (matBccApp[bcc, 1] == 1) { # ça veut dire qu'on ne prend pas le deuxième ??
              resultLine <- infoLines[[index]]
              eventsmat[indexNames, 1] <- as.numeric(resultLine$variantLengthUp)
              eventsmat[indexNames, 2:NCOL(eventsmat)] <- resultLine$variantCountsUp
              eventsnames[indexNames] <- resultLine$eventName
              psiInfo[indexNames, ] <- resultLine$psiInfoUp
              
              eventsmat[indexNames + 1, 1] <- as.numeric(resultLine$variantLengthLow)
              eventsmat[indexNames + 1, 2:NCOL(eventsmat)] <- resultLine$variantCountsLow
              eventsnames[indexNames + 1] <- resultLine$eventName
              psiInfo[indexNames + 1, ] <- resultLine$psiInfoLow
              class(eventsmat) <- "numeric"
              indexNames <- indexNames + 2
            }
        }
        events.df <- data.frame(events.names = eventsnames,
                               events.mat = eventsmat)
    }
    
    ## update col names
    colnames(events.df) <- c("events.names", "events.length", paste("counts", 
        seq_len(length(colnames(events.df)) - 2), sep=""))
    
    close(fpath)
    psiInfo <- data.frame(events.names = eventsnames, as.data.frame(psiInfo))
    
    output <- list(countsEvents=events.df, psiInfo=psiInfo, 
        exonicReadsInfo=exonicReads, k2rgFile=fileNameK2RG)
    class(output) <- c("list", "countsData")
    return(output)
}

print.countsData <- function(x, ...) {
    print(x$countsEvents)
}
