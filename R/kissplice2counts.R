kissplice2counts <- function(fileName, counts=2, pairedEnd=FALSE, order=NULL,
                            exonicReads=TRUE, k2rg=FALSE, keep=c("All"),
                            remove=NULL) {
    
    ######## check function inputs
    
    if(!file.exists(fileName)) {
        stop("Input error: user's file does not exist.")
    }
    if(!counts%in%c(0,1,2)){
        stop("Input error: 'counts' option is not equal to 0, 1 or 2.")
    }
    if(!is.logical(pairedEnd)){
        stop("Input error: 'pairedEnd' option must be a boolean.")
    }
    if(!is.logical(exonicReads)){
        stop("Input error: 'exonicReads' option must be a boolean.")
    }
    if(!is.logical(k2rg)){
        stop("Input error: 'k2rg' option must be a boolean.")
    }
    if(!is.null(order) & !is.vector(order, mode="numeric")) {
        stop("Input error: 'order' option must be a numeric vector.")
    }
    if(!pairedEnd & !is.null(order)) {
        stop("Input error: 'order' option can be set only if 
            pairedEnd is TRUE.")
    }
    
    if(k2rg & keep!=c("All")) {
        if(!is.vector(keep)) {
            stop("Input error: 'keep' option must be a vector.")
        }

        keep <- unique(keep)
        values <- c("deletion", "insertion", "IR", "ES", "altA", "altD", 
            "altAD", "-", "indel")
        
        # vector of unauthorized 'keep' values
        notinValues <- keep[!keep%in%values]
        if(length(notinValues) > 0) {
            notinValues <- paste(notinValues, collapse = ", ")
            stop("Input error: 'keep' must be a vector of values in 
    c('deletion', 'insertion', 'IR', 'ES', 'altA', 'altD', 'altAD', 
    '-', 'indel').")
        }
    }
    
    if(k2rg & !is.null(remove)) {
        if(!is.vector(remove)) {
            stop("Input error: 'remove' must be a vector.")
        }

        if(!"ES"%in%keep) {
            remove <- unique(remove)
            values <- c("deletion", "insertion", "IR", "ES", "altA", "altD", 
                        "altAD", "-", "indel", "MULTI")
            
            # vector of unauthorized 'remove' values
            notinValues <- remove[!remove%in%values]
            if(length(notinValues) > 0) {
                notinValues <- paste(notinValues, collapse = ", ")
                stop("Input error: 'remove' must be a vector of values in 
        c('deletion', 'insertion', 'IR', 'ES', 'altA', 'altD', 'altAD', 
        '-', 'indel', 'MULTI').")
            }
        }
        
        
        if(keep!=c("All") & "ES"%in%keep) {
            values <- c("altA", "altD", "altAD", "MULTI",
                       "MULTI_altA", "MULTI_altD", "MULTI_altAD")
            
            # vector of unauthorized 'remove' values
            notinValues <- remove[!remove%in%values]
            if(length(notinValues) > 0) {
                notinValues <- paste(notinValues, collapse = ", ")
                stop("Input error: if 'keep' contains 'ES', 'remove' must be a 
    vector of values in: c('altA', 'altD', 'altAD', 'MULTI', 'MULTI_altA',
    MULTI_altD', 'MULTI_altAD').")
            }
        } else if(keep!=c("All")) {
            stop("Input error: 'keep' and 'remove' can not be used together,
    unless 'keep' contains 'ES' (see the vignette for more informations).")
        }
    }
    
    ######## check options compatibility
  
    splitFileName <- strsplit(fileName, "_type_")[[1]]
    l <- length(splitFileName)
    if(l > 1) {
      fileType <- substr(splitFileName[l], 1, 1)
      if(fileType == "0" & counts != 0) {
        counts <- 0
        warning("Changing 'counts' value to 0 for consistency with inputed file type.")
      }
    }
    
    if (counts == 1 & exonicReads == TRUE) { 
        ## when counts=1 set automatically exonicReads=FALSE
        exonicReads <- FALSE
        warning("Changing 'exonicReads' value to FALSE for consistency with
    counts=1.")
    }
    if (counts == 0 & exonicReads == FALSE) { 
        ## when counts=0 set automatically exonicReads=TRUE
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
        
        
        if (any(grepl("branching_nodes", lines[1]))) {
            indexStart <- 6
        } else {
            indexStart <- 5
        }
        
        ## get all the informations for all lines
        infoLines <- lapply(lines, .getInfoLine, counts, pairedEnd, order, 
                            exonicReads, indexStart)
        
        eventsnames <- unlist(lapply(infoLines, function(X) X$eventName))
        
        nbLines <- length(infoLines)
        eventsmat <- cbind(
            unlist(lapply(infoLines, function(X) X$variantLength)),
            matrix(unlist(lapply(infoLines, function(X) X$variantCounts)), 
                nbLines, byrow = TRUE))
            
        infoFirstLine <- infoLines[[1]]
        psiInfo <- matrix(NA, nbLines, length(infoFirstLine$psiInfo))
        for (index in seq_len(nbLines))
            psiInfo[index, ] <- infoLines[[index]]$psiInfo
        
    } else {
        fileNameK2RG <- fileName
        POS <- 3
        EVENT <- 5
        PARA <- 13
        EVENTNAME <- 16
        keepEvents <- .wantedEvents(keep, remove)
        
        lines <- readLines(fpath)
        lines <- lines[!startsWith(lines, "#")]
        lines <- strsplit(x = lines, split = "\t", fixed = TRUE)
        
        ## keep only events of the selected type (keepEvents)
        keptLines <- lapply(lines, function(X) {if(!X[PARA] %in% c("one_path_mapped","not_mapped","missed_type2") & X[EVENT] %in% keepEvents) X})
        keptLines <- keptLines[!vapply(keptLines, is.null, isTRUE(1))]
        keptLines <- keptLines[!duplicated(unlist(lapply(keptLines, 
            function(X) X[EVENTNAME])))]
        lEvents <- unlist(lapply(keptLines, function(X) X[EVENTNAME]))
        
        nbLines <- length(keptLines)
        
        ## get all the informations for all lines
        infoLines <- lapply(keptLines, .getInfoLineK2rg, counts, pairedEnd,
                            order, exonicReads)
        
        ## events.names
        eventsnames <- rep(NA, nbLines * 2)
        eventsnames[seq(1, nbLines * 2, by = 2)] <- unlist(lapply(infoLines, 
            function(X) X$eventName))
        eventsnames[seq(2, nbLines * 2, by = 2)] <- unlist(lapply(infoLines, 
            function(X) X$eventName))
        
        ## events.mat
        infoFirstLine <- infoLines[[1]]
        eventsmat <- matrix(NA, nbLines * 2, 
            length(infoFirstLine$variantCountsUp) + 1)
        eventsmat[seq(1, nbLines * 2, by = 2), ] <- cbind(
            unlist(lapply(infoLines, function(X) X$variantLengthUp)),
            matrix(unlist(lapply(infoLines, function(X) X$variantCountsUp)), 
                nbLines, byrow = TRUE))
        eventsmat[seq(2, nbLines * 2, by = 2), ] <- cbind(
            unlist(lapply(infoLines, function(X) X$variantLengthLow)),
            matrix(unlist(lapply(infoLines, function(X) X$variantCountsLow)), 
                nbLines, byrow = TRUE))
        
        ## psiInfo
        psiInfo <- matrix(NA, nbLines * 2, length(infoFirstLine$psiInfoUp))
        # For loop that can not be optimized at the moment because of NULL 
        # value yielded by the '.countsSet' function
        for (index in seq_len(nbLines)) {
            psiInfo[(index - 1)*2 + 1, ] <- infoLines[[index]]$psiInfoUp
            psiInfo[(index - 1)*2 + 2, ] <- infoLines[[index]]$psiInfoLow
        }
        
        
    }
    
    eventsdf <- data.frame(events.names = eventsnames, events.mat = eventsmat)
    colnames(eventsdf) <- c("events.names", 
                            "events.length", 
                            paste("counts", 
                                seq_len(length(colnames(eventsdf)) - 2), 
                                sep=""))
    
    psiInfo <- data.frame(events.names = eventsnames, as.data.frame(psiInfo))
    
    output <- list(countsEvents = eventsdf, 
                    psiInfo = psiInfo, 
                    exonicReadsInfo = exonicReads, 
                    k2rgFile = fileNameK2RG)
    class(output) <- c("list", "countsData")
    
    close(fpath)
    return(output)
}

print.countsData <- function(x, ...) {
    print(x$countsEvents)
}
