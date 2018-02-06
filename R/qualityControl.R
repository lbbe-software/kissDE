qualityControl <- function(countsData, conditions, storeFigs=FALSE, 
                            returnPCAdata=FALSE) {
    
  # Check parameters
  if(!is.logical(returnPCAdata)) {
    stop("Input error: 'returnPCAdata' must be a boolean.")
  }
  
  ## Other checks in .readAndPrepareData 
  
    if (storeFigs == FALSE) {
        pathToFigs <- NA
    } else {
        if (isTRUE(storeFigs)) {
            pathToFigs <- paste0(tempdir(), "/kissDEFigures")
        } else {
            pathToFigs <- storeFigs
        }
    }
    
    if (!is.na(pathToFigs)) {
        if(!dir.exists(pathToFigs))
            dir.create(pathToFigs, recursive = TRUE)
        message("Figures are stored in ", pathToFigs)
        if (isTRUE(storeFigs)) {
            message("This directory is temporary.
  It will be removed when the R session is closed.")
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
    # remove the * in the condition vector
    if("*"%in%conditions){
        toRm <- which(conditions=="*")
        conditions <- conditions[-toRm]
    }
    conditionsNames <- sort(unique(conditions))
    n_conds <- length(conds)
    
    ###################################################
    ### select events with highest variance (on PSI)
    ###################################################
    countsData2 <- reshape(countsData[, c(1,(dimns):(dimns + n_conds))], 
        timevar="Path", idvar="ID", direction="wide")
    
    nr_cumsum <- c(1, 1 + cumsum(nr))
    nr_tot <- sum(nr)
    
    for (i in seq_len(n)) {
        nr_i <- seq_len(nr[i]) + nr_cumsum[i]
        denominator <- (countsData2[, nr_i] + countsData2[, nr_tot + nr_i])
        PSI <- countsData2[, nr_i] / denominator
        colnames(PSI) <- paste(conditionsNames[i], paste0("repl", seq_len(nr[i])), sep="_")
        countsData2 <- cbind(countsData2, PSI)
        
        nr_j <- seq_len(nr[i]) +  nr_cumsum[i]
        countsData2[, nr_tot * 2 + nr_j] <- do.call("cbind", lapply(nr_j,
               function(X) {
                   countsData2[intersect(which(countsData2[, X] < 10), which(countsData2[, nr_tot + X] < 10)), nr_tot * 2 + X] <- NA
                   countsData2[, nr_tot * 2 + X]
               }
        ))
    }
    
    countsData2$vars <- rowVars(
        as.matrix(countsData2[, ((nr_tot+1)*2):(nr_tot*2+1+n_conds)]),
        na.rm=TRUE)
    ## remove all NAs
    countsData2 <- countsData2[complete.cases(
        countsData2[, ((nr_tot+1)*2):(nr_tot*2+1+n_conds)]), ]
    ntop <- min(500, dim(countsData2)[1])
    selectntop <- order(countsData2$vars, decreasing=TRUE)[seq_len(ntop)]
    countsData2Selected <- countsData2[selectntop,]
    
    
    ###################################################
    ### code chunk number 3: heatmap
    ###################################################
    if (storeFigs == FALSE) {
        heatmap.2(
            as.matrix(as.dist(1 - cor(
                countsData2Selected[, 
                    ((nr_tot+1)*2):(nr_tot*2+1+n_conds)]))), 
            margins=c(10, 10), cexRow=1, cexCol=1, 
            density.info="none", trace="none")
        par(ask=TRUE)
    } else {
        filename <- paste(pathToFigs, "/heatmap.png", sep="")
        png(filename)
        heatmap.2(as.matrix(as.dist(1 - cor(
            countsData2Selected[, 
                ((nr_tot+1)*2):(nr_tot*2+1+n_conds)]))), 
            margins=c(10, 10), cexRow=1, cexCol=1, 
            density.info="none", trace="none")
        void <- dev.off()
    }
    
    ###################################################
    ### PCA plot
    ###################################################
    pca <- prcomp(t(
        countsData2Selected[, ((nr_tot+1)*2):(nr_tot*2+1+n_conds)]))
    pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
    pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
    pc1lab <- paste0("PC1 (", as.character(pc1var), "%)")
    pc2lab <- paste0("PC2 (", as.character(pc2var), "%)")
    # create the data.frame for ggplot
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=sort(conditions))
    if (storeFigs == FALSE) {
        p <- ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + 
            geom_point(size=3) + xlab(pc1lab) + ylab(pc2lab)
        print(p)
    } else {
        filename <- paste(pathToFigs, "/pca.png", sep="")
        p <- ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + 
            geom_point(size=3) + xlab(pc1lab) + ylab(pc2lab)
        ggsave(filename, plot = p, device = "png", width = 7, height = 7)
    }
    
    if (returnPCAdata == TRUE){
        return(d)
    }
}
