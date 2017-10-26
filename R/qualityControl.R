qualityControl <- function(countsData, conditions, storeFigs=FALSE) {
    
    if (storeFigs == FALSE) {
        pathToFigs <- NA
    } else {
        if (isTRUE(storeFigs)) {
            pathToFigs <- paste0(getwd(), "/kissDEFigures")
        } else {
            pathToFigs <- storeFigs
        }
    }
    
    if (!is.na(pathToFigs)) {
        if(!dir.exists(pathToFigs))
            dir.create(pathToFigs, recursive = TRUE)
        message(paste("Figures are stored in", pathToFigs))
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
    conditionsNames <- sort(unique(conditions))
    
    ###################################################
    ### select events with highest variance (on PSI)
    ###################################################
    countsData2 <- reshape(countsData[, c(1,(dimns):(dimns + length(conds)))], 
        timevar="Path", idvar="ID", direction="wide")
    for (i in seq_len(n)) {
        for (j in seq_len(nr[i])) {
            countsData2$PSI <- countsData2[, (1+j+sum(nr[0:(i-1)]))] /
                (countsData2[, (1+j+sum(nr[0:(i-1)]))] + 
                    countsData2[, (1+sum(nr)+j+sum(nr[0:(i-1)]))])
            ## replace PSI by NA if count less or equal to 10 reads 
            ## for the 2 isoforms
            indexNA <- intersect(
                which(countsData2[, (1+j+sum(nr[0:(i-1)]))] < 10), 
                which(countsData2[, (1+sum(nr)+j+sum(nr[0:(i-1)]))] < 10))
            countsData2$PSI[indexNA] <- NA
            colnames(countsData2)[(sum(nr)*2+1)+j+sum(nr[0:(i-1)])] <- 
                paste(conditionsNames[i], paste0("repl", j), sep="_")
        }
    }
    countsData2$vars <- apply(
        as.matrix(countsData2[, ((sum(nr)+1)*2):(sum(nr)*2+1+length(conds))]),
        1, var, na.rm=TRUE)
    ntop <- min(500, dim(countsData2)[1])
    selectntop <- order(countsData2$vars, decreasing=TRUE)[seq_len(ntop)]
    countsData2Selected <- countsData2[selectntop,]
    ## remove all NAs
    countsData2Selected <- countsData2Selected[complete.cases(
        countsData2Selected[, ((sum(nr)+1)*2):(sum(nr)*2+1+length(conds))]), ]
    
    ###################################################
    ### code chunk number 3: heatmap
    ###################################################
    if (storeFigs == FALSE) {
        heatmap.2(
            as.matrix(as.dist(1 - cor(
                countsData2Selected[, 
                    ((sum(nr)+1)*2):(sum(nr)*2+1+length(conds))]))), 
            margins=c(10, 10), cexRow=1, cexCol=1, 
            density.info="none", trace="none")
        par(ask=TRUE)
    } else {
        filename <- paste(pathToFigs, "/heatmap.png", sep="")
        png(filename)
        heatmap.2(as.matrix(as.dist(1 - cor(
            countsData2Selected[, 
                ((sum(nr)+1)*2):(sum(nr)*2+1+length(conds))]))), 
            margins=c(10, 10), cexRow=1, cexCol=1, 
            density.info="none", trace="none")
        void <- dev.off()
    }
    
    ###################################################
    ### PCA plot
    ###################################################
    pca <- prcomp(t(
        countsData2Selected[, ((sum(nr)+1)*2):(sum(nr)*2+1+length(conds))]))
    fac <- factor(sort(conditions))
    colorpalette <- c("#192823", "#DD1E2F", "#EBB035", "#06A2CB", 
        "#218559", "#D0C6B1")
    colors <- colorpalette[seq_len(n)]
    pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
    pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
    pc1lab <- paste0("PC1 (", as.character(pc1var), "%)")
    pc2lab <- paste0("PC2 (", as.character(pc2var), "%)")
    if (storeFigs == FALSE) {
        par(oma=c(2, 1, 1, 1))
        plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, 
            xlab=pc1lab, ylab=pc2lab, main="PCA plot")
        par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
        plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
        legend("bottom", legend=levels(fac), xpd=TRUE, horiz=TRUE, 
            inset=c(0, 0), bty="n", pch=20, col=colors)
    } else {
        filename <- paste(pathToFigs, "/pca.png", sep="")
        png(filename)
        par(oma=c(2, 1, 1, 1))
        plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, 
            xlab=pc1lab, ylab=pc2lab, main="PCA plot")
        par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
        plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
        legend("bottom", legend=levels(fac), xpd=TRUE, horiz=TRUE, 
            inset=c(0, 0), bty="n", pch=20, col=colors)
        void <- dev.off()
    }
}
