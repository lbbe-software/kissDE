writeOutputKissDE <- function(resDiffExprVariant, output, adjPvalMax=1, 
                        dPSImin=0, writePSI=FALSE) {
    
    ######## check function inputs
    
    if(length(output)>1 | !is.character(output)) {
        stop("Input error: 'output' must be a character.")
    }
    
    if ((length(adjPvalMax) > 1) | (!is.double(adjPvalMax)) | 
        (adjPvalMax < 0) | (adjPvalMax > 1)) {
        stop("Input error: 'adjPvalMax' must be a double between 0 and 1.")
    }
    
    if (length(dPSImin) > 1 | !is.double(dPSImin) | dPSImin < 0 | dPSImin > 1) {
        stop("Input error: 'dPSImin' must be a double between 0 and 1. 
    'dPSImin=0.1' will catch all dPSI from -1 to -0.1 and 
    all dPSI from 0.1 to 1.")
    }
    
    if(!is.logical(writePSI)) {
        stop("Input error: 'writePSI' must be a boolean.")
    }
    
    ########
    
    k2rgFile <- resDiffExprVariant$k2rgFile
    
    resDiffExprVariant$k2rgRes <- NA
    
    if (writePSI) {
        .writePSITable(resDiffExprVariant, adjPvalMax, dPSImin, output)
    } else {
        if (is.null(k2rgFile)) {
            .writeTableOutput(resDiffExprVariant$finalTable, adjPvalMax, 
                dPSImin, output)
        } else {
            .writeMergeOutput(resDiffExprVariant, k2rgFile, adjPvalMax, 
                dPSImin, output)
            resDiffExprVariant$k2rgRes <- output
        }
    }
    
    saveRDS(resDiffExprVariant,paste(output,"rds",sep="."))
}
