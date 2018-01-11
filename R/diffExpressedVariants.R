diffExpressedVariants <- function(countsData, conditions, pvalue=1, 
                                    filterLowCountsVariants=10, 
                                    flagLowCountsConditions=10,
                                    technicalReplicates=FALSE) {
    
    ######## check function inputs
    
    if(length(pvalue)>1 | !is.double(pvalue) | pvalue<0 | pvalue>1) {
        stop("Input error: 'pvalue' must be a double between 0 and 1.")
    }
    
    if((length(filterLowCountsVariants)>1) | 
        (filterLowCountsVariants-round(filterLowCountsVariants+0.5)!=0) | 
        (filterLowCountsVariants<0)) {
            stop("Input error: 'filterLowCountsVariants' must be a 
    positive integer.")
    }
    
    if((length(flagLowCountsConditions)>1) | 
        (flagLowCountsConditions-round(flagLowCountsConditions+0.5)!=0) | 
        (flagLowCountsConditions<0)) {
            stop("Input error: 'flagLowCountsConditions' must be a 
    positive integer.")
    }
    
    if(!is.logical(technicalReplicates)) {
        stop("Input error: 'technicalReplicates' must be a boolean.")
    }
    
    ########
    
    message("Pre-processing the data...")
    chunk0 <- tryCatch({.readAndPrepareData(countsData, conditions)
        #### chunk 0 var ####
        ## chunk0$countsData
        ## chunk0$conditions
        ## chunk0$dim
        ## chunk0$n
        ## chunk0$nr
        ## chunk0$sortedconditions
        ## chunk0$ASSBinfo
    }, error=function(err) {
        return(NA)
        stop(err)
    })
    
    if (!is.na(chunk0[1])) {  # no error in chunk 0
        ## in case counts option in kissplice2counts is at 1 or 2, 
        ## we have info about junction counts (ASSB), that will
        ## be useful to correct the computation of delta psi in the end.
        ## They are stored here.
        ASSBinfo <- chunk0$ASSBinfo  
        if (!is.null(ASSBinfo)) {
            li <- seq_len(NROW(ASSBinfo))
            li <- li[which(li%%2 != 0)] ## select odd indexes
            ASSBinfo <- ASSBinfo[li, ]
        }
        message("Trying to fit models on data...")
        message("This can be a time-consuming step, so do not hesitate to have 
    a look at the very well-written vignette !")
        chunk1 <- tryCatch({.modelFit(chunk0$countsData, chunk0$n, chunk0$nr, 
                                        ASSBinfo, filterLowCountsVariants,
                                        technicalReplicates)
            #### chunk 1 var ####
            ## chunk1$pALLGlobalPhi.glm.nb 
            ## chunk1$sing.events
            ## chunk1$dataPart3
            ## chunk1$ASSBinfo
            ## chunk1$allEventtables
            ## chunk1$length
            ## chunk1$phi
            ## chunk1$dispData
        }, error=function(err) {
            return(NA)
            stop(err, " An error occured, unable to fit models on data.")
        }) 
    } else {  # error in chunk 0
        chunk1 <- NA
    }
    
    if (!is.na(chunk1[1])) {  # no error in chunk 1 nor in chunk 0
        message("Computing pvalues...")
        chunk2 <- tryCatch({.bestModelandSingular(chunk1$pALLGlobalPhi.glm.nb, 
                                                    chunk1$sing.events, 
                                                    chunk1$dataPart3, 
                                                    chunk1$allEventtables, 
                                                    pvalue, 
                                                    chunk1$phi,
                                                    chunk0$nr, 
                                                    chunk1$dispData)
            #### chunk 2 var ####  
            ## chunk2$noCorrectPVal
            ## chunk2$correctedPVal
            ## chunk2$signifVariants
        }, error=function(err) {
            message(err, " Returning only 'resultFitNBglmModel' 
    and 'sing.events'")
            return(list(resultFitNBglmModel=chunk1$pALLGlobalPhi.glm.nb, 
                        sing.events=chunk1$sing.events))
        })
    } else {
        chunk2 <- NA
    }
    
    if (!is.na(chunk2[1])) {  # no error during chunk1
        if (length(chunk2) > 2) {  # no error during chunk2
            message("Computing size of the effect and last cutoffs...")
            class(chunk2$correctedPVal) <- c("pval", 
                                                class(chunk2$correctedPVal))
            class(chunk2$noCorrectPVal) <- c("pval", 
                                                class(chunk2$noCorrectPVal))
            chunk3 <- tryCatch({
                sizeOfEffect <- .sizeOfEffectCalc(chunk2$signifVariants, 
                                                chunk1$ASSBinfo, chunk0$n, 
                                                chunk0$nr, 
                                                chunk0$sortedconditions, 
                                                flagLowCountsConditions, 
                                                chunk1$lengths, 
                                                countsData$exonicReadsInfo)
                return(list(finalTable=sizeOfEffect$signifVariants.sorted, 
                            correctedPVal=chunk2$correctedPVal, 
                            uncorrectedPVal=chunk2$noCorrectPVal, 
                            resultFitNBglmModel=chunk1$pALLGlobalPhi.glm.nb,
                            `f/psiTable`=sizeOfEffect$psiTable,
                            k2rgFile=countsData$k2rgFile))
            }, error=function(err) {
                message(err, " Returning only 'resultFitNBglmModel' and 
                    pvalues tab")
                return(list(correctedPVal=chunk2$correctedPVal,
                            uncorrectedPVal=chunk2$noCorrectPVal,
                            resultFitNBglmModel=chunk1$pALLGlobalPhi.glm.nb))
            })
        } else { # error in chunk 2 does not allow to compute chunk 3
            return(chunk2)
        }
    } else {
        return(NA)
    }
}

print.pval <- function(x, ...) {
    x <- format.pval(x)
    print(x, quote = FALSE)
}
