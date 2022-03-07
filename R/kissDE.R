kissDE <- function(fileName, conditions, output, counts = 2, pairedEnd = FALSE, order = NULL, exonicReads = TRUE, k2rg = FALSE, keep = c("All"), remove = NULL, pvalue = 1, 
                   filterLowCountsVariants = 10, flagLowCountsConditions = 10, technicalReplicates = FALSE, nbCore = 1, adjPvalMax = 1, dPSImin = 0, writePSI = FALSE, doQualityControl = TRUE, resultsInShiny=TRUE) {
  count <- kissplice2counts(fileName, counts, pairedEnd, order, exonicReads, k2rg, keep, remove)
  if(doQualityControl) {
    fileSplit <- strsplit(output,split = "/")[[1]]
    fold <- paste(fileSplit[1:(length(fileSplit)-1)],collapse = "/")
    qualityControl(countsData = count, conditions = conditions, storeFigs = fold)
  }
  res <- diffExpressedVariants(count, conditions, pvalue, filterLowCountsVariants, flagLowCountsConditions, technicalReplicates, nbCore)
  
  writeOutputKissDE(res, output, adjPvalMax, dPSImin)
  if(writePSI) {
    writeOutputKissDE(res, paste(output, "PSIs", sep = "."), adjPvalMax, dPSImin, writePSI)
  }
  
  if(resultsInShiny) {
    exploreResults(rdsFile = paste(output, "rds", sep = "."))
  }
}
