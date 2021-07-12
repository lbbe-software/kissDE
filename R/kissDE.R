kissDE <- function(fileName, conditions, output, counts=2, pairedEnd=FALSE, order=NULL, exonicReads=TRUE, k2rg=FALSE, keep=c("All"), remove=NULL, pvalue=1, filterLowCountsVariants=10, flagLowCountsConditions=10, technicalReplicates=FALSE, nbCore=1, adjPvalMax=1, dPSImin=0, writePSI=FALSE) {
  count=kissplice2counts(fileName, counts, pairedEnd, order, exonicReads, k2rg, keep, remove)
  res=diffExpressedVariants(count,conditions, pvalue, filterLowCountsVariants, flagLowCountsConditions, technicalReplicates, nbCore)
  writeOutputKissDE(resDiffExprVariant, output, adjPvalMax, dPSImin)
  if(writePSI) {
    writeOutputKissDE(resDiffExprVariant, paste(output,"PSIs",sep="."), adjPvalMax, dPSImin, writePSI)
  }
}
