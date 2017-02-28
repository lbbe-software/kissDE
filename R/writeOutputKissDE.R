writeOutputKissDE <- function(resDiffExprVariant, adjPvalMax=1, dPSImin=0, output, writePSI=FALSE) {
	if (adjPvalMax > 1 || adjPvalMax < 0) {
		print("ERROR: Invalid pvalMax (0 <= pvalMax <= 1).")
		return
	}
	if (dPSImin > 1 || dPSImin < 0) {
		print("ERROR: Invalid dPSImin (0 <= dPSImin <= 1). A dPSImin = 0.1 will catch all dPSI from -1 to -0.1 and all dPSI from 0.1 to 1.")
		return
	}
	k2rgFile <- resDiffExprVariant$k2rgFile
	if (writePSI){
		.writePSITable(resDiffExprVariant, adjPvalMax, dPSImin, output)
	} else{
		if (is.null(k2rgFile)) {
			.writeTableOutput(resDiffExprVariant$finalTable, adjPvalMax, dPSImin, output)
		}
		else {
			.writeMergeOutput(resDiffExprVariant, k2rgFile, adjPvalMax, dPSImin, output)
		}
	}
}