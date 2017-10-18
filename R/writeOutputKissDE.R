writeOutputKissDE <- function(resDiffExprVariant, output, adjPvalMax=1, 
                              dPSImin=0, writePSI=FALSE) {
	if (length(adjPvalMax) > 1 | !is.double(adjPvalMax) | adjPvalMax < 0 | adjPvalMax > 1) {
		stop("Input error : adjPvalMax option must be a double between 0 and 1.")
	}
	
	if (length(dPSImin) > 1 | !is.double(dPSImin) | dPSImin < 0 | dPSImin > 1) {
		stop("Input error : dPSImin option must be a double between 0 and 1. A dPSImin = 0.1 
            will catch all dPSI from -1 to -0.1 and all dPSI from 0.1 to 1.")
	}
  if(!is.logical(writePSI)) {
    stop("Input error : writePSI option must be a boolean.")
  }
  if(length(output)>1 | !is.character(output)) {
    stop("Input error : output option must be a character.")
  }
	
	k2rgFile <- resDiffExprVariant$k2rgFile
	
	if (writePSI){
		.writePSITable(resDiffExprVariant, adjPvalMax, dPSImin, output)
	} else{
		if (is.null(k2rgFile)) {
			.writeTableOutput(resDiffExprVariant$finalTable, adjPvalMax, dPSImin, 
												output)
		}
		else {
			.writeMergeOutput(resDiffExprVariant, k2rgFile, adjPvalMax, dPSImin, 
												output)
		}
	}
}