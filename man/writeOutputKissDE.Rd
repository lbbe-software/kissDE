\name{writeOutputKissDE}
\alias{writeOutputKissDE}

\title{Write an output file for \code{KissDE}'s \code{diffExpressedVariants} result.}

\description{If you were working with a fasta file as input for \code{kissplice2counts}, this function will output a tab-delimited table containing one alternative splicing event/SNV per line. The fields are : the ID of the path, the variable part length, the counts of each variant for each conditions, the adjusted p-value, the dPSI and a boolean indicating if the event/SNV was correctly cover by reads. If you were working with a kissplice2refgenome file, this function will add five columns to the kissplice2refgenome's file, containing the following \code{KissDE} results : normalized counts, normalized psi, adjusted p-value, delta-PSI, warnings.
}

\usage{
  writeOutputKissDE(resDiffExprVariant, adjPvalMax = 1, dPSImin = 0, output)
}

\arguments{
  \item{resDiffExprVariant}{a list, returned by \code{diffExpressedVariants}.}
  
  \item{adjPvalMax}{a double, indicating the maximum adjusted p-value allowed to output an event/SNV.}
  
  \item{dPSImin}{a double, indicating the minimum absolute value of a dPSI allowed to output an event/SNV.}
  
  \item{output}{a character indicating the name of the file to save \code{writeOutputKissDE} output.}
}

\value{
  None.
}

\examples{
  kissplice2refgenome_file <- system.file("extdata", "output_k2rg_alt_splicing.txt", 
    package = "kissDE")
  mySplicingconditions <- c("C1", "C1", "C2", "C2")
  counts <- kissplice2counts(fileName = kissplice2refgenome_file, pairedEnd = TRUE, k2rg = TRUE)
  res <- diffExpressedVariants(countsData = counts, conditions = mySplicingconditions)
  # writeMergOutput(finalTable = res, output = "path/to/my/output/file")
  # writeMergOutput(finalTable = res, output = "path/to/my/output/file",
  # adjPvalMax = 0.05, dPSImin = 0.1)
}
  
  

