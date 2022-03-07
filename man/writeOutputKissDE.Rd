\name{writeOutputKissDE}
\alias{writeOutputKissDE}

\title{Create and store the output of the \code{\link{diffExpressedVariants}}
function in a file and in a rds object.}


\description{If a \code{KisSplice} fasta file was used as input for the 
analysis, \code{\link{writeOutputKissDE}} will output a tab-delimited file
containing one alternative splicing event/SNV per line. The columns are: the 
ID of the variation, the variable part length, the counts of each variant for 
each condition, the adjusted p-value (FDR), the deltaPSI and a boolean 
indicating if the splicing event/SNV was sufficiently expressed (controled by 
the \code{flagLowCountsConditions} option from the
\code{\link{diffExpressedVariants}} function).

If a \code{KisSplice2RefGenome} file was used as input for the analysis,
this function will add five columns to the \code{KisSplice2RefGenome}
file, with the following \code{KissDE} results: normalized counts,
PSI computed from normalized counts, adjusted p-value, deltaPSI and a boolean 
indicating if the splicing event/SNV was sufficiently expressed in at least 
half of the conditions (controled by the \code{flagLowCountsConditions} option 
from the \code{\link{diffExpressedVariants}} function).

In both cases, an rds object is saved in the output folder, that can be used to 
explore the results of kissDE through a Shiny application with the 
\code{\link{exploreResults}} function.
}


\usage{
writeOutputKissDE(resDiffExprVariant, output, adjPvalMax = 1, dPSImin = 0, 
    writePSI = FALSE)
}


\arguments{
\item{resDiffExprVariant}{a list, returned by 
\code{\link{diffExpressedVariants}}.}

\item{output}{a character indicating the path and file name to save 
\code{\link{writeOutputKissDE}} output.}

\item{adjPvalMax}{a double indicating the threshold for adjusted p-value. Only
SNVs/splicing events with an adjusted p-value lower than this threshold will be
kept in the output file.}

\item{dPSImin}{a double indicating the threshold for the deltaPSI. Only
SNVs/splicing events having an absolute value of deltaf/deltaPSI
higher than this threshold will be kept in the output file.}

\item{writePSI}{a boolean indicating if the user wants the f/PSI table to be
printed (\code{TRUE}) instead of the final table (\code{FALSE}, default).}
}


\value{
None.
}


\examples{
kissplice2refgenome_file <- system.file("extdata", 
    "output_k2rg_alt_splicing.txt", package="kissDE")
mySplicingconditions <- c("C1", "C1", "C2", "C2")
counts <- kissplice2counts(fileName=kissplice2refgenome_file, counts=2,
    pairedEnd=TRUE, k2rg=TRUE)
# res <- diffExpressedVariants(countsData=counts, 
#   conditions=mySplicingconditions)
# writeOutputKissDE(res, output="results.tsv")
# writeOutputKissDE(res, output="significants_results.tsv", 
#   adjPvalMax=0.05, dPSImin=0.1)
# writeOutputKissDE(res, output="psi_results.tsv", adjPvalMax=0.05,
#   dPSImin=0.1, writePSI=TRUE)
}