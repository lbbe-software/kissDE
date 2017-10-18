\name{writeOutputKissDE}
\alias{writeOutputKissDE}

\title{Create and store the output of the \code{diffExpressedVariants} function 
in a readable file.}


\description{If a \code{KisSplice} fasta file was used in input of the 
analysis, \code{writeOutputKissDE} will output a tab-delimited table
containing one alternative splicing event/SNV per line. The fields are: the ID
of the variation, the variable part length, the counts of each variant for 
each conditions, the adjusted p-value, the deltaPSI and a boolean indicating 
if the splicing event/SNV was correctly expressed in at least half of the
conditions.

If a \code{KisSplice2RefGenome} file was used in input of the analysis,
this function will add five columns to the \code{KisSplice2RefGenome}
file, containing the following \code{KissDE} results: normalized counts,
normalized PSI, adjusted p-value, deltaPSI and a boolean indicating if the
splicing event/SNV was correctly expressed in at least half of the conditions.

\code{writePSI} must be set to \code{TRUE} to write the table containing the 
f/PSI values.
}


\usage{
writeOutputKissDE(resDiffExprVariant, adjPvalMax = 1, dPSImin = 0, output, 
    writePSI = FALSE)
}


\arguments{
\item{resDiffExprVariant}{a list, returned by \code{diffExpressedVariants}.}

\item{adjPvalMax}{a double indicating the threshold for adjusted p-value. Only
SNVs/splicing events with an adjusted p-value lower than this threshold will be
written in the output file.}

\item{dPSImin}{a double indicating the threshold for the deltaPSI. Only
SNVs/splicing events which has the absolute value of their deltaf/deltaPSI
higher than this threshold will be written in the output file.}

\item{output}{a character indicating the name of the file to save 
\code{writeOutputKissDE} output.}

\item{writePSI}{a boolean indicating if the user wants the f/PSI table to be
printed instead of the final table.}
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
#	    conditions=mySplicingconditions)
# writeOutputKissDE(res, output="results.tsv")
# writeOutputKissDE(res, output="significants_results.tsv", 
#			adjPvalMax=0.05, dPSImin=0.1)
# writeOutputKissDE(res, output="psi_results.tsv", adjPvalMax=0.05,
#     dPSImin=0.1, writePSI=TRUE)
}