\name{writeOutputKissDE}
\alias{writeOutputKissDE}

\title{Create and store the output of the \code{diffExpressedVariants} function 
in a readable file.}


\description{If a fasta file was used as input in the \code{kissplice2counts} 
function which the result was used as input in the \code{diffExpressedVariants} 
function, \code{writeOutputKissDE} will output a tab-delimited table
containing one alternative splicing event/SNV per line. The fields are: the ID
of the path, the variable part length, the counts of each variant for each
conditions, the adjusted p-value, the dPSI and a boolean indicating if the 
event/SNV was correctly cover by reads.

If you were working with a \code{KisSplice2RefGenome} file, this function will
add five columns to the \code{KisSplice2RefGenome}'s file, containing the
following \code{KissDE} results: normalized counts, normalized psi, adjusted 
p-value, delta-PSI, warnings.

If you want to write the table containing the PSI/f values, you should set
\code{writePSI} to \code{TRUE}. 
}


\usage{
writeOutputKissDE(resDiffExprVariant, adjPvalMax = 1, dPSImin = 0, output, 
    writePSI = FALSE)
}


\arguments{
\item{resDiffExprVariant}{a list, returned by \code{diffExpressedVariants}.}

\item{adjPvalMax}{a double indicating the maximum adjusted p-value allowed to
output an event/SNV.}

\item{dPSImin}{a double indicating the minimum absolute value of a dPSI
allowed to output an event/SNV.}

\item{output}{a character indicating the name of the file to save 
\code{writeOutputKissDE} output.}

\item{writePSI}{a boolean indicating if the user wants the PSI/f table to be
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
# writeOutputKissDE(res, output="path/to/my/output/file")
# writeOutputKissDE(res, output="path/to/my/output/file", 
#			adjPvalMax=0.05, dPSImin=0.1)
# writeOutputKissDE(res, output="path/to/my/output/file", adjPvalMax=0.05,
#     dPSImin=0.1, writePSI=TRUE)
}