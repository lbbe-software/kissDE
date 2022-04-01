\name{exploreResults}
\alias{exploreResults}

\title{Print kissDE results in an interactive Shiny application}


\description{This function will read the rds file created by 
\code{\link{writeOutputKissDE}} and creat an interactive Shiny 
application allowing users to explore and plot the results of kissDE.}


\usage{
exploreResults(rdsFile)
}


\arguments{
\item{rdsFile}{a string indicating the path to the rds file outputed by
\code{\link{writeOutputKissDE}}.}
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
# exploreResults("results.tsv.rds")
}
