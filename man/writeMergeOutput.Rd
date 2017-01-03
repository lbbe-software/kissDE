\name{writeMergeOutput}
\alias{writeMergeOutput}

\title{Write an output file combining kissplice2refgenome and \code{KissDE} results.}

\description{Function that add five columns to the kissplice2refgenome's file, used in \code{KissDE}, containing the following \code{KissDE} results : normalized counts, normalized psi, adjusted p-value, delta-PSI, warnings. If the analysis with \code{KissDE} was made without a kissplice2refgenome file, then no output is created.
}

\usage{
  writeMergeOutput(finalTable, psiTable, output, k2rgFile)
}

\arguments{
  \item{finalTable}{a list, returned with the name \code{finalTable} as a result of \code{diffExpressedVariants}.}
  
  \item{psiTable}{a list, returned with the name \code{`f/psiTable`} as a result of \code{diffExpressedVariants}.}
  
  \item{output}{a character indicating the name of the file to save \code{writeMergeOutput} output.}
  
  \item{k2rgFile}{a character indicating the name of the kissplice2refgenome file to add columns in.}
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
  # writeMergOutput(finalTable = res$finalTable, psiTable = res$`f/psiTable`, 
  #  output = "path/to/my/output/file", k2rgFile = kissplice2refgenome_file)
}
  
  

