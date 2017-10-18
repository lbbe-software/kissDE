\name{kissDE-package}
\alias{kissDE-package}
\alias{kissDE}
\docType{package}

\title{Retrieves condition-specific variants in RNA-seq data}

\description{The \pkg{kissDE} package retrieves condition-specific variants in 
RNA-seq data. Each variation (SNPs, alternative-splicings, indels) is 
represented as a pair of variants. The quantification of each variant is 
summarized as a count, in each condition and each replicate where it was 
measured. The package tests for enrichment of a variant in a condition. Data 
counts are modelled using either a poisson or a negative binomial (with several 
ways of estimating the variance). Likelihood ratio tests are then performed 
using the GLM (Genralized Linear Model) framework.
}


\details{Main functions:

diffExpressedVariants(countsData, conditions, pvalue = 1, 
    filterLowCountsVariants = 10, flagLowCountsConditions = 10,
    technicalReplicates = FALSE)
    
qualityControl(countsData, conditions, storeFigs = FALSE)

kissplice2counts(fileName, counts = 0, pairedEnd = FALSE, order = NULL, 
    exonicReads = TRUE, k2rg = FALSE, keep = c("All"), remove = NULL)

writeOutputKissDE(resDiffExprVariant, adjPvalMax = 1, dPSImin = 0, output, 
    writePSI = FALSE)
}


\note{
Authors of the package: 
Clara Benoit-Pilven, Camille Marchet, Janice Kielbassa, Lilia Brinza, Audric 
Cologne and Vincent Lacroix all contributed code and ideas.

Contributors of the package: 
Franck Picard and Laurent Jacob provided statistical expertise for the models 
underlying kissDE.
Vincent Miele provided expertise for the development of the R package. 

Maintainer of the package: 
Aurélie Siberchicot
}

\keyword{package}
