\name{kissDE-package}
\alias{kissDE-package}
\alias{kissDE}
\docType{package}

\title{Retrieves condition-specific variants in RNA-seq data}

\description{Retrieves condition-specific variants in RNA-seq data. Each variation (SNP, indel, alternative splicing) is represented as a pair of variants. The quantification of each variant is summarised as a count, in each condition and each replicate where it was measured. The package tests for enrichment of a variant in a condition. Data counts are modelled using either a poisson or a negative binomial (with several ways of estimating the variance). Likelihood ratio tests are then performed using the glm framework.
}

\details{Main functions:

  diffExpressedVariants(countsData, conditions, storeFigs = FALSE, pvalue = 0.05, 
  filterLowCountsVariants = 10, flagLowCountsConditions = 10, discoSNP = FALSE)

  qualityControl(countsData, conditions, storeFigs = FALSE)
  
  kissplice2counts(fileName, counts = 0, pairedEnd = FALSE, order = NULL,
  exonicReads = TRUE, discoSNP = FALSE, k2rg = FALSE)
  
  \tabular{ll}{
    Package: \tab \code{kissDE}\cr
    Type: \tab Package\cr
    Version: \tab 1.2.0\cr
    Date: \tab 2016-02-16\cr
    License: \tab GPL (>=2)\cr
  }
}

\author{
Lilia Brinza, Vincent Lacroix, Frank Picard, Vincent Miele and Camille Marchet\cr
Maintainers: marchetcamille@gmail.com
}

\keyword{package}
