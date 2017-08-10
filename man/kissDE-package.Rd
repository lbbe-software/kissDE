\name{kissDE-package}
\alias{kissDE-package}
\alias{kissDE}
\docType{package}

\title{Retrieves condition-specific variants in RNA-seq data}

\description{Retrieves condition-specific variants in RNA-seq data. Each 
variation (SNPs, alternative-splicings, indels) is represented as a pair of
variants. The quantification of each variant is summarized as a count, in each
condition and each replicate where it was measured. The package tests for
enrichment of a variant in a condition. Data counts are modelled using either a
poisson or a negative binomial (with several ways of estimating the variance).
Likelihood ratio tests are then performed using the GLM (Genralized Linear 
Model) framework.
}

\details{Main functions:

  diffExpressedVariants(countsData, conditions, storeFigs=FALSE, pvalue=0.05,
		filterLowCountsVariants=10, flagLowCountsConditions=10, discoSNP=FALSE)

  qualityControl(countsData, conditions, storeFigs=FALSE)
  
  kissplice2counts(fileName, counts=0, pairedEnd=FALSE, order=NULL,
  	exonicReads=TRUE, discoSNP=FALSE, k2rg=FALSE)
}

\keyword{package}
