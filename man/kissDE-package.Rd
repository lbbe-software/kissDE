\name{kissDE-package}
\alias{kissDE-package}
\alias{kissDE}
\docType{package}
\title{Retrieves condition-specific variants in RNA-seq data.}
\description{Retrieves condition-specific variants in RNA-seq data. Each variation (SNP, indel, alternative splicing) is represented as a pair of variants. The quantification of each variant is summarised as a count, in each condition and each replicate where it was measured. The package tests for enrichment of a variant in a condition. Data counts are modelled using either a poisson or a negative binomial (with several ways of estimating the variance). Likelihood ratio tests are then performed using the glm framework.
}
\details{
\tabular{ll}{
Package: \tab kissDE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-03\cr
License: \tab GPL\cr
}
}

\author{
Lilia Brinza, Vincent Lacroix, Frank Picard, Vincent Miele and Camille Marchet\cr
Maintainers: vincent.miele@univ-lyon1.fr, camille.marchet@inria.fr
}
\keyword{ package }

