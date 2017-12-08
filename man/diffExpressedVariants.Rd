\name{diffExpressedVariants}
\alias{diffExpressedVariants}

\title{Retrieve condition-specific variants in RNA-seq data}

\description{Function that retrieves condition-specific variants in RNA-seq 
data.}


\usage{
diffExpressedVariants(countsData, conditions, pvalue = 1, 
    filterLowCountsVariants = 10, flagLowCountsConditions = 10,
    technicalReplicates = FALSE)
}


\arguments{
\item{countsData}{a data frame containing the counts in the appropriate format
(see Details below).}

\item{conditions}{a character vector containing the experimental conditions.}

\item{pvalue}{a numerical value indicating the p-value threshold below which
the events will be kept in the final data frame.}

\item{filterLowCountsVariants}{a numerical value indicating the global variant 
count value (see Details below) below which events are filtered out in order 
to increase statistical power of the analysis. Both variant must have a read
coverage below this value in order to remove the event. This filter is done 
after the normalization and the overdispersion estimation.}

\item{flagLowCountsConditions}{a numerical value indicating the global 
condition count value (see Details below) below which we flag events as 
'lowCounts' in the final data frame. At least n-1 conditions (over n 
conditions) must have low counts to flag the event as 'lowCounts'.}

\item{technicalReplicates}{a boolean value indicating if the counts in 
\code{countsData} come from technical replicates only or not.}
}


\value{
\code{\link{diffExpressedVariants}} returns a list of 6 objects:
\item{finalTable}{a data frame containing the columns
    \itemize{
        \item \code{ID}: the variation identifier
        \item \code{Length_diff}: the size of the variable region
        \item \code{UP_Condi_Rj_Norm (resp LP_Condi_Rj_Norm)}: returns the
        normalized counts of the first variant (UP, resp. second variant: LP),
        for the condition i (\code{Condi}) and the replicate j (\code{Rj})
        \item \code{Adjusted_pvalue}: p-value adjusted for multiple testing with
        Benjamini & Hochberg method
        \item \code{Deltaf/DeltaPSI}: difference of relative abundance of
        variants across conditions. For instance if there are 2 conditions, 
        \code{deltaPSI} returns relative abudance in condition 2 - relative
        abundance in condition 1. Inclusion variant's counts are corrected for
        the length of the variant so that we do not overestimate the PSI value.
        \item \code{lowcounts}: a column that flag low counts in data. 
        If \code{TRUE}, at least n-1 conditions over n conditions have less 
        than 10 reads.
}}
\item{correctedPval}{a numeric vector containing p-values after correction 
for multiple testing}
\item{uncorrectedPVal}{a numeric vector containing p-values before correction
for multiple testing}
\item{resultFitNBglmModel}{a data frame containing the results of the fitting
of the model to the data}
\item{f/psiTable}{a data frame containing the allele frequency (f)/Percent 
Spliced In (PSI) of each replicate}
\item{k2rgFile}{a string containing either the \code{KisSplice2RefGenome} file
path and name or NULL if no \code{KisSplice2RefGenome} input file was given}
}


\details{
The \code{countsData} data frame must be formatted as follows: 
    \itemize{
        \item Column 1: names of the events
        \item Column 2: lengths (in bp) of the variants
        \item Column 3 to n: counts corresponding to each replicate of each 
            experimental condition of one variant
    }
Each row corresponds to one variant, thus an event correspond to two rows with
the longest variant (or inclusion variant) in the first row (thus denotated as 
upper path: UP) and the smallest variant (or exclusion variant) in the second 
row (thus denotated as lower path: LP).
This data frame can be obtained using \code{\link{kissplice2counts}}
function.\\
The global variant count is the minimal number of reads that cover one or the 
other variant across all the replicates (sum by variant).\\
The global condition count is the minimal number of reads that cover one or 
the other condition (sum by replicates for each conditions).}

\references{
Lopez-Maestre et al., 2016. Snp calling from rna-seq data without a reference 
genome: identification, quantification, differential analysis and impact on the 
protein sequence. \emph{Nucleic Acids Research}, \bold{44}(19):e148. 
\url{https://doi.org/10.1093/nar/gkw655}
}

\examples{
fpath1 <- system.file("extdata", "output_kissplice_SNV.fa", package = "kissDE")
mySNVcounts <- kissplice2counts(fpath1, pairedEnd = TRUE)
mySNVconditions <- c("EUR", "EUR", "TSC", "TSC")
# diffSNV <- diffExpressedVariants(mySNVcounts, mySNVconditions)

fpath2 <- system.file("extdata", "table_counts_alt_splicing.txt", 
package = "kissDE")
mySplicingconditions <- c("C1", "C1", "C2", "C2")
mySplicingcounts <- read.table(fpath2, header = TRUE)
# diffSplicing <- diffExpressedVariants(mySplicingcounts, mySplicingconditions)
}