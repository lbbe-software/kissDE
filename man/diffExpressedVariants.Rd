\name{diffExpressedVariants}
\alias{diffExpressedVariants}

\title{Retrieve condition-specific variants in RNA-seq data}

\description{Function that retrieves condition-specific variants in RNA-seq data.}

\usage{
  diffExpressedVariants(countsData, conditions, storeFigs = FALSE, pvalue = 0.05, 
  filterLowCountsVariants = 10, flagLowCountsConditions = 10, discoSNP = FALSE)
  }

\arguments{
	\item{countsData}{the output of the \code{kissplice2counts} function or a data frame containing in its first column the names of the events, in its second column the lengths of the events, and in the following columns, the counts corresponding to each replicate of each experimental condition of one variant. Each row corresponds with one variant.}
	
	\item{conditions}{a character vector containing the experimental conditions, of length \code{NCR} (see Details below).}
	
	\item{storeFigs}{a logical or a string indicating if the plots should be stored and in which directory.
	By default (\code{FALSE}), the figures are not stored but plotted one at time when Enter is pressed.
	If \code{storeFigs} is \code{TRUE}, the figures are stored in a \code{kissDEFigures} directory which is created in the current directory.
	If \code{storeFigs} is a path (string, e.g. '/path/to/figs'), a new directory called with this value is created to store the figures.
	Plots are stored in .png format. 
	If the \code{diffExpressedVariants} function is a part of an automatised worflow, we recommand to set this option to \code{TRUE} or to a personal value.}

	\item{pvalue}{a numerical value indicating the threshold below which the pvalues are kept.}
	
	\item{filterLowCountsVariants}{a numerical value indicating the value below which events are filtered out. After having computing the overdispersion, we screen out events which both variants have global counts < \code{filterLowCountsVariants} (sum of each replicate x condition) in order to gain strength for the models.}
	
	\item{flagLowCountsConditions}{a numerical value indicating the value below which we flag events as 'lowCounts' in the final data frame. This column aims to help seeing which event (variant 1 counts + variant 2 counts) has low counts (< \code{flagLowCountsBoundary = 10}) for at least n-1 conditions (over n).}
	
	\item{discoSNP}{a logical indicating if you are using discoSNP++'s output file instead of \code{KisSplice}'s file}
	
	\item{output}{a character indicating the name of the file to output the results.}
	}
	
\value{
  \code{diffExpressedVariants} returns a list of 4 objects:
    \item{finalTable}{a data frame containing the columns
      \itemize{
        \item \code{ID}: the variation ID
        \item \code{Length_diff}: the size of the first variant
        \item \code{UP_Condn_Rm_Norm (resp LP_Condn_Rm_Norm)}: returns the normalized counts for the first variant (UP, resp. second variant : LP), for the condition n (Condn) and the replicate m (Rm)
        \item \code{Adjusted_pvalue}: pvalue adjusted for multiple testing with FDR = 5\%
        \item \code{Deltaf/DeltaPSI}: difference of relative abundance of variants across conditions. For instance if there are 2 conditions, \code{deltaPsi} returns relative abudance in condition 2 - relative abundance in condition 1. If \code{counts} is at 1 or 2, inclusion variants are corrected for the length of the variant so that we do not overestimate the psi value. If \code{counts} is at 0, we correct the longer variant by diving by this correcting factor: (longer_variant_length + readLength - 2 *overlap + 1)/(lower_variant_lenght + readLength - 2 * overlap + 1) with readLength and overlap coming from the input options.
        \item \code{lowcounts}: a column to highlight low counts in data. If 1, at least one replicate shows less than 10 reads, 0 else.
	}}
	\item{correctedPval}{a numeric vector containing pvalues after correction for multiple testing}
	\item{uncorrectedPVal}{a numeric vector containing pvalues before correction for multiple testing}
	\item{resultFitNBglmModel}{a data frame containing the results of the fitting of models to the data}
	\item{psiTable}{a data frame containing the Percent Spliced In (PSI) of each replicates}
	
	\code{diffExpressedVariants} also returns a graphic representing the variance by the mean of event counts and the Poisson, Quasi-Poisson and Negative Binomial models.
}

\details{
  \code{NCR} is defined as the number of combinaisons Condition/Replicate. It is the length of the experimental design.
  If we define \eqn{C_i} the \eqn{i^{th}} experimental condition and \eqn{R_{ij}} the \eqn{j^{th}} replicate of the \eqn{i^{th}} experimental condition, \code{Ncr} is equal to \eqn{\sum_{i}\sum_{j} \left | C_iR_{ij} \right |}.
  For example, a balanced experimental design defined by 2 replicates in the first condition and 2 replicates in the second one has a length of 4 (\code{NCR} = 4). 
  An unbalanced experimental design defined by 2 replicates in the first condition and 3 replicates in the second one has a length of 5 (\code{NCR} = 5). 
}

\examples{
fpath1 <- system.file("extdata", "output_kissplice_SNP.fa", package = "kissDE")
mySNPcounts <- kissplice2counts(fpath1, pairedEnd = TRUE)
mySNPconditions <- c("EUR", "EUR", "TSC", "TSC")
# diffSNP <- diffExpressedVariants(mySNPcounts, mySNPconditions)

fpath2 <- system.file("extdata", "table_counts_alt_splicing.txt", package = "kissDE")
mySplicingconditions <- c("C1", "C1", "C2", "C2")
mySplicingcounts <- read.table(fpath2, header = TRUE)
# diffSplicing <- diffExpressedVariants(mySplicingcounts, mySplicingconditions)
}
