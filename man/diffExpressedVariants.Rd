\name{diffExpressedVariants}
\alias{diffExpressedVariants}

\title{Retrieve condition-specific variants in RNA-seq data}

\description{Function that retrieves condition-specific variants in RNA-seq data.}

\usage{
  diffExpressedVariants(countsData, conditions, storeFigs = FALSE, pvalue = 0.05, 
  filterLowCountsVariants = 10, flagLowCountsConditions = 10, discoSNP = FALSE)
  }

\arguments{
	\item{countsData}{a data.frame containing: in its first column the names of the events, in its second column the lengths of the events, in the following columns, the counts corresponding to each pair of condition/replicate}
	
	\item{conditions}{a character vector that gives the conditions, it has as much elements as there are pairs of condition/replicates}
	
	\item{storeFigs}{a logical or a string indicating if the plots should be stored and in which directory.
	By default (\code{FALSE}), the figures are not stored but plotted one at time when Enter is pressed.
	If \code{storeFigs} is \code{TRUE}, the figures are stored in a \code{kissDEFigures} directory which is created in the current directory.
	If \code{storeFigs} is a path (string, e.g. '/path/to/figs'), a new directory called with this value is created to store the figures.
	Plots are stored in .png format. 
	If the \code{diffExpressedVariants} function is a part of an automatised worflow, we recommand to set this option to \code{TRUE} or to a personal value.}

	\item{pvalue}{pvalue}
	
	\item{filterLowCountsVariants}{to change the default value under which events are filtered out. After having computing the overdispersion, we screen out events which both variants have global counts < [filterBoundary] (sum of each replicate x condition) in order to gain strength for the models.}
	
	\item{flagLowCountsConditions}{to change the default value under which we flag events as 'lowCounts' in the final data frame. This column aims to help seeing which event (variant 1 counts + variant 2 counts) has low counts (< [flagLowCountsBoundary = 10]) for at least n-1 conditions (over n).}
	
	\item{discoSNP}{put to TRUE if you are using discoSNP++'s output file instead of \code{KisSplice}'s file}
	}
	
\value{
	\code{diffExpressedVariants} returns several data.table:

	\item{finalTable}{with
	\itemize{
    \item \code{ID}: the variation ID
    \item \code{Length}: the size of the first variant
    \item \code{UP_Condn_Rm_Norm (resp LP_Condn_Rm_Norm)}: returns the normalised counts for the first variant (UP, resp. second variant : LP), for the condition n (Condn) and the replicate m (Rm)
    \item \code{Adjusted_pvalue}: pvalue adjusted for multiple testing with FDR = 5\%
    \item \code{Deltaf/DeltaPSI}: difference of relative abundance of variants across conditions. For instance if there are 2 conditions, deltaPsi returns relative abudance in condition 2 - relative abundance in condition 1. If \code{counts} is at 1 or 2, inclusion variants are corrected for the length of the variant so that we do not overestimate the psi value. If \code{counts} is at 0, we correct the longer variant by diving by this correcting factor: (longer_variant_length + readLength - 2 *overlap + 1)/(lower_variant_lenght + readLength - 2 * overlap + 1) with readLength and overlap coming from the input options.
    \item \code{Low_counts}: a column to highlight low counts in data. If 1, at least one replicate shows less than 10 reads, 0 else.
	}}
	
	\item{uncorrectedPVal}{with pvalues before correction for multiple testing}
	
	\item{correctedPval}{with pvalues after correction for multiple testing}
	
	\item{resultFitNBglmModel}{with the results of the fitting of models to the data}
}

% \examples{
% conditions <- c("treated","treated","untreated","untreated")
% file <- ("myresults.fa")
% counts <- kissplice2counts(file)
% events <- diffExpressedVariants(counts, conditions)

% % }
