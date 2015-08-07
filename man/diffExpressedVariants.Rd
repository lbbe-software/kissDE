\name{diffExpressedVariants}
\alias{diffExpressedVariants}
\title{Retrieve condition-specific variants in RNA-seq data.}
\description{Retrieve condition-specific variants in RNA-seq data. 
}
\usage{diffExpressedVariants(countsData, conditions, storeFigs=FALSE, pathFigs="None", pvalue=0.05, filterLowCountsVariants=10, flagLowCountsConditions=10, readLength=75, overlap=42, discoSNP=FALSE)}
\arguments{
	\item{countsData}{a data.frame containing : in its first column the names of the events, in its second column the lengths of the events, in the following columns, the counts corresponding to each pair of condition/replicate}
	\item{conditions}{a character vector that gives the conditions, it has as much elements as there are pairs of condition/replicates}
	\item{storeFigs=FALSE}{if this function is needed to be used into an automatised worflow, setting this option to TRUE will create a directory where plots are stored in .png format. By default (FALSE) the figures are plotted one at time when Enter is pressed. }
	\item{pathFigs="None"}{to give a path (string, e.g. '/path/to/figs') where to store the figures. By default the figure directory is called 'kissDEFigures' and created in the currend directory.}
	\item{filterLowCountsVariants=10}{to change the default value under which events are filtered out. After having computing the overdispersion, we screen out events which both variants have global counts < [filterBoundary](sum of each replicate x condition) in order to gain strength for the models.}
	\item{flagLowCountsConditions=10}{to change the default value under which we flag events as 'lowCounts' in the final data frame. This column aims to help seeing which event (variant 1 counts + variant 2 counts) has low counts (< [flagLowCountsBoundary=10])for at least n-1 conditions (over n).}
	\item{readLength=75}{length of the reads from the data set used. This size is used to correct the delta PSI computation.}
	\item{overlap=42}{mapping seed length, usually k+1 (with k size of the k-mers in KisSplice). The value is used to correct the delta PSI computation.}
	}
	\item{discoSNP=FALSE}{ put to TRUE if you are using discoSNP++'s output file instead of KisSplice's file}

\details{
}
\value{
	code{diffExpressedVariants} returns several data.table. 
	finalTable with :
	\item{ID}{the variation ID}
	\item{Length}{the size of the first variant}
	\item{UP_Condn_Rm_Norm}{(resp LP_Condn_Rm_Norm) returns the normalised counts for the first variant (UP, resp. second variant : LP), for the condition n (Condn) and the replicate m (Rm)}
	\item{Adjusted_pvalue}{pvalue adjusted for multiple testing with FDR=5\%}
	\item{Deltaf/DeltaPSI}{difference of relative abundance of variants across conditions. For instance if there are 2 conditions, deltaPsi returns relative abudance in condition 2 - relative abundance in condition 1. If --counts is at 1 or 2, inclusion variants are corrected for the length of the variant so that we do not overestimate the psi value. If --counts is at 0, we correct the longer variant by diving by this correcting factor: (longer_variant_length + readLength - 2*overlap + 1)/(lower_variant_lenght+ readLength - 2*overlap + 1) with readLength and overlap coming from the input options.}
	\item{Low_counts}{a column to highlight low counts in data. If 1, at least one replicate shows less than 10 reads, 0 else.}
	uncorrectedPVal with pvalues before correction for multiple testing, correctedPval with pvalues after correction for multiple testing, and resultFitNBglmModel with the results of the fitting of models to the data.

}
\source{
}
\author{
}


\references{
}
% \examples{
% conditions <- c("treated","treated","untreated","untreated")
% file <- ("myresults.fa")
% counts <- kissplice2counts(file)
% events <- diffExpressedVariants(counts, conditions)

% % }
\keyword{}
