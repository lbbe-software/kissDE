\name{diffExpressedVariants}
\alias{diffExpressedVariants}
\title{Retrieve condition-specific variants in RNA-seq data.}
\description{Retrieve condition-specific variants in RNA-seq data. 
}
\usage{diffExpressedVariants(countsData, conditions, storeFigs=FALSE, pvalue=0.05)}
\arguments{
	\item{countsData}{a data.frame containing : in its first column the names of the events, in its second column the lengths of the events, in the following columns, the counts corresponding to each pair of condition/replicate}
	\item{conditions}{a character vector that gives the conditions, it has as much elements as there are pairs of condition/replicates}
	\item{storeFigs=FALSE}{if this function is needed to be used into an automatised worflow, setting this option to TRUE will create a Figures directory where plots are stored in .png format. By default (FALSE) the figures are plotted one at time when Enter is pressed. }
	\item{pvalue=0.05}{pvalue to select significant variants can be set, its value is 0.05 by default.}
}
\details{
}
\value{
	code{diffExpressedVariants} returns a data.table with :
	\item{ID}{the variation ID}
	\item{Length}{the size of the first variant}
	\item{UP_Condn_Rm_Norm}{(resp LP_Condn_Rm_Norm) returns the normalised counts for the first variant (UP, resp. second variant : LP), for the condition n (Condn) and the replicate m (Rm)}
	\item{pALLGlobalPhi.glm}{adjusted pvalue}
	\item{Deltaf/DeltaPSI}{difference of relative abundance of variants across conditions. For instance if there are 2 conditions, deltaPsi returns relative abudance in condition 2 - relative abundance in condition 1.}
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
