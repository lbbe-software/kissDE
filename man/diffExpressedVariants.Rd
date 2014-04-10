\name{diffExpressedVariants}
\alias{diffExpressedVariants}
\title{Retrieve condition-specific variants in RNA-seq data.}
\description{Retrieve condition-specific variants in RNA-seq data. 
}
\usage{diffExpressedVariants(countsData, conditions)}
\arguments{
	\item{countsData}{a data.frame containing : in its first column the names of the events, in its second column the lengths of the events, in the following columns, the counts corresponding to each pair of condition/replicate}
	\item{conditions}{a character vector that gives the conditions, it has as much elements as there are pairs of condition/replicates}
}
\details{
}
\value{
	code{diffExpressedVariants} returns a data.table with :
	\item{ID}{the variation ID}
	\item{Length}{the size of the first variant}
	\item{UP_Condn_Rm_Norm}{(resp LP_Condn_Rm_Norm) returns the normalised counts for the first variant (UP, resp. second variant : LP), for the condition n (Condn) and the replicate m (Rm)}
	\item{pALLGlobalPhi.glm}{adjusted pvalue}
	\item{Deltaf/DeltaPSI}{difference of relative abundance of variants across conditions}
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
