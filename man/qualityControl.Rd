\name{qualityControl}
\alias{qualityControl}

\title{Some helpful plots about the experimental plan}
\description{Function that ensures that the experimental plan passed to diffExpressedEvents has no error, by providing some helpful plots about the data. This function is meant to be used before launching diffExpressedEvents, to validate the data.}
\usage{
qualityControl(countsData, conditions)
}
\arguments{
    \item{countsData}{a data.frame containing : in its first column, the lengths of the events, in the following columns, the counts corresponding to each pair of condition/replicate}
    \item{conditions}{a character vector that gives the conditions, it has as much elements as there are pairs of condition/replicates}
}
\details{
}
\value{
}
\references{
}
\author{
}
\note{
}
\seealso{
}
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
% }

