\name{qualityControl}
\alias{qualityControl}

\title{Some helpful plots about the experimental plan}
\description{Function that ensures that the experimental plan passed to diffExpressedEvents has no error, by providing some helpful plots about the data. This function is meant to be used before launching diffExpressedEvents, to validate the data.}
\usage{
qualityControl(countsData, conditions, storeFigs=FALSE)
}
\arguments{
    \item{countsData}{a data.frame containing : in its first column, the lengths of the events, in the following columns, the counts corresponding to each pair of condition/replicate}
    \item{conditions}{a character vector that gives the conditions, it has as much elements as there are pairs of condition/replicates}
    \item{storeFigs=FALSE}{if this function is needed to be used into an automatised worflow, setting this option to TRUE will create a directory where plots are stored in .png format. By default (FALSE) the figures are plotted one at time when Enter is pressed. }
    \item{pathFigs="None"}{to give a path (string, e.g. '/path/to/figs') where to store the figures. By default the figure directory is called 'kissDEFigures' and created in the currend directory.}
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

