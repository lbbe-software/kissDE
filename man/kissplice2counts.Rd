\name{kissplice2counts}
\alias{kissplice2counts}
\title{Pre-treatment of KisSplice outputs.}
\description{Converts KisSplice .fa outputs to data.frame objects that can be used by others functions of the package.
}
\usage{
kissplice2counts(fileName)
}
\arguments{
  \item{fileName}{ a character vector that gives the path to .fa file that must be pre-treated.}
}
\details{
}
\value{kissplice2counts returns a data.frame object that contains several columns :
  \item{events.names}{the name of the event based on KisSplice notation}
  \item{V1}{contains the length of the event}
  \item{V2 ... to Vn}{contain the counts corresponding to the replicates of the conditions}
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

