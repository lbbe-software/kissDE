\name{qualityControl}
\alias{qualityControl}

\title{Quality control plots}


\description{Function that ensures that the experimental plan passed to 
\code{diffExpressedVariants} has no error, by providing some helpful plots 
about the data. This function is meant to be used before launching 
\code{diffExpressedVariants}, to validate the data.}


\usage{
qualityControl(countsData, conditions, storeFigs = FALSE)
}


\arguments{
\item{countsData}{the output of the \code{kissplice2counts} function or a data
frame containing the counts in the appropriate format (see Details below).}

\item{conditions}{a character vector that gives the conditions. It has as 
many elements as there are samples in total.}

\item{storeFigs}{a logical or a string indicating if the plots should be 
stored and in which directory.	By default (\code{FALSE}), the figures are not
stored but plotted one at a time when Enter is pressed. If \code{storeFigs} is
\code{TRUE}, the figures are stored in a \code{kissDEFigures} directory which
is created in the working directory. If \code{storeFigs} is a path (string, 
e.g. '/path/to/figs'), this directory is created to store the figures. Plots 
are stored in .png format. If the \code{qualityControl} function is a part 
of an automatised worflow, we recommand to set this option to \code{TRUE} or 
to a user defined value.}
}


\value{
The figures are stored or displayed in the R session.
}

\details{
\code{countsData} input must be formatted as follows: in its first column the 
names of the events, in its second column the lengths of the events, and in 
the following columns, the counts corresponding to each replicate of each 
experimental condition of one variant. Each row corresponds to one variant.
}

\seealso{
\code{\link{diffExpressedVariants}}
}

\examples{
fpath <- system.file("extdata", "output_kissplice_SNV.fa", package="kissDE")
mySNVcounts <- kissplice2counts(fpath, pairedEnd=TRUE)
mySNVconditions <- c("EUR", "EUR", "TSC", "TSC")
qualityControl(mySNVcounts, mySNVconditions)
}