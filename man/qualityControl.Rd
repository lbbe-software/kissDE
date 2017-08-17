\name{qualityControl}
\alias{qualityControl}

\title{Some helpful plots about the experimental plan}


\description{Function that ensures that the experimental plan passed to 
\code{diffExpressedVariants} has no error, by providing some helpful plots 
about the data. This function is meant to be used before launching 
\code{diffExpressedVariants}, to validate the data.}


\usage{
qualityControl(countsData, conditions, storeFigs = FALSE)
}


\arguments{
\item{countsData}{the output of the \code{kissplice2counts} function or a data
frame containing in its first column the names of the events, in its second
column the lengths of the events, and in the following columns, the	counts
corresponding to each pair of conditions/replicates}

\item{conditions}{a character vector that gives the conditions, it has as 
many elements as there are replicates in total.}

\item{storeFigs}{a logical or a string indicating if the plots should be 
stored and in which directory.	By default (\code{FALSE}), the figures are not
stored but plotted one at time when Enter is pressed. If \code{storeFigs} is
\code{TRUE}, the figures are stored in a \code{kissDEFigures} directory which
is created in the current directory. If \code{storeFigs} is a path (string, 
e.g. '/path/to/figs'), a new directory called with this value is created to
store the figures. Plots are stored in .png format. If the 
\code{qualityControl} function is a part of an automatised worflow, we 
recommand to set this option to \code{TRUE} or to a personal value.}
}


\value{
The figures are stored or displayed in the R session.
}


\examples{
fpath <- system.file("extdata", "output_kissplice_SNV.fa", package="kissDE")
mySNVcounts <- kissplice2counts(fpath, pairedEnd=TRUE)
mySNVconditions <- c("EUR", "EUR", "TSC", "TSC")
qualityControl(mySNVcounts, mySNVconditions)
}