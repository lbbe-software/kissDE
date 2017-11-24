\name{qualityControl}
\alias{qualityControl}

\title{Quality control plots}


\description{Function that provide some helpful plots about the data to
ensure that the experimental plan passed to \code{diffExpressedVariants} 
is correct and to control the quality of the data.
This function should be run before launching \code{diffExpressedVariants}, 
to validate the data.}


\usage{
qualityControl(countsData, conditions, storeFigs = FALSE, returnPCAdata = FALSE)
}


\arguments{
\item{countsData}{the output of the \code{kissplice2counts} function or a data
frame containing the counts in the appropriate format (see Details below).}

\item{conditions}{a character vector that gives the conditions' order. It has as 
many elements as there are samples in total.}

\item{storeFigs}{a logical or a string indicating if the plots should be 
stored and in which directory. If the \code{qualityControl} function is a part 
of an automatised worflow, we recommand to set this option to \code{TRUE} (or 
to a user defined value). Thus, the figures will be stored in a 
\code{kissDEFigures} directory which is created in the working directory. If 
\code{storeFigs} is a path (string, e.g. '/path/to/figs'), this directory is created 
to store the figures. Plots are stored in .png format. By default (\code{FALSE}), 
the figures are not stored but plotted one at a time.}

\item{returnPCAdata}{a logical indicating if the data frame used in the PCA
analysis should be returned. By default (\code{FALSE}), the data frame is not
returned.}
}


\value{
The figures are saved or displayed in the R session.
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