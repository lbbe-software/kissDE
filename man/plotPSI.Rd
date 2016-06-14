\name{plotPSI}
\alias{plotPSI}

\title{Plot the PSI}

\description{Function to plot the PSI for each replicates. This function is meant to be used after launching \code{diffExpressedVariants}, to visualize the significants events.}

\usage{
  plotPSI(diffVariants, conditions, thresholdPvalue = 0.05, thresholdDeltaPSI = 0.10, 
  show = "PSI", storeFigs = FALSE)
}

\arguments{
  \item{diffVariants}{the output of the \code{diffExpressedVariants} function}
  
  \item{conditions}{a character vector that gives the conditions, it has as many elements as there are replicates in total.}
  
  \item{thresholdPvalue}{a numerical value between 0 and 1. 
  Indicate the p-value ajusted under which we consider the event as significant. 
  By default (\code{0.05}), only the events with a p-value adjusted lower than 0.05 are shown as significant on the plot.}
  
  \item{thresholdDeltaPSI}{a numerical value between 0 and 1. 
  Indicate the DeltaPSI over which we consider the event as significant. 
  By default (\code{0.10}), only the events with a DeltaPSI higher than 0.10 and lower than -0.10 are shown as significant on the plot.}
  
  \item{show}{a string value with 3 differents possibilities : "PSI" (the default), "mean" and "both".
  This argument indicate what should be plot : only the value of PSI for each replicates ("PSI"), only the mean of PSI by condition ("mean") or the both ("both").}
  
	\item{storeFigs}{a logical or a string indicating if the plots should be stored and in which directory.
	By default (\code{FALSE}), the figures are not stored but plotted one at time when Enter is pressed.
	If \code{storeFigs} is \code{TRUE}, the figures are stored in a \code{kissDEFigures} directory which is created in the current directory.
	If \code{storeFigs} is a path (string, e.g. '/path/to/figs'), a new directory called with this value is created to store the figures.
	Plots are stored in .png format. 
	If the \code{plotPSI} function is a part of an automatised worflow, we recommand to set this option to \code{TRUE} or to a personal value.}
}


\value{
The figure is stored or displayed in the R session.
}


\examples{
fpath <- system.file("extdata", "table_counts_alt_splicing.txt", package = "kissDE")
mySplicingconditions <- c("C1", "C1", "C2", "C2")
mySplicingcounts <- read.table(fpath, header = TRUE)
# diffSplicing <- diffExpressedVariants(mySplicingcounts, mySplicingconditions, pvalue = 1)
# plotPSI(diffSplicing, mySplicingconditions)
}