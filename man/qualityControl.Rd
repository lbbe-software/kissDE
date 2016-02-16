\name{qualityControl}
\alias{qualityControl}

\title{Some helpful plots about the experimental plan}

\description{Function that ensures that the experimental plan passed to diffExpressedEvents has no error, by providing some helpful plots about the data. This function is meant to be used before launching diffExpressedEvents, to validate the data.}

\usage{
  qualityControl(countsData, conditions, storeFigs = FALSE)
}

\arguments{
  \item{countsData}{a data.frame containing: in its first column, the lengths of the events, in the following columns, the counts corresponding to each pair of condition/replicate}
  
  \item{conditions}{a character vector that gives the conditions, it has as much elements as there are pairs of condition/replicates}
  
	\item{storeFigs}{a logical or a string indicating if the plots should be stored and in which directory.
	By default (\code{FALSE}), the figures are not stored but plotted one at time when Enter is pressed.
	If \code{storeFigs} is \code{TRUE}, the figures are stored in a \code{kissDEFigures} directory which is created in the current directory.
	If \code{storeFigs} is a path (string, e.g. '/path/to/figs'), a new directory called with this value is created to store the figures.
	Plots are stored in .png format. 
	If the \code{diffExpressedVariants} function is a part of an automatised worflow, we recommand to set this option to \code{TRUE} or to a personal value.}
}

% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
% }

