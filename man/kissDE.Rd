\name{kissDE}
\alias{kissDE}

\title{Run the whole kissDE analysis}


\description{This function will sequentially call the 
\code{\link{kissplice2counts}}, \code{\link{diffExpressedVariants}} and 
\code{\link{writeOutputKissDE}} functions in order to run a complete 
kissDE analysis. It may also call the \code{\link{qualityControl}} and 
the \code{\link{exploreResults}} functions.}


\usage{
kissDE(fileName, conditions, output, counts = 2, pairedEnd = FALSE,
       order = NULL, exonicReads = TRUE, k2rg = FALSE, keep = c("All"), remove = NULL,
       pvalue = 1, filterLowCountsVariants = 10, flagLowCountsConditions = 10,
       technicalReplicates = FALSE, nbCore = 1, adjPvalMax = 1, dPSImin = 0,
       writePSI = TRUE, doQualityControl = TRUE, resultsInShiny=TRUE)
}


\arguments{
\item{fileName}{a string indicating the path to the \code{KisSplice} 
(\code{.fa}) or the \code{KisSplice2RefGenome} (tab-delimited) file.}

\item{conditions}{a character vector containing the experimental conditions.}
  
\item{output}{a character indicating the path and file name to save 
\code{\link{writeOutputKissDE}} output.}
  
\item{counts}{an interger (0, 1 or 2) corresponding to the \code{KisSplice}
\code{counts} option used (see Details below).}
  
\item{pairedEnd}{a logical indicating if the data is paired-end (\code{FALSE},
default). If set to \code{TRUE}, the sum of the counts from the pair of reads
will be computed.
It can be used along with \code{counts} option. By default, it is
assumed that, in the \code{KisSplice} command line, two reads of the same pair 
has been inputed as following each other. If it is not the case, see option 
\code{order}.}
  
\item{order}{a numeric vector indicating the actual order of the corresponding 
paired reads in the columns of the \code{KisSplice} output such that they can 
be summed. This option goes along with \code{pairedEnd = TRUE}, if the read 
pairs are not in the expected order (see \code{pairedEnd} option).
It has as many elements as there are samples in total.
For more information on this parameter, see Details in 
\code{\link{kissplice2counts}}.}
       
\item{exonicReads}{a logical indicating if exonic/intronic read counts will be 
kept (\code{TRUE}, default) or discareded (\code{FALSE}). This option only 
works if \code{counts = 2}.}

\item{k2rg}{a logical indicating if the input file is a 
\code{KisSplice2RefGenome} (\code{TRUE}) output or a 
\code{KisSplice} (\code{FALSE}, default) output file.}

\item{keep}{a character vector listing the names of the events to be kept for
the statistical test (for \code{k2rg = TRUE}, analyses all of the events by 
default). The test will be more sensitive the selected events. Event(s) name(s) 
must be part of this list: deletion, insertion, indel, IR, ES, altA, altD, altAD, alt, 
- (for unclassified events). For more information on this parameter, see Details 
in \code{\link{kissplice2counts}}.}

\item{remove}{a character vector listing the names of the events to remove
for the statistical test (for \code{k2rg = TRUE}, does not remove any event by 
default). The test will be more sensitive for the non-selected events.
Event(s) name(s) must be part of this list: deletion, insertion, indel, IR, ES, 
altA, altD, altAD, alt, - (for unclassified events), MULTI.
This option can not be used along with the \code{keep} option, unless ES is
one of the events to be kept. In this case, the \code{remove} option will
work on specific ES events.
For more information on this parameter, see Details in 
\code{\link{kissplice2counts}}.}
                           
\item{pvalue}{a numerical value indicating the p-value threshold below which
the events will be kept in the final data frame.}

\item{filterLowCountsVariants}{a numerical value indicating the global variant 
count value (see Details below) below which events are filtered out in order 
to increase statistical power of the analysis. Both variant must have a read
coverage below this value in order to remove the event. This filter is done 
after the normalization and the overdispersion estimation.}

\item{flagLowCountsConditions}{a numerical value indicating the global 
condition count value (see Details below) below which we flag events as 
'lowCounts' in the final data frame. At least n-1 conditions (over n 
conditions) must have low counts to flag the event as 'lowCounts'.}

\item{technicalReplicates}{a boolean value indicating if the counts in 
\code{countsData} come from technical replicates only or not.}

\item{nbCore}{an integer indicating the number of cores to use for the model 
fitting step.}
                           
\item{adjPvalMax}{a double indicating the threshold for adjusted p-value. Only
SNVs/splicing events with an adjusted p-value lower than this threshold will be
kept in the output file.}

\item{dPSImin}{a double indicating the threshold for the deltaPSI. Only
SNVs/splicing events having an absolute value of deltaf/deltaPSI
higher than this threshold will be kept in the output file.}

\item{writePSI}{a boolean indicating if the user wants the f/PSI table to be
printed (\code{TRUE}, default) along with the final table (\code{FALSE}).}

\item{doQualityControl}{a boolean indicating if the user wants quality control plots to be written in the output folder (\code{TRUE}, default) or not (\code{FALSE}). See details in \code{\link{qualityControl}}.}

\item{resultsInShiny}{a boolean indicating if the user wants the results to be printed in a Shiny application (\code{TRUE}, default) or not (\code{FALSE}). See details in \code{\link{exploreResults}}.}
}


\value{
None.
}


\examples{
kissplice2refgenome_file <- system.file("extdata", 
    "output_k2rg_alt_splicing.txt", package="kissDE")
mySplicingconditions <- c("C1", "C1", "C2", "C2")
#kissDE(fileName=kissplice2refgenome_file, conditions=mySplicingconditions, 
#       output="results.tsv", counts=2, pairedEnd=TRUE, k2rg=TRUE)
}
