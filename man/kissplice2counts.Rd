\name{kissplice2counts}
\alias{kissplice2counts}

\title{Pre-treatment of \code{KisSplice} outputs}

\description{Function that converts \code{KisSplice} \code{.fa} outputs to data frame that can be used by others functions of the \code{KissDE} package.
}

\usage{
  kissplice2counts(fileName, counts = 0, pairedEnd = FALSE, order = NULL, 
  exonicReads = TRUE, discoSNP = FALSE, k2rg = FALSE)
}

\arguments{
  \item{fileName}{a string indicating the path to the \code{.fa} file that must be pre-treated.}
  
  \item{counts}{an interger (0, 1 or 2) corresponding with the \code{KisSplice} \code{counts} option used.
    By default, as in \code{KisSplice}, it is turned to 0, assuming there is no special counting option.
    If the data used has been produced using \code{counts} 1 or 2 in \code{KisSplice}, this option sums the counts on the junctions 
    provided by \code{counts} to obtain counts by condition. In this example, there is a part of a \code{KisSplice} output with \code{counts} = 2:
    
    |AS1_0|SB1_0|S1_0|ASSB1_0|AS2_27|SB2_41|S2_0|ASSB2_21| 
    
    AS3_1|SB3_0|S3_0|ASSB3_0|AS4_7|SB4_8|S4_0|ASSB4_2
    
    In a regular \code{KisSplice} output, it would be
    
    |C1_0|C2_47|C3_0|C4_14 (with 47 = 41+27-21 and 14 = 1+7+8-2)}
  
  \item{pairedEnd}{a logical indicating if the data is paired-end (not by default), this function realises the sum of the counts from the two reads of the pair. 
    It can be used along with \code{counts} option. 
    By default, it is assumed that the data is "in order", which means two reads of the same pair has been put next to each other in \code{KisSplice} output.
    If it is not the case, see option \code{order}.}
  
  \item{order}{a numeric vector indicating the actual order of the reads so that they can be summed.
    This option goes with \code{pairedEnd} option. 
    If the data is not organised, the reads corresponding to a pair have not been input one after another, we need to explicitely give the order. 
    For instance, if there are two pairs and if the input in \code{KisSplice} has been
    
    -r sample1_reads1.fa -r sample2_reads1.fa -r sample1_reads2.fa 
    
    -r sample2_reads2.fa, 
    
    the input is not organised as two sets of reads coming from the same pair (1 or 2) which are not next to each other. 
    The vector order to give would be \code{order = c(1, 2, 1, 2)}.}
  
  \item{exonicReads}{a logical indicating if exonic reads will be taken into account in the counts. 
    It is valid to use the \code{kissplice2counts} function only if the \code{counts} option is at 2 as the information about junction and exonic reads is needed.}
  
  \item{discoSNP}{a logical indicating if you are using discoSNP++'s output file instead of \code{KisSplice}'s file}
  
  \item{k2rg}{a logical indicating if you are using kissplice2refgenome's output file instead of \code{KisSplice}'s file}
}


\value{\code{kissplice2counts} returns a list of 3 objects:
  \item{countsEvents}{a data frame containing several columns: a first column (\code{events.names}) containing the name of the event based on \code{KisSplice} notation, a second one (\code{X1}) containing the length of the event, and several others (\code{X2} to \code{Xn}) containing the counts corresponding to the replicates of the conditions.}
  \item{discoInfo}{a logical indicating if you are using discoSNP++'s output file instead of \code{KisSplice}'s file}
  \item{psiInfo}{a data frame containing only one column (\code{events.names}) containing the name of the event based on \code{KisSplice} notation, as in the first column of the \code{countsEvents} output.}
}


\examples{
  fpath <- system.file("extdata", "output_kissplice_SNP.fa", package = "kissDE")
  mySNPcounts <- kissplice2counts(fpath, pairedEnd = TRUE)
  names(mySNPcounts)
}
