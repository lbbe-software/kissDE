\name{kissplice2counts}
\alias{kissplice2counts}

\title{Pre-treatment of \code{KisSplice} outputs}

\description{Converts \code{KisSplice} \code{.fa} outputs to data.frame objects that can be used by others functions of the package.
}

\usage{
  kissplice2counts(fileName, counts = 0, pairedEnd = FALSE, order = NULL, 
  exonicReads = TRUE, discoSNP = FALSE, k2rg = FALSE)
}

\arguments{
  \item{fileName}{a character vector that gives the path to \code{.fa} file that must be pre-treated.}
  
  \item{counts}{an option to fill in if \code{KisSplice} \code{counts} option was on. 
    By default, as in \code{KisSplice} it is turned to 0, assuming there is no special counting option. 
    If the data used has been produced using \code{counts} 1 or 2, this option sums the counts on the junctions 
    provided by \code{counts} to obtain counts by condition. In this example, there is a part of a \code{KisSplice} output with \code{counts} = 2:
    
    |AS1_0|SB1_0|S1_0|ASSB1_0|AS2_27|SB2_41|S2_0|ASSB2_21| 
    
    AS3_1|SB3_0|S3_0|ASSB3_0|AS4_7|SB4_8|S4_0|ASSB4_2
    
    In a regular \code{KisSplice} output it would be
    
    |C1_0|C2_47|C3_0|C4_14 (with 47 = 41+27-21 and 14 = 1+7+8-2)}
  
  \item{pairedEnd}{if the data is paired-end (not by default), this function realises the sum of the counts from the two reads of the pair. 
    It can be used along with counts option. By default, it is assumed that the data is "in order", which means two reads of the same pair 
    has been put next to each other in \code{KisSplice} output. If it is not the case, see option \code{order}.}
  
  \item{order}{this option goes with \code{pairedEnd} option. If the data is not organised, the reads corresponding to a pair have not 
    been input one after another, we need to explicitely give the order. 
    The actual order of the reads is given in a numeric vector so that they can be summed. 
    For instance, if there are two pairs and if the input in \code{KisSplice} has been
    
    -r sample1_reads1.fa -r sample2_reads1.fa -r sample1_reads2.fa 
    
    -r sample2_reads2.fa, 
    
    the input is not organised as two sets of reads coming from the same pair (1 or 2) are not next to each other. 
    The vector order to give would be \code{order = c(1, 2, 1, 2)}.}
  
  \item{exonicReads}{if this option is turned to \code{FALSE}, exonic reads will no be taken into account in the counts. 
    It is valid to use this function only if the  \code{counts} option is at 2 as the information about junction and exonic reads is needed.}
  
  \item{discoSNP}{put to \code{TRUE} if you are using discoSNP++'s output file instead of \code{KisSplice}'s file}
  
  \item{k2rg}{put to \code{TRUE} if you are using kissplice2refgenome's output file instead of \code{KisSplice}'s file}
}

\value{\code{kissplice2counts} returns a data.frame object that contains several columns:
  \item{events.names}{the name of the event based on \code{KisSplice} notation}
  
  \item{V1}{contains the length of the event}
  
  \item{V2 ... to Vn}{contain the counts corresponding to the replicates of the conditions}
}

% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
% }
