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
  \item{counts=0}{ an option to set up if KisSplice --counts option was on. By default, as in KisSplice it is turned to 0, assuming there is no special counting option. If the data used has been produced using --counts 1 or 2, this option sums the counts on the junctions provided by --counts to obtain counts by condition. 
   In this example, there is a part of a KisSplice output with --counts =2 : |AS1_0|SB1_0|S1_0|ASSB1_0|AS2_27|SB2_41|S2_0|ASSB2_21|AS3_1|SB3_0|S3_0|ASSB3_0|AS4_7|SB4_8|S4_0|ASSB4_2
   In a regular KisSplice output it would be :
   |C1_0|C2_47|C3_0|C4_14 (with 47 = 41+27-21 and 14 = 1+7+8-2)}
  \item{pairedEnd=FALSE}{ if the data is paired-end (not by default), this function realises the sum of the counts from the two reads of the pair. It can be used along with counts option. By default, it is assumed that the data is "in order", which means two reads of the same pair has been put next to each other in KisSplice output. If it is not the case, see option order}
 \item{order=NULL}{ this option goes with pairedEnd option. If the data is not organised : the reads corresponding to a pair have not been input one after another, the actual order of the reads can be given in a vector so that they can be summed. For instance, if there are two pairs and if the input in KisSplice has been :
 -r first_set_pair1.fa -r first_set_pair2.fa -r second_set_pair1.fa -r second_set_pair2.fa
 the input is not organised as two sets of reads coming from the same pair (1 or 2) are not next to each other. The vector order to give would be :
 order=c("p1","p2","p1","p2").}
 \item{exonicReads=TRUE}{ if this option is turned to FALSE, exonic reads will no be taken into account in the counts. It is valid to use this function only if the counts option is at 1 or 2 as the information about junction and exonic reads is needed.}


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

