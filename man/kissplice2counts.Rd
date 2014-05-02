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
  \item{counts=FALSE}{ an option to set up if KisSplice --counts option was on. This option sums the counts on the junctions provided by --counts to obtain counts by condition. It is also useful if the data is paired-end, as it sums the total counts for the pair of reads. 
   In this example, there is a part of a KisSplice output with --counts =2 : |AS1_0|SB1_0|S1_0|ASSB1_0|AS2_27|SB2_41|S2_0|ASSB2_21|AS3_1|SB3_0|S3_0|ASSB3_0|AS4_7|SB4_8|S4_0|ASSB4_2
   Reads of the junctions number 1 and 2 are paired-ends reads, so are numbers 3 and 4. So in a regular KisSplice output it would be :
   |C1_47|C2_14 (with 47 = 41+27-21 and 14 = 1+7+8-2)}
  \item{countsWanted=NULL}{ a vector to provide if option counts is active, to let know what counts are desired. In the latter example, c("R1","R1","R2","R2") would be provided, to precise that reads of the junctions number 1 and 2 come from the same replicate of the experiment and must be summed, same things for 3 and 4.}
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

