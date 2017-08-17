\name{kissplice2counts}
\alias{kissplice2counts}

\title{Pre-treatment of \code{KisSplice} outputs}

\description{Function that converts \code{KisSplice} \code{.fa} outputs and
\code{KisSplice2RefGenome} outputs to data frame that can be used by others
functions of the \code{KissDE} package.
}


\usage{
kissplice2counts(fileName, counts = 0, pairedEnd = FALSE, order = NULL, 
exonicReads = TRUE, k2rg = FALSE, keep = c("All"), remove = NULL)
}


\arguments{
\item{fileName}{a string indicating the path to the \code{.fa} file that must
be pre-treated.}

\item{counts}{an interger (0, 1 or 2) corresponding with the \code{KisSplice}
\code{counts} option used. By default, as in \code{KisSplice}, it is turned
to 0, assuming there is no special counting option. If the data used has been
produced using \code{counts} 1 or 2 in \code{KisSplice}, this option sums the
counts on the junctions provided by \code{counts} to obtain counts by
condition. In this example, there is a part of a \code{KisSplice} output with
\code{counts} = 2:

|AS1_0|SB1_0|S1_0|ASSB1_0|AS2_27|SB2_41|S2_0|ASSB2_21| 

AS3_0|SB3_0|S3_0|ASSB3_0|AS4_7|SB4_8|S4_0|ASSB4_2

In a regular \code{KisSplice} output, it would be

|C1_0|C2_47|C3_1|C4_13 (with 47 = 27+41-0-21 and 13 = 7+8-0-2)}

\item{pairedEnd}{a logical indicating if the data is paired-end (not by
default), this function realises the sum of the counts from the two reads of
the pair. It can be used along with \code{counts} option. By default, it is
assumed that the data is "in order", which means two reads of the same pair
has been put next to each other in \code{KisSplice} output. If it is not the
case, see option \code{order}.}

\item{order}{a numeric vector indicating the actual order of the reads so
that they can be summed. This option goes with \code{pairedEnd} option. If
the data is not organised, the reads corresponding to a pair have not been
input one after another, we need to explicitely give the order. For instance,
if there are two pairs and if the input in \code{KisSplice} has been
-r sample1_reads1.fa -r sample2_reads1.fa -r sample1_reads2.fa 
-r sample2_reads2.fa, 
the input is not organised as two sets of reads coming from the same pair (1
or 2) which are not next to each other. The vector order to give would be
\code{order = c(1, 2, 1, 2)}.}

\item{exonicReads}{a logical indicating if exonic reads will be taken into
account in the counts. 
It is valid to use the \code{kissplice2counts} function only if the 
\code{counts} option is at 2 as the information about junction and exonic
reads is needed.}

\item{k2rg}{a logical indicating if you are using 
\code{KisSplice2RefGenome}'s output file instead of \code{KisSplice}'s file}

\item{keep}{a character vector indicating the names of the events to keep for
the statistical test. The test will be more sensitive for those events, the
other ones will not be tested.
Event(s) name(s) must be: deletion, insertion, IR, ES, altA, altD, altAD, 
alt, unclassified, unclassifiedSNP, -. If you only want to work with intron
retention events and altAD events, the vector should be 
\code{keep = c("IR","altAD")}.
If you want to work on specific ES events (i.e., ES_MULTI), the vector should
be \code{keep = c("ES")} and the \code{remove} option should be used too.
By default, all events are kept. This option goes only with \code{k2rg}
option.}

\item{remove}{a character vector indicating the names of the events to remove
for the statistical test. The test will be more sensitive for the other
events. Event(s) name(s) must be: deletion, insertion, IR, ES, altA, altD,
altAD, alt, unclassified, unclassifiedSNP, -, MULTI. If you want to work on
all events except insertions and deletions, the vector should be 
\code{remove = c("insertion","deletion")}.
This option can not be used if the \code{keep} option is used, unless ES is
one of the event(s) to be kept. In this case, the \code{remove} option will
work on specific ES events. If you want to work on ES_MULTI, then you should
do \code{keep=c("ES")} with \code{remove = c("alt","altA","altD","altAD")}.
By default, none events are removed.
This option goes only with \code{k2rg} option.}
}


\value{\code{kissplice2counts} returns a list of 4 objects:
\item{countsEvents}{a data frame containing several columns: a first column
(\code{events.names}) containing the name of the event based on 
\code{KisSplice} notation, a second one (\code{events.length}) containing the
length of the event, and several others (\code{counts1} to \code{countsN})
containing the counts corresponding to the replicates of the conditions.}

\item{psiInfo}{a data frame containing only one column (\code{events.names})
containing the name of the event based on \code{KisSplice} notation, as in
the first column of the \code{countsEvents} output.}

\item{exonicReadsInfo}{a logical indicating if you are using exonicReads.}

\item{k2rgFile}{a string containing the \code{KisSplice2RefGenome} file's
name. If you are working with a file coming from KisSplice and not 
\code{KisSplice2RefGenome}, this object is equal to NULL.}
}


\details{
The \code{kissplice2counts} function returns a list of 5 elements but only
the \code{countsEvents} one is showed when the output is printed.
}


\examples{
fpath <- system.file("extdata", "output_kissplice_SNV.fa", package="kissDE")
mySNVcounts <- kissplice2counts(fpath, pairedEnd=TRUE)
head(mySNVcounts)
names(mySNVcounts)
}