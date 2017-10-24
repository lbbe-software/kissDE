\name{kissplice2counts}
\alias{kissplice2counts}

\title{Pre-treatment of \code{KisSplice} outputs}

\description{Function that converts \code{KisSplice} (\code{.fa}) output and
\code{KisSplice2RefGenome} output to data frame that can be used by other
functions of the \code{KissDE} package.
}


\usage{
kissplice2counts(fileName, counts = 0, pairedEnd = FALSE, order = NULL, 
    exonicReads = TRUE, k2rg = FALSE, keep = c("All"), remove = NULL)
}


\arguments{
\item{fileName}{a string indicating the path to the \code{KisSplice} 
(\code{.fa}) or the \code{KisSplice2RefGenome} file that must be converted.}

\item{counts}{an interger (0, 1 or 2) corresponding with the \code{KisSplice}
\code{counts} option used (see Details below).}

\item{pairedEnd}{a logical indicating if the data is paired-end (not by
default). This function computes the sum of the counts from the pair of reads.
It can be used along with \code{counts} option. By default, it is
assumed that the data is "in order", which means two reads of the same pair
has been put next to each other in \code{KisSplice} output. If it is not the
case, see option \code{order}.}

\item{order}{a numeric vector indicating the actual order of the reads so
that they can be summed. This option goes with \code{pairedEnd} option.
For more information on this parameter, see Details below.}

\item{exonicReads}{a logical indicating if exonic reads will be taken into
account in the counts. This option only works if \code{counts} is set to 2.}

\item{k2rg}{a logical indicating if the input file is a 
\code{KisSplice2RefGenome} output file instead of a \code{KisSplice} file}

\item{keep}{a character vector indicating the names of the events to keep for
the statistical test. The test will be more sensitive for those events, the
other ones will not be tested. Event(s) name(s) must be: deletion, insertion,
IR, ES, altA, altD, altAD, alt, unclassified.
By default, all events are kept.
This option works only when \code{k2rg} option is set to \code{TRUE}.
For more information on this parameter, see Details below.}

\item{remove}{a character vector indicating the names of the events to remove
for the statistical test. The test will be more sensitive for the other events.
Event(s) name(s) must be: deletion, insertion, IR, ES, altA, altD, altAD, alt,
unclassified, MULTI.
This option can not be used if the \code{keep} option is used, unless ES is
one of the events to be kept. In this case, the \code{remove} option will
work on specific ES events.
By default, not any events are removed.
This option works only when \code{k2rg} option is set to \code{TRUE}.
For more information on this parameter, see Details below.}
}


\value{\code{kissplice2counts} returns a list of 4 objects:
\item{countsEvents}{a data frame containing several columns: a first column
(\code{events.names}) containing the name of the event based on 
\code{KisSplice} notation, a second one (\code{events.length}) containing the
length of the event, and several others (\code{counts1} to \code{countsN})
containing the counts corresponding to the replicates of the conditions.}

\item{psiInfo}{a data frame containing information to compute the PSI values.
This data frame is used only when \code{counts} is different from 0.}

\item{exonicReadsInfo}{a logical indicating if \code{exonicReads} are used.}

\item{k2rgFile}{a string containing the \code{KisSplice2RefGenome} file's
name. It is equal to NULL if the input file comes from KisSplice and not from
\code{KisSplice2RefGenome}.}
}


\details{
The \code{kissplice2counts} function returns a list of 4 elements 
(\code{countsEvents}, \code{psiInfo}, \code{exonicReadsInfo} and
\code{k2rgFile}) but only \code{countsEvents} is showed when the output is
printed.


The \code{counts} parameter:

By default, as in \code{KisSplice}, the \code{counts} option is set
to 0, assuming there is no special counting option. Below, an example of the
upper path counts format output by
\code{KisSplice} when \code{counts} is set to 2:

|AS1_0|SB1_0|S1_0|ASSB1_0|AS2_27|SB2_41|S2_0|ASSB2_21|

AS3_0|SB3_0|S3_0|ASSB3_0|AS4_7|SB4_8|S4_0|ASSB4_2.

In a regular \code{KisSplice} output (\code{counts = 0}), it would be:

|C1_0|C2_47|C3_1|C4_13 (with 47 = 27+41+0-21 and 13 = 7+8+0-2)


The \code{order} parameter:

If the reads corresponding to a pair have not been input one after another,
the order need to be explicitely given to the \code{kissplice2counts} function.
For instance, if there are two pairs and if the input in \code{KisSplice} has 
been -r sample1_reads1.fa -r sample2_reads1.fa -r sample1_reads2.fa 
-r sample2_reads2.fa, the input is not organised as two sets of reads coming 
from the same pair (1 or 2) which are not next to each other. The vector order
to give would be \code{order = c(1, 2, 1, 2)}.


The \code{keep} and \code{remove} parameters:

The options \code{keep} and \code{remove} allow the user to select the type of
alternative splicing events that should be analysed. For example, to work only
with intron retention events and alternative donor events, the vector should
be \code{keep = c("IR","altD")}. To work on all events except insertions and
deletions, the vector should be \code{remove = c("insertion","deletion")}.
To work specifically on single exon skipping (ES) events (only
one exon can be included or excluded), both \code{keep}
and \code{remove} options must be used. The \code{keep} option should
be set to \code{c("ES")} and the \code{remove} option should be set to
\code{c("alt","altA","altD","altAD","MULTI")}.
}


\examples{
fpath <- system.file("extdata", "output_kissplice_SNV.fa", package="kissDE")
mySNVcounts <- kissplice2counts(fpath, pairedEnd=TRUE)
names(mySNVcounts)
str(mySNVcounts)
head(mySNVcounts$countsEvents)
}