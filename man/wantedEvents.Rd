\name{wantedEvents}
\alias{wantedEvents}

\title{Give the events to be used in the \code{KissDE} analysis with k2rg.}

\description{Function that converts \code{keep} and \code{remove} vectors into 
one vector containing the events to work with.
}

\usage{
  wantedEvents(keep = c("All"), remove = NULL)
}

\arguments{
  \item{keep}{a character vector indicating the names of the events to keep for
  the statistical test. The test will be more sensitive for those events, the
  other ones will not be tested. Event(s) name(s) must be: deletion, 
  insertion, IR, ES, altA, altD, altAD, alt, unclassified, unclassifiedSNP, -.
  By default, all events are kept.}
    
  \item{remove}{a character vector indicating the names of the events to remove
  for the statistical test. The test will be more sensitive for the other 
  events. Event(s) name(s) must be: deletion, insertion, IR, ES, altA, altD, 
  altAD, alt, unclassified, unclassifiedSNP, -, MULTI.
  This option can not be used if the \code{keep} option is used, unless ES is 
  one of the event(s) to be kept. In this case, the \code{remove} option will 
  work on specific ES events. By default, none events are removed.}
}

\value{
The vector containing the names of the events to work with.
}

\examples{
wantedEvents(keep=c("IR", "altAD"))
wantedEvents(remove=c("insertion", "deletion"))
wantedEvents(keep=c("altA", "ES"), remove=c("alt", "altA", "altD", "altAD"))
}