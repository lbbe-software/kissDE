\name{diffExpressedEvents}
\alias{diffExpressedEvents}
\title{Differentially expressed events in RNA-seq data}
\description{Retrieves differentially expressed events in RNA-seq data.
}
\usage{diffExpressedEvent(n,nr,file)}
\arguments{
	\item{n}{}
	\item{nr}{}
	\item{file}{}
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\value{
	code{diffExpressedEvents} returns a data.table with :
	\item{ID}{the event ID}
	\item{Length}{the size of the sequence}
	\item{UP_Condn_Rm_Norm}{(resp LP_Condn_Rm_Norm) returns the counts for the upper path (UP, resp. lower path : LP), for the condition n (Condn) and the replicate m (Rm)}
	\item{pALLGlobalPhi.glm}{}
	\item{Deltaf/DeltaPSI}{}
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\author{
  Lilia Brinza and Vincent Lacroix
  Maintainers: vincent.miele@univ-lyon1.fr, camille.marchet@inria.fr
}


\references{
%%  ~~ possibly secondary sources and usages ~~
}
% \examples{

% n<-2 # two biological conditions
% m<-c(3,2) # three replicates for each condition
% data.file<-"my_input_file.txt"
% events<-diffExpressedEvents(n,m,data.file)

% }
\keyword{}
