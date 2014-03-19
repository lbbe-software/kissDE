\name{diffExpressedEvents}
\alias{diffExpressedEvents}
\title{Retrieve condition-specific variants in RNA-seq data.}
\description{Retrieve condition-specific variants in RNA-seq data. 
}
\usage{diffExpressedEvents(n,nr,file)}
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
	\item{ID}{the variation ID}
	\item{Length}{the size of the first variant}
	\item{UP_Condn_Rm_Norm}{(resp LP_Condn_Rm_Norm) returns the normalised counts for the first variant (UP, resp. second variant : LP), for the condition n (Condn) and the replicate m (Rm)}
	\item{pALLGlobalPhi.glm}{adjusted pvalue}
	\item{Deltaf/DeltaPSI}{difference of relative abundance of variants across conditions}
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\author{
}


\references{
%%  ~~ possibly secondary sources and usages ~~
}
% \examples{

% n<-2 # two biological conditions
% m<-c(3,2) # three replicates for each condition
% data.file<-"my_input_file.fa"
% events<-diffExpressedEvents(n,m,data.file)

% }
\keyword{}
