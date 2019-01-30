\name{summary.inla.climate}
\alias{summary.inla.climate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summary-function for the \code{inla.climate} model
}
\description{
  Displays a clear summary of the information contained in the list returned by the \code{inla.climate} function.
}
\usage{
## S3 method for class 'inla.climate'
\method{summary}{inla.climate}(object, digits=4L,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
     List returned by \code{inla.climate} function. Also works if the TCR sampling procedure has been skipped.
}
\item{digits}{
Number of digits to print.
}
\item{...}{
Other arguments.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Eirik Myrvoll-Nilsen \href{eirik.myrvoll-nilsen@uit.no}{eirik.myrvoll-nilsen@uit.no}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{inla.climate}}
}
\examples{
if(require("INLA",quietly=TRUE)){
  library(longmemo)
  n <- 300
  H <- 0.8
  m <- 4
  Qco2 <- 3.8
  forcing <- arima.sim(n=n, list(ar=c(0.6)),sd = 1)
  fgnsim <- simFGN0(n,H)

  fgnmean <- numeric(n)
  struktur <- (0.5+seq(0,n-1,length.out=n))^(0.8-3/2)
  for(i in 1:n){
      fgnmean[i] <- rev(struktur[1:i])\%*\%forcing[1:i]
  }
  y <- fgnsim + fgnmean

  result.climate <- inla.climate(data=y,forcing=forcing,m=m,Qco2=Qco2)
  summary(result.climate)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
