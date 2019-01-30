\name{print.summary.inla.climate}
\alias{print.summary.inla.climate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Print function for a \code{summary.inla.climate} object
}
\description{
  Displays the return object from \code{summary.inla.climate} to an easily readable format.
}
\usage{
## S3 method for class 'inla.climate'
\method{print}{summary.inla.climate}(x, digits=4L,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Object of type \code{summary.inla.climate}.
}
\item{digits}{
Number of digits to print.
}
\item{...}{
Other arguments.
}
}

\author{
  Eirik Myrvoll-Nilsen \href{eirik.myrvoll-nilsen@uit.no}{eirik.myrvoll-nilsen@uit.no}
}
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
