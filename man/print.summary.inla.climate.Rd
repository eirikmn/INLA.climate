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
  data(GISS_E2_R)
  result.climate <- inla.climate(data=GISS_E2_R$Temperature,forcing=GISS_E2_R$Forcing)
  print(summary(result.climate))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
