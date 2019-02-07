\name{print.inla.climate}
\alias{print.inla.climate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Print-function for the \code{inla.climate} model
}
\description{
  Prints key information regarding an INLA.climate fit and the TCR sampling procedure.
}
\usage{
\method{print}{inla.climate}(x, digits=4L,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     \code{inla.climate} object returned by \code{inla.climate} function. Also works if the TCR sampling procedure has been skipped.
}
\item{digits}{
Number of digits to print.
}
\item{...}{
Other arguments
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
  print(result.climate)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
