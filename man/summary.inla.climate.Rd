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

\value{
Returns an object of class \code{summary.inla.climate}, a list with components:
\item{call}{
The call for \code{object}.
}
\item{cpu.used}{
Time spent on various tasks for \code{object}.
}
\item{hyperpar}{
Summary statistics on the hyperparameters of \code{object}.
}\item{TCR}{
Summary statistics on the transient climate response of \code{object}, should there be any.
}\item{tcr.samples}{
Number of samples used to estimate the transient climate response of \code{object}, if such a procedure has been performed.
}
\item{mu.full.Bayesian}{
Boolean describing if a complete sampling procedure for the forcing response was computed for \code{object}.
}\item{mu.samples}{
Number of samples used to estimate the forcing response of \code{object}, if such a procedure has been performed.
}\item{neffp}{
Expected effective number of parameters in \code{object}.
}\item{mlik}{
The log marginal likelihood of \code{object}.
}\item{linear.predictor}{
The linear predictor of \code{object}, if it has been computed.
}\item{digits}{
The number of significant digits used for \code{summary.inla.climate}.
}\item{random.names}{
The names of the random effects of \code{object}.
}\item{family}{
The family assigned to the likelihood of \code{object}.
}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
  summary(result.climate)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
