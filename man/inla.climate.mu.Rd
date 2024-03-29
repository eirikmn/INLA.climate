\name{inla.climate.mu}
\alias{inla.climate.mu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Temperature response estimation
}
\description{
Computes Bayesian inference about temperature response of \code{forcing} given parameter information contained in \code{result} obtained from Monte Carlo simulations generated by sampling from the posterior marginal distributions obtained by INLA.
}
\usage{
inla.climate.mu(result, forcing, quick=FALSE, T0.corr = NULL, nsamples = 100000, seed = 1234, print.progress = FALSE, model="fgn")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
     An \code{inla.climate} or an \code{inla} object which can be accessed using \code{results$inla.result}, where \code{results} is the returned object from the \code{inla.climate} function.
}
  \item{forcing}{
     The forcing for which the temperature response will be obtained.
}
  \item{quick}{
    Boolean variable. If \code{FALSE}, all Monte Carlo simulations will be stored. If \code{TRUE}, the quantiles will not be computed.
}
\item{T0.corr}{
  The temperature dataset analysed with INLA is shifted by \code{T0 = mean(data[1:20])} and the results generated by \code{inla.climate.mu} be shifted accordingly. If this variable is set equal to \code{T0}, the results will be aligned with the original input.
}
  \item{nsamples}{
     The number of Monte Carlo simulations used in obtaining transient climate response inference.
}
  \item{seed}{
     Seed used for random number generator.
}
  \item{print.progress}{
    Prints progression if \code{TRUE}.
}
  \item{model}{
    String which determines which model to be used when computing the temperature response. "fgn" and "arfima" produces a long-range dependent response, whilst "ar1" computes it as a sum of exponentials.
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list of class \code{inla.climate.mu} which contains the following objects:
\item{\code{mu.mean}}{The mean value of the temperature response due to forcing.}
\item{\code{mu.sd}}{The standard deviation of the temperature response due to forcing.}
\item{\code{mu.quant0.025}}{The 2.5\% quantile value of the temperature response due to forcing.}
\item{\code{mu.quant0.5}}{The 50\% quantile of the temperature response due to forcing.}
\item{\code{mu.quant0.975}}{The 97.5\% quantile of the temperature response due to forcing.}
\item{\code{samples}}{List containing the Monte Carlo simulations of the TCR and all hyperparameters.}
\item{\code{time}}{The time spent doing the Monte Carlo simulations.}
}
\author{
  Eirik Myrvoll-Nilsen \href{eirik.myrvoll-nilsen@uit.no}{eirik.myrvoll-nilsen@uit.no}
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{inla.climate}}
}
\examples{
if(require("INLA",quietly=TRUE)){
  data(GISS_E2_R)
  result.climate <- inla.climate(data=GISS_E2_R$Temperature,forcing=GISS_E2_R$Forcing)
  result.mu <- inla.climate.mu(result.climate,forcing = GISS_E2_R$Forcing,quick = FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }% use one of  RShowDoc("KEYWORDS")
\keyword{ model }% __ONLY ONE__ keyword per line
