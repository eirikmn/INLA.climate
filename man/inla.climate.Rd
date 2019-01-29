\name{inla.climate}
\alias{inla.climate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An INLA climate forcing model for global mean surface temperature
}
\description{
Fits forcing data to a given temperature dataset using INLA. The undocumented forcing variables not included in \code{forcing} are collectively assumed a stochastic fGn process. Also uses Monte Carlo simulations to obtain Bayesian inference about the transient climate response.
}
\usage{
inla.climate(data, forcing, Qco2 = NULL,m = 4, stoc = "fgn", print.progress = FALSE,
   inla.options = list(compute.dic=TRUE, nruns=1, initialtheta=NULL, stepLength = 0.01, verbose = FALSE),
   tcr.options = list(mcmcsamples = 100000, seed = 1234))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     Global mean surface temperature. Must be a numeric vector (for now).
}
  \item{forcing}{
     The documented forcing data as a sum of all contributive forcing variables. Must be a numeric vector.
}
  \item{Qco2}{
     Qco2 The slope of the forcing when assuming a CO2 doubling in the future. Used for computing the transient climate response. Qco2 = NULL will skip the transient climate response sampling procedure.
}
  \item{m}{
     m The number of AR(1) processes used in the fGn approximation. I recommend having at least four (up to six, depending on \code{stoc} and \code{lagmax}). The model will perform slower, but more accurately if this is increased.
}
  \item{stoc}{
    The long memory model used to describe stochastic forcing. Currently, only fractional Gaussian noise (\code{stoc="fgn"}) and ARFIMA(0,d,0) (\code{stoc="arfima"}) are supported.
  }
  \item{print.progress}{
    Prints progression if \code{TRUE}.
  }

  \item{inla.options}{
    A list of options for the INLA fitting procedure.
    \code{compute.dic}
    If \code{TRUE} the deviance information criterion (DIC) is computed.
    \code{nruns}
    The number of times the INLA fitting procedure should be executed. Values obtained from convergence is used as initial values in the subsequent procedure. Increasing this is likely to improve converge.
    \code{initialtheta}
    The initial values used for the hyperparameters (in internal scaling) in the Newton-Raphson numerical scheme. Accuracy can sometimes be improved by repeating the algorithm with starting values set to the values obtained from the previous execution. For returned object \code{result}, these can be accessed by     \code{result$inla.result$mode$theta}.
    \code{inla.options$stepLength} sets the steplength for the numerical scheme. A poor convergence can sometimes be improved by adjusting this.
    \code{inla.options$verbose} prints the details of the INLA numerical procedure if left \code{TRUE}.
    \code{inla.options$Hprior} the prior used for the Hurst exponent. If left \code{NULL} a penalise complexity prior will be generated (somewhat time consuming).
    \code{inla.options$Hpriorparam} denotes the parameters to be used if the PC prior for the Hurst exponent should be generated.
  }
  \item{tcr.options}{
  A list of options for the transient climate response sampling procedure.
  \code{tcr.options$mcmcsamples} denotes the number of generated Monte Carlo simulations for the TCR.
  \code{seed} is used as initial input for the random number generator.
  }
}
\details{
If \code{length(data)>length(forcing)} or \code{data} contains \code{NA} values the missing values will be predicted using INLA. For return value \code{results}, a summary can be displayed using \code{summary.inla.climate(results)}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{\code{inla.result}}{An object of class "inla" returned by the \code{inla}-function. See \code{help(inla)} for more details.}
 \item{\code{hyperparam}}{A list containing the marginal posterior distribution, its mean, standard deviation and the 0.025, 0.5 and 0.975 quantiles of each hyperparameters.}
  \item{\code{latent}}{Contains the marginal posterior mean, standard deviation and the 0.025, 0.5 and 0.975 quantiles of the latent Gaussian variables fitted to the data. Also includes predictions if there are any.}
  \item{\code{time}}{Time spent for convergence. \code{time$fGn} gives the time spent fitting the temperature to an fGn, \code{time$TCR} gives the time spent doing the Monte Carlo sampling for the transient climate response and \code{time$total} gives total time elapsed time for the \code{inla.climate} function.}
  \item{\code{TCR}}{A list containing the mean, standard deviation and 0.025, 0.5 and 0.975 quantiles obtained from the transient climate response simulation procedure. \code{TCR$samples} gives a list of the Monte Carlo samples for the TCR and the hyperparameters.}

\item{\code{misc}}{
A list containing the call of the \code{inla.climate} function, as well as the input variables \code{C}, \code{m}, \code{stoc}, \code{mcmcsamples}, \code{seed}. Also includes the value for the constant term $T_0$, defined by the mean of the 20 first values of \code{y}.
}
\item{\code{log.mlikelihood}}{
The log marginal-likelihood of the fit obtained by INLA.
}
\item{\code{dic}}{The Deviance Information Criterion (DIC) for the fit, obtained by INLA.}
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
\code{\link{inla.climate.tcr}},\code{\link{inla.climate.mu}}
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
  summary.inla.climate(result.climate)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}% use one of  RShowDoc("KEYWORDS")
