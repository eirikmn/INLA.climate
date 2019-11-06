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
inla.climate(data, forcing, Qco2 = NULL,compute.mu=NULL, stepLength = 0.01, restart.inla = FALSE, m = 4, model = "fgn", formula=NULL, print.progress = FALSE,
   inla.options = list(),
   tcr.options = list(),
   mu.options = list(),
   ar1.options = list())
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
  \item{compute.mu}{
    Decides if and how the \code{forcing} temperature response should be explicitly computed. Doing so will significantly increase runtime. If method \code{compute.mu=1} is chosen, the samples will not be stored and the quantiles will not be computed. If \code{compute.mu=2} is chosen, full Bayesian inference is returned. Method \code{2} requires significantly more memory.
}
  \item{stepLength}{
    Numerical value which sets the steplength for the numerical scheme within INLA. A poor convergence can sometimes be improved by adjusting this.
  }
  \item{restart.inla}{
    Boolean variable indicating whether INLA should restart at the solution found from the first convergence of the Newton-Rhapson algorithm. This can sometimes improve convergence.
  }
  \item{m}{
     m The number of AR(1) processes used in the fGn approximation. I recommend having at least four (up to six, depending on \code{model} and \code{lagmax}). The model will perform slower, but more accurately if this is increased.
}
  \item{model}{
    The model used to describe stochastic forcing. For long memory responses the fractional Gaussian noise (\code{model="fgn"}) and ARFIMA(0,d,0) (\code{model="arfima"}) are supported. The first-order autoregressive model (\code{model="ar1"}) is also supported.
  }
  \item{formula}{
    The formula for the additive predictor. If \code{formula=NULL} the default option of forced and unforced temperature response, with no additional trend is selected.
  }
  \item{print.progress}{
    Prints progression if \code{TRUE}.
  }

  \item{inla.options}{
    A list of options for the \code{R-INLA} fitting procedure. See \link[INLA]{control.inla}, \link[INLA]{control.compute}, \link[INLA]{control.family} and \link[INLA]{control.predictor} for detailed description of options. There are, however, some default values that has been overridden in \code{INLA.climate}:
    \code{num.threads = 1},
    \code{control.compute\$dic = TRUE} if there are no \code{NA} values in \code{data},
    \code{control.compute\$config = TRUE},
    \code{control.inla\$reordering = "metis"},
    \code{control.family\$hyper\$prec\$initial = 12},
    \code{control.family\$hyper\$prec\$fixed = TRUE}
  }
  \item{tcr.options}{
  A list of options for the transient climate response sampling procedure. The following options can be set:
  \code{nsamples} denotes the number of Monte Carlo simulations used to estimate the TCR. Default value is 100000.
  \code{seed} is used as initial input for the random number generator. Default value is 1234.
  }
  \item{mu.options}{
  A list of options for the forcing response sampling procedure. The following options can be set:
  \code{full.Bayesian} is a boolean variable that indicates if full Bayesian inference is to be computed. If \code{TRUE}, all samples for each point in time is stored. This will increase runtime and memory consumption considerably. If \code{FALSE}, only the posterior mean and standard deviation for each point is computed. Default value is \code{FALSE}.
  \code{nsamples} denotes the number of Monte Carlo simulations used to estimate the forcing response. Default value is 100000.
  \code{seed} is used as initial input for the random number generator. Default value is 1234.
  }
  \item{ar1.options}{
  A list of options for the AR(1) weights and parameter evaluation. The following options can be set:
  \code{nsamples} denotes the number of Monte Carlo simulations used to estimate the forcing response. Default value is 100000.
  \code{seed} is used as initial input for the random number generator. Default value is 1234.
  }
}
\details{
If \code{length(data)>length(forcing)} or \code{data} contains \code{NA} values the missing values will be predicted using \link[INLA]{inla}. For return value \code{results}, a summary can be displayed using \code{summary.inla.climate(results)}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{\code{inla.result}}{An object of class "inla" returned by the \code{inla}-function. See \link[INLA]{inla} for more details.}
 \item{\code{hyperparam}}{A list containing the marginal posterior distribution, its mean, standard deviation and the 0.025, 0.5 and 0.975 quantiles of each hyperparameters.}
  \item{\code{latent.field}}{This list contains the marginal posterior mean, standard deviation and the 0.025, 0.5 and 0.975 quantiles of the latent Gaussian variables fitted to the data and the underlying AR(1) components. Also includes predictions if there are any.}
  \item{\code{time}}{Time spent for convergence. \code{time$fGn} gives the time spent fitting the temperature to an fGn, \code{time$TCR} gives the time spent doing the Monte Carlo sampling for the transient climate response and \code{time$total} gives total time elapsed time for the \code{inla.climate} function.}
  \item{\code{TCR}}{A list containing the mean, standard deviation and 0.025, 0.5 and 0.975 quantiles obtained from the transient climate response simulation procedure. \code{TCR$samples} gives a list of the Monte Carlo samples for the TCR and the hyperparameters.}

\item{\code{misc}}{
A list containing the call of the \code{inla.climate} function, as well as the input variables \code{data}, \code{forcing}, \code{C}, \code{m}, \code{model}, and the options set for the INLA procedure and, if they are computed, the TCR and forcing response sampling schemes. It also includes the value for the constant term $T_0$, defined by the mean of the 20 first values of \code{y}.
}
\item{\code{log.mlikelihood}}{
The log marginal-likelihood of the fit obtained by INLA.
}
\item{\code{dic}}{The Deviance Information Criterion (DIC) for the fit, obtained by INLA.}
}
\author{
Eirik Myrvoll-Nilsen \href{eirik.myrvoll-nilsen@uit.no}{eirik.myrvoll-nilsen@uit.no}
}
\seealso{
\code{\link[INLA]{inla}},\code{\link{inla.climate.mu}}
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
\keyword{ts}% use one of  RShowDoc("KEYWORDS")
\keyword{model}
