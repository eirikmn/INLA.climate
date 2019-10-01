\name{rgeneric.ar1}
\alias{rgeneric.ar1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rgeneric model for first order autoregressive forcing model for INLA.
}
\description{
A script to allow fractional Gaussian noise adapted to the climate forcing model to be passed into the \code{R-INLA} framework.
}
\usage{
rgeneric.ar1(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cmd}{
Functions that specifies the latent Gaussian field and the prior distributions of the hyperparameters.
}
  \item{theta}{
The hyperparameters of the latent Gaussian model in logarithmic scaling.
}
}
\details{
This function sets up the radiative forcing long memory process as an \code{R-INLA} \code{rgeneric} model. It is used internally by \code{INLA.climate} and is not intended to be used separately by the user.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
When used as an input argument for \code{inla.rgeneric.define} this will return a model eligible to be used within the R-INLA framework (see example).
}
\author{
Eirik Myrvoll-Nilsen <\email{eirikmn91@gmail.com}>
}
\seealso{
\code{\link[INLA]{rgeneric}}
}
\examples{

\dontrun{
  if(require("INLA",quietly=TRUE)){
  data(GISS_E2_R)
  n <- dim(GISS_E2_R)[1]
  
  
  model.approx <- inla.rgeneric.define(rgeneric.ar1,
                                        n=n,forcing=GISS_E2_R$Forcing)
                                        
  formula = y ~ -1+ f(idy, model=model.approx)
  result.approx = inla(formula,family="gaussian", data=data.frame(y=GISS_E2_R$Temperature,idy=1:n),
                        num.threads=1, control.inla=list(reordering="metis"),
                        control.family = list(hyper = list(prec = list(initial = 12, fixed=TRUE))))
}
}


}
