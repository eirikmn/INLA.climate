\name{inla.climate.ar1}
\alias{inla.climate.ar1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AR(1) weights and parameter evaluation
}
\description{
Computes Bayesian inference about the weights and first-order autocorrelation parameters of the AR(1) components based on the posterior distribution obtained from INLA.
}
\usage{
inla.climate.ar1(result, m, nsamples = 100000, seed = 1234, print.progress = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
     An \code{inla.climate} or an \code{inla} object which can be accessed using \code{results$inla.result}, where \code{results} is the returned object from the \code{inla.climate} function.
}
  \item{m}{
     The number of AR(1) components used to model the temperature response.
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
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
If \code{result} is a list of class \code{inla.climate}, then this function returns the same object, but appends the following objects:
\item{\code{result$ar1}}{A list containing the mean, standard deviation, quantiles and samples for the weights and first lag correlation parameter. If \code{m} is equal to one, the weight is equal to one and the posterior distribution for the first-lag correlation parameter is given instead of the samples.}
\item{\code{result$time$ar1}}{The number of Monte Carlo simulations used. They are only used if \code{m>1}.}
\item{\code{result$misc$ar1.options$nsamples}}{The number of Monte Carlo simulations used. They are only used if \code{m>1}.}
\item{\code{result$misc$ar1.options$seed}}{The seed used for the random number generator.}

If \code{result} is a list of class \code{inla}, then a list containing the elements of \code{result$ar1} and \code{result$time$ar1} above is returned.

}
\references{
%% ~put references to the literature/web site here ~
Fredriksen, H.B., Rypdal, M. (2017)
  \emph{Long-Range Persistence in Global Surface Temperatures Explained by Linear Multibox Energy Balance Models}.
  Journal of Climate 30 (18), 7157–7168.
}
\author{
  Eirik Myrvoll-Nilsen \href{eirik.myrvoll-nilsen@uit.no}{eirik.myrvoll-nilsen@uit.no}
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{inla.climate}}, \code{\link{inla.climate.mu}}
}
\examples{
if(require("INLA",quietly=TRUE)){
  data(GISS_E2_R)
  result.climate <- inla.climate(data=GISS_E2_R$Temperature,forcing=GISS_E2_R$Forcing,m=3,model="ar1")
  result.tcr <- inla.climate.ar1(result.climate$inla.result,m=3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }% use one of  RShowDoc("KEYWORDS")
\keyword{ models }% __ONLY ONE__ keyword per line
