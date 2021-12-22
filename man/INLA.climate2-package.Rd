\name{INLA.climate2-package}
\alias{INLA.climate2-package}
\alias{INLA.climate2}
\docType{package}
\title{
\packageTitle{INLA.climate2}
}
\description{
\packageDescription{INLA.climate2}
}
\details{

%The DESCRIPTION file:
\packageDESCRIPTION{INLA.climate2}
\packageIndices{INLA.climate2}
%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{INLA.climate2}

Maintainer: \packageMaintainer{INLA.climate2}
}
\references{
Rue, H., Martino, S. and Chopin, N. (2009)
  \emph{Approximate Bayesian inference for latent Gaussian models using integrated nested Laplace approximations (with discussion)}.
  Journal of the Royal Statistical Society, Series B.
%~~ Literature or other references for background information ~~
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%~~ documentation directory ~~
\keyword{ package }
\seealso{
\code{\link[INLA]{inla}}
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{

if(require("INLA",quietly=TRUE)){
  data(GISS_E2_R)

  m <- 4
  Qco2 <- 3.8
  compute.mu <- "full"

  result.climate <- inla.climate(data = GISS_E2_R$Temperature, forcing = GISS_E2_R$Forcing,
                          Qco2 = Qco2, compute.mu = compute.mu, m = m )

  summary.inla.climate
  summary(result.climate)
  print(result.climate$time)
  print(result.climate$hyperparam$means)

  plot(result.climate)
}


}
