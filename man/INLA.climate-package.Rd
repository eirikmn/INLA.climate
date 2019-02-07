\name{INLA.climate-package}
\alias{INLA.climate-package}
\alias{INLA.climate}
\docType{package}
\title{
\packageTitle{INLA.climate}
}
\description{
\packageDescription{INLA.climate}
}
\details{

%The DESCRIPTION file:
\packageDESCRIPTION{INLA.climate}
\packageIndices{INLA.climate}
%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{INLA.climate}

Maintainer: \packageMaintainer{INLA.climate}
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
  print(result.climate$time)
  print(result.climate$hyperparam$means)
  summary(result.climate)
  plot(result.climate)
}


}
