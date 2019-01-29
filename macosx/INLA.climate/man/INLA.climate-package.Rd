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

%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
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
  summary(result.climate)
  plot(result.climate)

}


}
