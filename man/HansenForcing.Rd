\name{HansenForcing}
\alias{HansenForcing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hansen forcing dataset.
}
\description{
This is the Hansen forcing dataset from Hansen et al. (2010).
}
\usage{
data(HansenForcing)
}
%- maybe also 'usage' for other objects documented here.
\format{The data is a list that contains these objects:
  \describe{
    \item{Year}{
      Time index denoting the year of the observations.
  }
    \item{Forcing}{
      Annual adjusted radiative forcing.
  }
  }
}
\references{
Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

}
\examples{
# Load data
data(HansenForcing, package = "INLA.climate")

#Plot temperature
plot(HansenForcing$Year,HansenForcing$Temperature)

#Plot forcing
plot(HansenForcing$Year,HansenForcing$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
