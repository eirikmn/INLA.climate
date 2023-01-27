\name{HansenForcingGHG}
\alias{HansenForcingGHG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hansen forcing dataset.
}
\description{
This is the greenhouse gas component of an updated version of the Hansen forcing dataset from Hansen et al. (2010).
}
\usage{
data(HansenForcingGHG)
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
data(HansenForcingGHG, package = "INLA.climate2")

#Plot forcing
plot(HansenForcingGHG$Year,HansenForcingGHG$Forcing)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
