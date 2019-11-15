\name{NorESM1_M}
\alias{NorESM1_M}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NorESM1-M temperature and forcing dataset.
}
\description{
This is the NorESM1-M temperature dataset. Also includes radiative forcing data obtained by combining the NorESM1-M forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.8.
}
\usage{
data(NorESM1_M)
}
%- maybe also 'usage' for other objects documented here.
\format{The data is a list that contains these objects:
  \describe{
    \item{Year}{
      Time index denoting the year of the observations.
  }
  \item{Temperature}{
      Annual global mean surface temperature.
  }
    \item{Forcing}{
      Annual adjusted radiative forcing.
  }
  }
}
\examples{
# Load data
data(NorESM1_M, package = "INLA.climate")

#Plot temperature
plot(NorESM1_M$Year,NorESM1_M$Temperature)

#Plot forcing
plot(NorESM1_M$Year,NorESM1_M$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
