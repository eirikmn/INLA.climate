\name{BCC_CSM1_1m}
\alias{BCC_CSM1_1m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BCC-CSM1-1(m) temperature and forcing dataset.
}
\description{
This is the BCC-CSM1-1(m) temperature dataset. Also includes radiative forcing data obtained by combining the BCC-CSM1-1(m) forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.6.
}
\usage{
data(BCC_CSM1_1m)
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
data(BCC_CSM1_1m, package = "INLA.climate")

#Plot temperature
plot(BCC_CSM1_1m$Year,BCC_CSM1_1m$Temperature)

#Plot forcing
plot(BCC_CSM1_1m$Year,BCC_CSM1_1m$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
