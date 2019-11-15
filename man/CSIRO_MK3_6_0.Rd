\name{CSIRO_MK3_6_0}
\alias{CSIRO_MK3_6_0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CSIRO-MK3-6-0 temperature and forcing dataset.
}
\description{
This is the CSIRO-MK3-6-0 temperature dataset. Also includes radiative forcing data obtained by combining the CSIRO-MK3-6-0 forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 2.6.
}
\usage{
data(CSIRO_MK3_6_0)
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
data(CSIRO_MK3_6_0, package = "INLA.climate")

#Plot temperature
plot(CSIRO_MK3_6_0$Year,CSIRO_MK3_6_0$Temperature)

#Plot forcing
plot(CSIRO_MK3_6_0$Year,CSIRO_MK3_6_0$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
