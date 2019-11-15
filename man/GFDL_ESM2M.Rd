\name{GFDL_ESM2M}
\alias{GFDL_ESM2M}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GFDL-ESM2M temperature and forcing dataset.
}
\description{
This is the GFDL-ESM2M temperature dataset. Also includes radiative forcing data obtained by combining the GISS-E2-R forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.4.
}
\usage{
data(GFDL_ESM2M)
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
data(GFDL_ESM2M, package = "INLA.climate")

#Plot temperature
plot(GFDL_ESM2M$Year,GFDL_ESM2M$Temperature)

#Plot forcing
plot(GFDL_ESM2M$Year,GFDL_ESM2M$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
