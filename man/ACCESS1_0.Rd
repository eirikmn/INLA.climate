\name{GISS_E2_R}
\alias{GISS_E2_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ACCESS1.0 temperature and forcing dataset.
}
\description{
This is the ACCESS1.0 temperature dataset. This Also includes radiative forcing data obtained by combining the ACCESS1.0 forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.0.
}
\usage{
data(ACCESS1_0)
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
data(ACCESS1_0, package = "INLA.climate")

#Plot temperature
plot(ACCESS1_0$Year, ACCESS1_0$Temperature)

#Plot forcing
plot(ACCESS1_0$Year, ACCESS1_0$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
