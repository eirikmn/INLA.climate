\name{CCSM4}
\alias{CCSM4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CCSM4 temperature and forcing dataset.
}
\description{
This is the CCSM4 temperature dataset. Also includes radiative forcing data obtained by combining the GISS-E2-R forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.6.
}
\usage{
data(CCSM4)
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
data(CCSM4, package = "INLA.climate")

#Plot temperature
plot(CCSM4$Year,GISS_E2_R$Temperature)

#Plot forcing
plot(CCSM4$Year,GISS_E2_R$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
