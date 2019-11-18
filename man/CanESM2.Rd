\name{CanESM2}
\alias{CanESM2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CanESM2 temperature and forcing dataset.
}
\description{
This is the CanESM2 temperature dataset. Also includes radiative forcing data obtained by combining the CanESM2 forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.8.
}
\usage{
data(CanESM2)
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
\references{
Chylek, Petr & Li, Jiangnan & Dubey, Manvendra & Wang, Muyin & Lesins, G.. (2011). Observed and model simulated 20th century Arctic temperature variability: Canadian Earth System Model CanESM2. Atmospheric Chemistry and Physics Discussions. 11. 22893-22907. 10.5194/acpd-11-22893-2011.

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

}
\examples{
# Load data
data(CanESM2, package = "INLA.climate")

#Plot temperature
plot(CanESM2$Year, CanESM2$Temperature)

#Plot forcing
plot(CanESM2$Year, CanESM2$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
