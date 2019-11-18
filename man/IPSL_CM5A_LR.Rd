\name{IPSL_CM5A_LR}
\alias{IPSL_CM5A_LR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
IPSL-CM5A-LR temperature and forcing dataset.
}
\description{
This is the IPSL-CM5A-LR temperature dataset. Also includes radiative forcing data obtained by combining the IPSL-CM5A-LR forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.1.
}
\usage{
data(IPSL_CM5A_LR)
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
Dufresne, Jean-Louis & Foujols, Marie-Alice & Denvil, S. & Caubel, Arnaud & Marti, Olivier & Aumont, Olivier & Balkanski, Yves & Bekki, Slimane & Bellenger, Hugo & Benshila, Rachid & Bony, Sandrine & Bopp, Laurent & Braconnot, Pascale & Brockmann, Patrick & Cadule, Patricia & Cheruy, F. & Codron, Francis & Cozic, Anne & Cugnet, David & Vuichard, N.. (2013). Climate change projections using the IPSL-CM5 earth system model: from CMIP3 to CMIP5. Climate Dynamics. 40. 2123–2165. 10.1007/s00382-012-1636-1. 

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

}
\examples{
# Load data
data(IPSL_CM5A_LR, package = "INLA.climate")

#Plot temperature
plot(IPSL_CM5A_LR$Year,IPSL_CM5A_LR$Temperature)

#Plot forcing
plot(IPSL_CM5A_LR$Year,IPSL_CM5A_LR$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
