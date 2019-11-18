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
\references{
Bi, Daohua & Dix, Martin & Marsland, S. & O'Farrell, Siobhan & Rashid, Harun & Uotila, P. & Hirst, & Kowalczyk, Eva & Golebiewski, & Sullivan, Arnold & Yan, Yizhong & Hannah, & Franklin, Charmaine & Sun, Zhongmin & Vohralik, & Watterson, & Fiedler, R. & Collier, Mark & Puri, Kamal. (2012). The ACCESS coupled model: Description, control climate and evaluation. Australian Meteorological and Oceanographic Journal. 63. 41-64. 10.22499/2.6301.004. 

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

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
