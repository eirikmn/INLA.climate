\name{GFDL_CM3}
\alias{GFDL_CM3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GFDL_CM3 temperature and forcing dataset.
}
\description{
This is the GFDL-CM3 temperature dataset. Also includes radiative forcing data obtained by combining the GISS-E2-R forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.0.
}
\usage{
data(GFDL_CM3)
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
Donner, Leo & Wyman, Bruce & Hemler, Richard & Horowitz, Larry & Ming, Yi & Zhao, Ming & Golaz, J-C & Ginoux, Paul & Lin, S-J & Schwarzkopf, M. & Austin, John & Alaka, Ghassan & Cooke, William & Delworth, Thomas & Freidenreich, S. & Gordon, C & Griffies, Stephen & Held, Isaac & Hurlin, William & Zeng, Fanrong. (2011). The Dynamical Core, Physical Parameterizations, and Basic Simulation Characteristics of the Atmospheric Component AM3 of the GFDL Global Coupled Model CM3. Journal of Climate. 24. 3484-3519. 10.1175/2011JCLI3955.1. 

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

}
\examples{
# Load data
data(GFDL_CM3, package = "INLA.climate")

#Plot temperature
plot(GFDL_CM3$Year,GFDL_CM3$Temperature)

#Plot forcing
plot(GFDL_CM3$Year,GFDL_CM3$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
