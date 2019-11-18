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
\references{
Rotstayn, Leon & Jeffrey, Stephen & Collier, Mark & Dravitzki, S. & Hirst, A. & Syktus, Jozef & Wong, K.. (2012). Aerosol- and greenhouse gas-induced changes in summer rainfall and circulation in the Australasian region: A study using single-forcing climate simulations. Atmos. Chem. Phys.. 12. 6377-6404. 10.5194/acp-12-6377-2012. 

Jeffrey, Stephen & Rotstayn, Leon & Collier, Mark & Dravitzki, S & Hamalainen, Carlo & Moeseneder, Christian & Wong, K & Syktus, Jozef. (2013). Australia's CMIP5 submission using the CSIRO-Mk3.6 model. Australian Meteorological and Oceanographic Journal. 63. 1-13. 10.22499/2.6301.001. 

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

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
