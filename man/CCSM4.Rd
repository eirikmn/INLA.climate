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
\references{
Gent, P.R., G. Danabasoglu, L.J. Donner, M.M. Holland, E.C. Hunke, S.R. Jayne, D.M. Lawrence, R.B. Neale, P.J. Rasch, M. Vertenstein, P.H. Worley, Z. Yang, and M. Zhang, 2011: The Community Climate System Model Version 4. J. Climate, 24, 4973–4991, https://doi.org/10.1175/2011JCLI4083.1

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

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
