\name{MPI_ESM_LR}
\alias{MPI_ESM_LR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MPI-ESM-LR temperature and forcing dataset.
}
\description{
This is the MPI-ESM-LR temperature dataset. Also includes radiative forcing data obtained by combining the MPI-ESM-LR forcing and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 4.1.
}
\usage{
data(MPI_ESM_LR)
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
Giorgetta, M. A., et al. ( 2013), Climate and carbon cycle changes from 1850 to 2100 in MPI‐ESM simulations for the Coupled Model Intercomparison Project phase 5, J. Adv. Model. Earth Syst., 5, 572– 597, doi:10.1002/jame.20038.

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

}

\examples{
# Load data
data(MPI_ESM_LR, package = "INLA.climate")

#Plot temperature
plot(MPI_ESM_LR$Year,MPI_ESM_LR$Temperature)

#Plot forcing
plot(MPI_ESM_LR$Year,MPI_ESM_LR$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
