\name{GFDL_ESM2G}
\alias{GFDL_ESM2G}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GFDL-ESM2G temperature and forcing dataset.
}
\description{
This is the GFDL-ESM2G temperature dataset. Also includes radiative forcing data obtained by combining the GISS-E2-R forcing data from Forster et al. (2013) and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.1.
}
\usage{
data(GFDL_ESM2G)
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
Dunne, J. & John, Jasmin & Adcroft, Alistair & Griffies, Stephen & Hallberg, Robert & Shevliakova, Elena & Ronald, Stouffer & Cooke, William & Dunne, Krista & Harrison, Matthew & Krasting, John & Malyshev, Sergey & Milly, Paul & Phillips, P. & Sentman, Lori & Samuels, Bonita & Spelman, Michael & Winton, Michael & Wittenberg, Andrew & Zadeh, Niki. (2012). GFDL’s ESM2 Global Coupled Climate–Carbon Earth System Models. Part I: Physical Formulation and Baseline Simulation Characteristics. Journal of Climate. 25. 6646-6665. 10.1175/JCLI-D-11-00560.1. 

Dunne, J. & John, Jasmin & Shevliakova, Elena & Ronald, Stouffer & Krasting, John & Malyshev, Sergey & Milly, Paul & Sentman, Lori & Adcroft, Alistair & Cooke, William & Dunne, Krista & Griffies, Stephen & Hallberg, Robert & Harrison, Matthew & Levy, Hiram & Wittenberg, Andrew & Phillips, Peter & Zadeh, Niki. (2013). GFDL’s ESM2 Global Coupled Climate–Carbon Earth System Models. Part II: Carbon System Formulation and Baseline Simulation Characteristics*. Journal of Climate. 26. 2247-2267. 10.1175/JCLI-D-12-00150.1. 

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

}

\examples{
# Load data
data(GFDL_ESM2G, package = "INLA.climate")

#Plot temperature
plot(GFDL_ESM2G$Year,GFDL_ESM2G$Temperature)

#Plot forcing
plot(GFDL_ESM2G$Year,GFDL_ESM2G$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
