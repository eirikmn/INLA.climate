\name{GISS_E2_R}
\alias{GISS_E2_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GISS-E2-R temperature and forcing dataset.
}
\description{
This is the GISS-E2-R temperature dataset, submitted by the Goddard Institute for Space Studies (GISS) for the CMIP5 data repository. Also includes radiative forcing data obtained by combining the GISS-E2-R forcing data from Forster et al. (2013) and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.8.
}
\usage{
data(GISS_E2_R)
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
\source{
National Aeronautics and Space Administration (NASA) - Goddard Institute for Space Studies (GISS) \href{data.giss.nasa.gov/gistemp}{https://data.giss.nasa.gov/gistemp/}.
}
\references{
GISTEMP Team, 2018: GISS Surface Temperature Analysis (GISTEMP). NASA Goddard Institute for Space Studies. Dataset accessed 20YY-MM-DD at https://data.giss.nasa.gov/gistemp/.

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

}
\examples{
# Load data
data(GISS_E2_R, package = "INLA.climate")

#Plot temperature
plot(GISS_E2_R$Year,GISS_E2_R$Temperature)

#Plot forcing
plot(GISS_E2_R$Year,GISS_E2_R$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
