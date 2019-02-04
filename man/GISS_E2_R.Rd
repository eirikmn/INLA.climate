\name{GISS_E2_R}
\alias{GISS_E2_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GISS-E2-R temperature and forcing dataset.
}
\description{
This is the GISS-E2-R temperature dataset, submitted by the Goddard Institute for Space Studies (GISS) for the CMIP5 data repository. Also includes radiative forcing data from ... .  The forcing data has been modified by taking the moving averages over a 17 year period in order to make it more similar the adjusted forcing of Forster (SITER NOE).
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

Hansen, J., R. Ruedy, M. Sato, and K. Lo, 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.
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
