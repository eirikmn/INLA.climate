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
