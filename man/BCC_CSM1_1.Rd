\name{BCC_CSM1_1}
\alias{BCC_CSM1_1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BCC-CSM1-1 temperature and forcing dataset.
}
\description{
This is the BCC-CSM1-1 temperature dataset. Also includes radiative forcing data obtained by combining the BCC-CSM-1 forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.2.
}
\usage{
data(BCC_CSM1_1)
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
data(BCC_CSM1_1, package = "INLA.climate")

#Plot temperature
plot(BCC_CSM1_1$Year, BCC_CSM1_1$Temperature)

#Plot forcing
plot(BCC_CSM1_1$Year, BCC_CSM1_1$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
