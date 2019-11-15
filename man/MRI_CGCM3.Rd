\name{MRI_CGCM3}
\alias{MRI_CGCM3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MRI-CGCM3 temperature and forcing dataset.
}
\description{
This is the MRI-CGCM3 temperature dataset. Also includes radiative forcing data obtained by combining the MRI-CGCM3 forcing data from Forster et al. (2013) and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.2.
}
\usage{
data(MRI_CGCM3)
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
data(MRI_CGCM3, package = "INLA.climate")

#Plot temperature
plot(MRI_CGCM3$Year,MRI_CGCM3$Temperature)

#Plot forcing
plot(MRI_CGCM3$Year,MRI_CGCM3$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
