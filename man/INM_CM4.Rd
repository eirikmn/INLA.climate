\name{INM_CM4}
\alias{INM_CM4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
INM-CM4 temperature and forcing dataset.
}
\description{
This is the INM-CM4 temperature dataset. Also includes radiative forcing data obtained by combining the INM-CM4 forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.0.
}
\usage{
data(INM_CM4)
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
data(INM_CM4, package = "INLA.climate")

#Plot temperature
plot(INM_CM4$Year,INM_CM4$Temperature)

#Plot forcing
plot(INM_CM4$Year,INM_CM4$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
