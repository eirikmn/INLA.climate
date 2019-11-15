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
