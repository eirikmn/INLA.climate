\name{MIROC5}
\alias{MIROC5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MIROC5 temperature and forcing dataset.
}
\description{
This is the MIROC5 temperature dataset. Also includes radiative forcing data obtained by combining the MIROC5 forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 4.1.
}
\usage{
data(MIROC5)
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
Watanabe, M., T. Suzuki, R. O’ishi, Y. Komuro, S. Watanabe, S. Emori, T. Takemura, M. Chikira, T. Ogura, M. Sekiguchi, K. Takata, D. Yamazaki, T. Yokohata, T. Nozawa, H. Hasumi, H. Tatebe, and M. Kimoto, 2010: Improved Climate Simulation by MIROC5: Mean States, Variability, and Climate Sensitivity. J. Climate, 23, 6312–6335, https://doi.org/10.1175/2010JCLI3679.1

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

}
\examples{
# Load data
data(MIROC5, package = "INLA.climate")

#Plot temperature
plot(MIROC5$Year,MIROC5$Temperature)

#Plot forcing
plot(MIROC5$Year,MIROC5$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
