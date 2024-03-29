\name{MIROC_ESM}
\alias{MIROC_ESM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MIROC_ESM temperature and forcing dataset.
}
\description{
This is the MIROC-ESM temperature dataset. Also includes radiative forcing data obtained by combining the MIROC-ESM forcing and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 4.3.
}
\usage{
data(MIROC_ESM)
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
Watanabe, Shingo & Hajima, Tomohiro & Sudo, Kengo & Nagashima, Tatsuya & Takemura, Toshihiko & Okajima, H. & Nozawa, Toru & Kawase, H. & Abe, Manabu & Yokohata, Tokuta & Ise, Toshifumi & Sato, Hisashi. (2011). MIROC-ESM: model description and basic results of CMIP5-20c3m experiments. Geoscientific Model Development Discussions. 4. 10.5194/gmdd-4-1063-2011.

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

}
\examples{
# Load data
data(MIROC_ESM, package = "INLA.climate")

#Plot temperature
plot(MIROC_ESM$Year,MIROC_ESM$Temperature)

#Plot forcing
plot(MIROC_ESM$Year,MIROC_ESM$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
