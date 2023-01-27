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
\references{
Yukimoto, S., Adachi, Y., Hosaka, M., Sakami, T., Yoshimura, H., Hirabara, M., Tanaka, T. Y., Shindo, E., Tsujino, H., Deushi, M., Mizuta,465R., Yabu, S., Obata, A., Nakano, H., Koshiro, T., Ose, T., and Kitoh, A. 2012: A New Global Climate Model of the Meteorological ResearchInstitute: MRI-CGCM3 –Model Description and Basic Performance–, J. Meteorolog. Soc. Jpn. Ser. II, 90A, 23–64.

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

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
