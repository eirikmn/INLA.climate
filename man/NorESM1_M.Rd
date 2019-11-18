\name{NorESM1_M}
\alias{NorESM1_M}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NorESM1-M temperature and forcing dataset.
}
\description{
This is the NorESM1-M temperature dataset. Also includes radiative forcing data obtained by combining the NorESM1-M forcing data and Hansen et al. (2010) such that the 18-yr moving averages are equal. The forcing slope coefficient when assuming a 1\% annual increase of CO2, which is used to estimate the TCR, is 3.8.
}
\usage{
data(NorESM1_M)
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
Bentsen, M. & Bethke, I. & Debernard, Jens & Iversen, Trond & Kirkevåg, Alf & Seland, Øyvind & Drange, Helge & Roelandt, Caroline & Seierstad, I. & Hoose, Corinna & Kristjánsson, J.. (2012). The Norwegian Earth System Model, NorESM1-M - Part 1: Description and basic evaluation. Geoscientific Model Development Discussions. 5. 2843-2931. 10.5194/gmdd-5-2843-2012. 

Iversen, T., Bentsen, M., Bethke, I., Debernard, J. B., Kirkevåg, A., Seland, Ø., Drange, H., Kristjansson, J. E., Medhaug, I., Sand, M., and Seierstad, I. A.: The Norwegian Earth System Model, NorESM1-M – Part 2: Climate response and scenario projections, Geosci. Model Dev., 6, 389–415, https://doi.org/10.5194/gmd-6-389-2013, 2013.

Hansen, J. and Ruedy, R. and Sato, M. and and Lo, K. 2010: Global surface temperature change, Rev. Geophys., 48, RG4004, doi:10.1029/2010RG000345.

Forster, P. M., T. Andrews, P. Good, J. M. Gregory, L. S. Jackson, and M. Zelinka, 2013: Evaluating adjusted forcing and modelspread for historical and future scenarios in the cmip5 generation of climate  models. Journal of Geophysical Research: Atmo-spheres, 118 (3), 1139–1150, doi:10.1002/jgrd.50174

}
\examples{
# Load data
data(NorESM1_M, package = "INLA.climate")

#Plot temperature
plot(NorESM1_M$Year,NorESM1_M$Temperature)

#Plot forcing
plot(NorESM1_M$Year,NorESM1_M$Forcing)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
