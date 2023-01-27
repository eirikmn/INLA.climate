\name{HadCRUT}
\alias{HadCRUT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HadCRUT4 temperature dataset.
}
\description{
This is the HadCRUT4 temperature dataset.
}
\usage{
data(HadCRUT)
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
  }
}
\references{
Jones, P.D., New, M., Parker, D.E., Martin, S. and Rigor, I.G., 1999: Surface air temperature and its variations over the last 150 years. Reviews of Geophysics 37, 173-199, doi:10.1029/1999RG900002

Jones, P.D., Osborn, T.J., Briffa, K.R., Folland, C.K., Horton, B., Alexander, L.V., Parker, D.E. and Rayner, N.A., 2001: Adjusting for sampling density in grid-box land and ocean surface temperature time series. J. Geophys. Res. 106, 3371-3380, doi:10.1029/2000JD900564

Jones, P.D., Lister, D.H., Osborn, T.J., Harpham, C., Salmon, M. and Morice, C.P., 2012: Hemispheric and large-scale land surface air temperature variations: an extensive revision and an update to 2010. Journal of Geophysical Research 117, D05127, doi:10.1029/2011JD017139.

Kennedy J.J., Rayner, N.A., Smith, R.O., Saunby, M. and Parker, D.E., 2011: Reassessing biases and other uncertainties in sea-surface temperature observations measured in situ since 1850 part 2: biases and homogenisation. Journal of Geophysical Research 116, D14104, doi:10.1029/2010JD015220

Osborn, T.J. and Jones, P.D., 2014: The CRUTEM4 land-surface air temperature data set: construction, previous versions and dissemination via Google Earth. Earth System Science Data 6, 61-68, doi:10.5194/essd-6-61-2014

Morice, C.P., Kennedy, J.J., Rayner, N.A. and Jones, P.D., 2012: Quantifying uncertainties in global and regional temperature change using an ensemble of observational estimates: the HadCRUT4 dataset. Journal of Geophysical Research, 117, D08101, doi:10.1029/2011JD017187
}
\examples{
# Load data
data(HadCRUT, package = "INLA.climate")

#Plot temperature
plot(HadCRUT$Year,HadCRUT$Temperature)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
\keyword{ts}
