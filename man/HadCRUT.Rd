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
