\name{h.map.maker}
\alias{h.map.maker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hurst exponent to weighted AR(1) sum mapping
}
\description{
Creates a mapping between the Hurst exponent of a fractional Gaussian noise and the weights and lag-one parameters of the specified AR(1) mixture approximation. The mapping functions are drawn by splines based on precomputed values included in the package.
}
\usage{
h.map.maker(m=4, lagmax=1000, model="fgn")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
The number of AR(1) processes used in the approximate model.
}
\item{lagmax}{
The number of lags of which the approximate autocorrelation functions are supposed to match the theoretical fGn one.
}
\item{model}{
The long memory process for which \code{H} is a parameter of. Currently, only fractional Gaussian noise (\code{model="fgn"}) and ARFIMA(0,d,0) (\code{model="arfima"},d=H-0.5) are supported.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list of functions. The \code{m} first maps the Hurst exponent to the lag-one correlation parameters, and the last \code{m} will map to the weights.
}
\references{
%% ~put references to the literature/web site here ~
Sørbye, S., Myrvoll-Nilsen, E. and Rue, H. (2018)
  \emph{An approximate fractional Gaussian noise model with O(n) computational cost}.
  Statistics and Computing.
}
\author{
Eirik Myrvoll-Nilsen \href{eirik.myrvoll-nilsen@uit.no}{eirik.myrvoll-nilsen@uit.no}
}
\seealso{
\code{\link{inla.climate}}
}
\examples{
m <- 4
lagmax <- 1000
H <- 0.8

h.mapping <- h.map.maker(m = m, lagmax = lagmax)

params <- numeric(m)
weights <- numeric(m)

for(i in 1:m){
  params[i] <- h.mapping[[i]](H)
  weights[i] <- h.mapping[[m+i]](H)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }% use one of  RShowDoc("KEYWORDS")
\keyword{ models }% __ONLY ONE__ keyword per line
