\name{compute.Hprior}
\alias{compute.Hprior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Penalise complexity prior for the Hurst exponent
}
\description{
Computes a penalise complexity prior function for the Hurst exponent in logarithmic scaling. Both fGn and ARFIMA(0,d,0) models are supported.
}
\usage{
compute.Hprior(n, upper, alpha, persistent=TRUE, model="fgn")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Length of the time series for which the prior is applied. Increasing this will theoretically improve the accuracy of the prior, but will be more costly to perform. Experiments have suggested that the prior is rather invariant of this parameter.
}
  \item{upper}{
The value of the penalise complexity parameter \code{upper}.
}
  \item{alpha}{
The value of the penalise complexity parameter \code{alpha}. This also denotes the probability that the Hurst exponent exceeds \code{upper}.
}
  \item{persistent}{
A boolean variable. If \code{TRUE}, the PC-prior will be created for only the persistent range \code{0.5<H<1}.
}
  \item{model}{
Determines which model the Hurst exponent belongs to. Currently, only fractional Gaussian noise (\code{model="fgn"}) and ARFIMA(0,d,0) (\code{model="arfima"}, d=H-0.5) models are supported.
}
}
\details{

The PC-prior parameters shape the density as the probability of the Hurst exponent exceeding \code{upper} is equal to \code{alpha}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns the PC-prior in logarithmic scaling as a function.
}
\references{
%% ~put references to the literature/web site here ~
Sørbye, S., Rue, H. (2017)
  \emph{Fractional Gaussian noise: Prior specification and model comparison}.
  Environmetrics, 29(5-6)
}
\author{
  Eirik Myrvoll-Nilsen \href{eirik.myrvoll-nilsen@uit.no}{eirik.myrvoll-nilsen@uit.no}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{inla.climate}}
}
\examples{
HH = seq(0.501,0.999,by=0.0001)
l.HH = log((HH-0.5)/(1-HH))

l.pcp.f = compute.Hprior(100,0.9,0.01,persistent=TRUE,model="fgn")
pcp = 0.5+0.5/(1+exp(-l.pcp.f(l.HH)))

plot(x=HH,y=pcp,type="l",xlab="H",ylab="Density")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
