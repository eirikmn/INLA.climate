\name{plot.inla.climate}
\alias{plot.inla.climate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot-function for the \code{inla.climate} model
}
\description{
  Plots results obtained from the INLA.climate fit and the TCR sampling procedure.
}
\usage{
## S3 method for class 'inla.climate'
\method{plot}{inla.climate}(x,  plot.random=TRUE,
                                plot.hyperpars=TRUE,
                                plot.tcr.samples=TRUE,
                                plot.mu = TRUE,
                                postscript = FALSE,
                                pdf = FALSE,
                                prefix = "inla.climate.plots/figure-",
                                ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     \code{inla.climate} object returned by \code{inla.climate} function. Also works if the TCR sampling procedure has been skipped.
}
\item{plot.random}{
Boolean variable indicating if marginal posterior means and 95\% credible intervals for the latent GMRF corresponding to the AR(1) mixture.
}
\item{plot.hyperpars}{Boolean indicating if the marginal posterior distributions of the hyperparameters should be plotted.}
\item{plot.tcr.samples}{
Boolean variable indicating if a histogram for the TCR samples should be plotted (should there be any).
}
\item{plot.mu}{
Boolean variable indicating if the forcing response that can be interpreted as the mean of the fGn should be plotted. The observations are also included in the plot.
}
\item{postscript}{
Boolean variable indicating if postscript files should be produced instead.
}
\item{pdf}{
Boolean variable indicating if pdf files should be produced instead.
}
\item{prefix}{
The prefix for created files. Additional numbering is added.
}
\item{...}{
Additional arguments to \code{postscripts()}, \code{pdf()} or \code{dev.new()}.
}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Eirik Myrvoll-Nilsen \href{eirik.myrvoll-nilsen@uit.no}{eirik.myrvoll-nilsen@uit.no}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{inla.climate}}
}
\examples{
if(require("INLA",quietly=TRUE)){
  library(longmemo)
  n <- 300
  H <- 0.8
  m <- 4
  Qco2 <- 3.8
  forcing <- arima.sim(n=n, list(ar=c(0.6)),sd = 1)
  fgnsim <- simFGN0(n,H)

  fgnmean <- numeric(n)
  struktur <- (0.5+seq(0,n-1,length.out=n))^(0.8-3/2)
  for(i in 1:n){
      fgnmean[i] <- rev(struktur[1:i])\%*\%forcing[1:i]
  }
  y <- fgnsim + fgnmean

  result.climate <- inla.climate(data=y,forcing=forcing,m=m,Qco2=Qco2)
  plot(result.climate)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
