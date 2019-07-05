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
\method{plot}{inla.climate}(x,  plot.random=TRUE,
                                plot.ar=TRUE,
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
Boolean variable indicating if marginal posterior means and 95\% credible intervals for the latent GMRF corresponding to the AR(1) mixture should be plotted.
}
\item{plot.ar}{
Boolean variable indicating if marginal posterior means and 95\% credible intervals for the underlying AR(1) components should be plotted.
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

\author{
  Eirik Myrvoll-Nilsen \href{eirik.myrvoll-nilsen@uit.no}{eirik.myrvoll-nilsen@uit.no}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{inla.climate}}
}
\examples{
if(require("INLA",quietly=TRUE)){
  data(GISS_E2_R)
  result.climate <- inla.climate(data=GISS_E2_R$Temperature,forcing=GISS_E2_R$Forcing)
  plot(result.climate)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
