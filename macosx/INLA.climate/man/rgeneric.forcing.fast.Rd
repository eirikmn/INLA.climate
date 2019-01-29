\name{rgeneric.forcing.fast}
\alias{rgeneric.forcing.fast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rgeneric model for fractional Gaussian noise forcing model for INLA.
}
\description{
A script to allow fractional Gaussian noise adapted to the climate forcing model to be passed into the \code{R-INLA} framework.
}
\usage{
rgeneric.forcing.fast(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cmd}{
%%     ~~Describe \code{cmd} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
}
\details{
This function is only used to 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
When used as an input argument for \code{inla.rgeneric.define} this will return a model eligible to be used within the R-INLA framework.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

\dontrun{
  if(require("INLA",quietly=TRUE)){
  library(longmemo)
  n <- 300
  H <- 0.8  
  m <- 4

  forcing <- arima.sim(n=n, list(ar=c(0.6)),sd = 1)
  fgnsim <- simFGN0(n,H)
  
  fgnmean <- numeric(n)
  struktur <- (0.5+seq(0,n-1,length.out=n))^(H-3/2)    
  for(i in 1:n){
      fgnmean[i] <- rev(struktur[1:i])\%*\%forcing[1:i]
  }
  y <- fgnsim + fgnmean
  
  Hprior <- compute.Hprior(n, 0.9, 0.1)

  funks <- h.map.maker(m=m)
  model.approx <- inla.rgeneric.define(rgeneric.forcing.fast,lprior.fun.H = Hprior,
                                        n=n,N=m,forcing=forcing,funks=funks)
  formula = y ~ -1+ f(idy, model=model.approx)
  result.approx = inla(formula,family="gaussian", data=data.frame(y=y,idy=1:n),num.threads=1,
                       control.inla=list(reordering="metis"),
                       control.family = list(hyper = list(prec = list(initial = 12, fixed=TRUE))))

}
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
