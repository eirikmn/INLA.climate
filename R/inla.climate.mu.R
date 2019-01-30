inla.climate.mu = function(result,forcing,nsamples=100000,full.Bayesian=FALSE,seed=1234,print.progress=FALSE){

  atch = tryCatch(attachNamespace("INLA"),error=function(x){})
  if(length(find.package("INLA",quiet=TRUE))==0){
    stop("This function requires INLA. Please install at www.R-INLA.org or by calling 'install.packages(\"INLA\", repos=c(getOption(\"repos\"), INLA=\"https://inla.r-inla-download.org/R/testing\"), dep=TRUE)' from R.")
  }
  
  if(print.progress){
    cat("Starting Monte Carlo sampling with n=",nsamples," simulations..\n",sep="")
  }
  set.seed(seed)
  inla.seed = as.integer(runif(1)*.Machine$integer.max)

  #n=result$misc$configs$contents$length[1]
  n=length(forcing)

  #x = inla.posterior.sample(nsamples,r,seed=inla.seed) #int.strategy=grid
  x = INLA::inla.hyperpar.sample(nsamples,result)
  hyperpars = matrix(NA,nrow=nsamples,ncol=4) #c(H,sf,F0,TCR)


  hyperpars[,1] = 0.5+0.5/(1+exp(-x[,2]))
  hyperpars[,2] = 1/sqrt(exp(x[,3]))
  hyperpars[,3] = x[,4]

   tid.start = proc.time()[[3]]
  #if(!is.loaded('Rc_mu')){
  #print('hallo')
  #  dyn.load(file.path(.Library,"INLA.climate/libs/Rc_mu.so"))
  #dyn.load('./src/colmeansr.so')
  #}

  if(full.Bayesian){
    mu.samples=matrix(NA,ncol=n,nrow=nsamples)
  }

  meansmc=numeric(n)
  xsumvec = numeric(n)
  x2sumvec =numeric(n)
  for(iter in 1:nsamples){

    res = .C('Rc_mu',mumeans=as.matrix(meansmc,ncol=1),as.double(forcing),as.integer(length(forcing)),
             as.double(hyperpars[iter,1]),as.double(hyperpars[iter,2]),as.double(hyperpars[iter,3]))

    if(full.Bayesian){
      mu.samples[iter,]=res$mumeans
    }

    xsumvec = xsumvec + res$mumeans
    x2sumvec = x2sumvec + res$mumeans^2


  }

  mu.mean = as.numeric(xsumvec/nsamples)

  mu.sd=as.numeric(sqrt( 1/(nsamples-1)*( x2sumvec -2*mu.mean*xsumvec + nsamples*mu.mean^2 ) ))

  if(full.Bayesian){
    mu.quant0.025 = numeric(n)
    mu.quant0.5 = numeric(n)
    mu.quant0.975 = numeric(n)

    for(iter in 1:n){

      dens = density(mu.samples[,iter])
      mu.quant0.025[iter]=INLA::inla.qmarginal(0.025,dens)
      mu.quant0.5[iter]=INLA::inla.qmarginal(0.5,dens)
      mu.quant0.975[iter]=INLA::inla.qmarginal(0.975,dens)
    }

  }

  tid.slutt = proc.time()[[3]]
  tid.mc=tid.slutt-tid.start

  if(print.progress){
    cat("Finished Monte Carlo sampling procedure in ",tid.mc," seconds\n",sep="")
  }

  ret = list(mu.mean=mu.mean, mu.sd = mu.sd)

  if(full.Bayesian){
    ret$mu.quant0.025=mu.quant0.025
    ret$mu.quant0.5=mu.quant0.5
    ret$mu.quant0.975=mu.quant0.975
    ret$samples=list(mu=mu.samples, H=hyperpars[,1],sigmaf=hyperpars[,2],F0=hyperpars[,3])
  }
  ret$time = tid.mc
  return(ret)

}
