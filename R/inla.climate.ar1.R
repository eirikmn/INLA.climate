inla.climate.ar1 = function(result,m=1,nsamples=100000,seed=1234,print.progress=FALSE){
  atch = tryCatch(attachNamespace("INLA"),error=function(x){})
  if(length(find.package("INLA",quiet=TRUE))==0){
    stop("This function requires INLA. Please install at www.R-INLA.org or by calling 'install.packages(\"INLA\", repos=c(getOption(\"repos\"), INLA=\"https://inla.r-inla-download.org/R/testing\"), dep=TRUE)' from R.")
  }
  
  if(class(result)=="inla.climate"){
    climate.res = result$inla.result
  }else if(class(result)=="inla"){
    climate.res = result
  }else{
    stop("Input 'result' not a valid class.")
  }
  
  
  if(print.progress){
    cat("Starting ar1 weights and first-lag correlation parameter Monte Carlo sampling with n = ",format(nsamples,scientific=F)," simulations..\n",sep="")
  }
  tid.start = proc.time()[[3]]
  set.seed(seed)
  inla.seed = as.integer(runif(1)*.Machine$integer.max)

  #n=climate.res$misc$configs$contents$length[1]

  
  #x = inla.posterior.sample(nsamples,r,seed=inla.seed) #int.strategy=grid
  
  if(m == 1){
    ww = 1
    pp = inla.tmarginal(function(x)1/(1+exp(-x)),climate.res$marginals.hyperpar$`Theta4 for idy`)
    zpp = inla.zmarginal(pp,silent=TRUE)
    ret = list(w = 1,
               p = list(mean = zpp$mean, sd = zpp$sd,
                         quant0.025=zpp$quant0.025,
                         quant0.5=zpp$quant0.5,
                         quant0.975=zpp$quant0.975,
                        density = pp) )
  }else{
    x = INLA::inla.hyperpar.sample(nsamples,climate.res)
    vv = cbind(rep(1,nsamples),x[,4:(2+m)])
    uu = x[,(2+m+1):(2*m+2)]
    ww=matrix(NA,nrow=nsamples,ncol=m)
    pp=matrix(NA,nrow=nsamples,ncol=m)
    for(k in 1:m){
      ww[,k] = exp(vv[,k])
      pp[,k] = 1/(1+rowSums(exp(-as.matrix(uu[,1:k]))))
    }
    ww = ww/rowSums(ww)
    
    
    ret = list()
    for(k in 1:m){
      mcfit = density(ww[,k])
      ret[[paste0("w",k)]] = list(mean=mean(ww[,k]), sd = sd(ww[,k]),
                                  quant0.025 = INLA::inla.qmarginal(0.025,mcfit),
                                  quant0.5 = INLA::inla.qmarginal(0.5,mcfit),
                                  quant0.975 = INLA::inla.qmarginal(0.975,mcfit),
                                  samples = ww[,k])
    }
    for(k in 1:m){
      mcfit = density(pp[,k])
      ret[[paste0("p",k)]] = list(mean=mean(pp[,k]), sd = sd(pp[,k]),
                                  quant0.025 = INLA::inla.qmarginal(0.025,mcfit),
                                  quant0.5 = INLA::inla.qmarginal(0.5,mcfit),
                                  quant0.975 = INLA::inla.qmarginal(0.975,mcfit),
                                  samples = pp[,k])
    }
  }
  
  tid.slutt = proc.time()[[3]]
  tid.mc=tid.slutt-tid.start

  if(print.progress){
    cat("Finished Monte Carlo sampling procedure in ",tid.mc," seconds\n",sep="")
  }

  
  if(class(result) == "inla.climate"){
    if(print.progress){
      print("Exporting inla.climate object")
    }
    result$ar1 = ret
    result$time$ar1 = tid.mc
    result$misc$ar1.options$nsamples = nsamples
    result$misc$ar1.options$seed = seed
    return(result)
  }else{
    if(print.progress){
      print("Exporting list object")
    }
    
    ret$time = tid.mc
    return(ret)
  }
}
