inla.climate.mu = function(result,forcing,quick=FALSE,T0.corr=0,nsamples=100000,seed=1234,
                           print.progress=FALSE,model="fgn"){

  catch = tryCatch(attachNamespace("INLA"),error=function(x){})
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
  if(is.null(T0.corr)){
    if(!is.null(result$climate.misc$T0)){
      T0.corr = result$climate.misc$T0
    }else{
      T0.corr=0
    }
  }
  
  
  
  if(print.progress){
    cat("Starting mu Monte Carlo sampling with n=",format(nsamples,scientific=F)," simulations..\n",sep="")
  }
  set.seed(seed)
  inla.seed = as.integer(runif(1)*.Machine$integer.max)

  #n=climate.res$misc$configs$contents$length[1]
  n=length(forcing)

  #x = inla.posterior.sample(nsamples,r,seed=inla.seed) #int.strategy=grid
  x = INLA::inla.hyperpar.sample(nsamples,climate.res)
  if(dim(x)[2]>4){
    model = "ar1"
  }
  if(model %in% c("fgn","arfima")){
    hyperpars = matrix(NA,nrow=nsamples,ncol=4) #c(H,sf,F0,TCR)
    hyperpars[,1] = 0.5+0.5/(1+exp(-x[,2]))
    hyperpars[,2] = 1/sqrt(exp(x[,3]))
    a=3
    hyperpars[,3] = x[,4]
    #hyperpars[,3] = -a+2*a/(1+exp(-x[,4]))
  }else{
    hyperpars = matrix(NA,nrow=nsamples,ncol=(2)) #c(H,sf,F0,w1,...,wm,L1,...,Lm)
    hyperpars[,1] = 1/sqrt(exp(x[,2]))
    a=3
    hyperpars[,2] = -a+2*a/(1+exp(-x[,3]))
    m = (dim(x)[2]-2)/2
    ar1.temp= inla.climate.ar1(climate.res,m=m,nsamples=nsamples,seed=seed,print.progress=print.progress)
    ww = matrix(NA,nrow=nsamples,ncol=m)
    LL = matrix(NA,nrow=nsamples,ncol=m)
    if(m == 1){
      ww = rep(1,nsamples)
      LL = inla.rmarginal(nsamples,ar1.temp$p$density)-1 #lambda
    }else{
      for(k in 1:m){
        ww[,k] = ar1.temp[[paste0("w",k)]]$samples
        LL[,k] = ar1.temp[[paste0("p",k)]]$samples-1 #lambda
      }
    }
    
  }
  


  

   tid.start = proc.time()[[3]]
  #if(!is.loaded('Rc_mu')){
  #print('hallo')
  #  dyn.load(file.path(.Library,"INLA.climate/libs/Rc_mu.so"))
  #dyn.load('./src/colmeansr.so')
  #}

  if(!quick){
    mu.samples=matrix(NA,ncol=n,nrow=nsamples)
  }

  meansmc=numeric(n)
  xsumvec = numeric(n)
  x2sumvec =numeric(n)
  for(iter in 1:nsamples){
    if(model %in% c("fgn","arfima")){
      
      res = .C('Rc_mu',mumeans=as.matrix(meansmc,ncol=1),as.double(forcing),as.integer(length(forcing)),
               as.double(hyperpars[iter,1]),as.double(hyperpars[iter,2]),as.double(hyperpars[iter,3]))
      
    }else if(model == "ar1"){
      
      if(!is.loaded('Rc_mu_ar1')){
        #dyn.load(file.path(.Library,"INLA.climate/libs/Rc_Q.so"))
        dyn.load(file.path("Rc_mu_ar1.so"))
      }
      if(m == 1){
        res = .C('Rc_mu_ar1',mumeans=as.matrix(meansmc,ncol=1),as.double(forcing),as.integer(length(forcing)),as.integer(m),
                 as.double(1),as.double(LL[iter]),as.double(hyperpars[iter,1]),
                 as.double(hyperpars[iter,2]))
      }else{
        res = .C('Rc_mu_ar1',mumeans=as.matrix(meansmc,ncol=1),as.double(forcing),as.integer(length(forcing)),as.integer(m),
                 as.double(ww[iter,]),as.double(LL[iter,]),as.double(hyperpars[iter,1]),
                 as.double(hyperpars[iter,2]))
      }
      
    }
    

    if(!quick){
      mu.samples[iter,]=res$mumeans
    }
    xsumvec = xsumvec + res$mumeans
    x2sumvec = x2sumvec + res$mumeans^2

  }

  mu.mean = as.numeric(xsumvec/nsamples)

  mu.sd=as.numeric(sqrt( 1/(nsamples-1)*( x2sumvec -2*mu.mean*xsumvec + nsamples*mu.mean^2 ) ))

  if(!quick){
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
    cat("Finished mu Monte Carlo sampling procedure in ",tid.mc," seconds\n",sep="")
  }

  ret = list(mean=mu.mean+T0.corr, sd = mu.sd)

  if(!quick){
    ret$quant0.025=mu.quant0.025+T0.corr
    ret$quant0.5=mu.quant0.5+T0.corr
    ret$quant0.975=mu.quant0.975+T0.corr
    if(model %in% c("fgn","arfima")){
      ret$samples=list(mu=mu.samples+T0.corr, H=hyperpars[,1],sigmaf=hyperpars[,2],F0=hyperpars[,3])
    }else if(model == "ar1"){
      ret$samples=list(mu=mu.samples+T0.corr, sigmaf=hyperpars[,1],F0=hyperpars[,2])
      if(m==1){
        ret$samples$p = LL+1
      }else{
        for(k in 1:m){
          ret$samples[[paste0("w",k)]] = ww[,k]
        }
        for(k in 1:m){
          ret$samples[[paste0("p",k)]] = LL[,k]+1
        }
      }
      
    }
    
  }
  
  if(class(result) == "inla.climate"){
    if(print.progress){
      print("Exporting inla.climate object")
    }
    result$mu = ret
    result$time$mu = tid.mc
    result$time$Total = result$time$Total + tid.mc
    result$misc$mu.options$nsamples = nsamples
    result$misc$mu.options$seed = seed
    if(quick){
      compute.mu = 2
    }else{
      compute.mu=1
    }
    result$misc$mu.options$compute.mu = compute.mu
    return(result)
  }else{
    if(print.progress){
      print("Exporting list object")
    }
    ret$time = tid.mc
    return(ret)
  }
  
}
