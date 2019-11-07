inla.climate.tcr = function(result,Qco2,nsamples=100000,seed=1234,
                            print.progress=FALSE,model="fgn"){
  
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
    cat("Starting TCR Monte Carlo sampling with n = ",format(nsamples,scientific=F)," simulations..\n",sep="")
  }
  set.seed(seed)
  inla.seed = as.integer(runif(1)*.Machine$integer.max)

  #n=climate.res$misc$configs$contents$length[1]

  
  #x = inla.posterior.sample(nsamples,r,seed=inla.seed) #int.strategy=grid
  x = INLA::inla.hyperpar.sample(nsamples,climate.res)
  if(dim(x)[2]>4){
    model = "ar1"
  }
  zmc = 1:80
  a=3
  
  if(model %in% c("fgn","arfima")){
    tcr.col = 4
    hyperpars = matrix(NA,nrow=nsamples,ncol=tcr.col) #c(H,sf,shift,TCR)
    
    hyperpars[,1] = 0.5+0.5/(1+exp(-x[,2]))
    hyperpars[,2] = 1/sqrt(exp(x[,3]))
    hyperpars[,3] = x[,4]#-a+2*a/(1+exp(-x[,4]))
  }else{
    tcr.col=3
    hyperpars = matrix(NA,nrow=nsamples,ncol=tcr.col) #c(H,sf,F0,w1,...,wm,L1,...,Lm)
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


  meansmc = numeric(80)
  for(iter in 1:nsamples){
    if(model %in% c("fgn","arfima")){
      #zzmc = hyperpars[iter,2]*(zmc+hyperpars[iter,3])
      res = .C('Rc_mu',mumeans=as.matrix(meansmc,ncol=1),as.double(zmc),as.integer(80),
               as.double(hyperpars[iter,1]),as.double(hyperpars[iter,2]),as.double(hyperpars[iter,3]))
    }else if(model == "ar1"){
      if(!is.loaded('Rc_mu_ar1')){
        #dyn.load(file.path(.Library,"INLA.climate/libs/Rc_Q.so"))
        dyn.load(file.path("Rc_mu_ar1.so"))
      }
      if(m == 1){
        res = .C('Rc_mu_ar1',mumeans=as.matrix(meansmc,ncol=1),as.double(zmc),as.integer(80),as.integer(m),
                 as.double(1),as.double(LL[iter]),as.double(hyperpars[iter,1]),
                 as.double(hyperpars[iter,2]))
      }else{
        res = .C('Rc_mu_ar1',mumeans=as.matrix(meansmc,ncol=1),as.double(zmc),as.integer(80),as.integer(m),
                 as.double(ww[iter,]),as.double(LL[iter,]),as.double(hyperpars[iter,1]),
                 as.double(hyperpars[iter,2]))
      }
      
      
    }
    
    #mu.eans = mu.cwrapper(as.double(zmc),as.integer(80),as.double(hyperpars[iter,1]),
    #                      as.double(hyperpars[iter,2]),as.double(hyperpars[iter,3]))
    # strukturmc = (0.5+seq(0,80-1,length.out=80))^(hyperpars[iter,1]-3/2)
    #  for(i in 1:80){
    #    meansmc[i] = rev(strukturmc[1:i])%*%zzmc[1:i]
    #  }
    Tres = Qco2/70 * res$mumeans #meansmc

    hyperpars[iter,tcr.col] = 1/20*sum(Tres[61:80])

  }
  mcfit = density(hyperpars[,tcr.col])
  tid.slutt = proc.time()[[3]]
  tid.mc=tid.slutt-tid.start

  if(print.progress){
    cat("Finished TCR Monte Carlo sampling procedure in ",tid.mc," seconds\n",sep="")
  }

  if(model %in% c("arfima","fgn")){
    ret = list(mean=mean(hyperpars[,tcr.col]),sd = sd(hyperpars[,tcr.col]),
               quant0.025=INLA::inla.qmarginal(0.025,mcfit),
               quant0.5=INLA::inla.qmarginal(0.5,mcfit),
               quant0.975=INLA::inla.qmarginal(0.975,mcfit),
               samples=list(
                 TCR=hyperpars[,tcr.col],H=hyperpars[,1],sigmaf=hyperpars[,2],shift=hyperpars[,3]))
  }else if(model == "ar1"){
    if(m==1){
      ret = list(mean=mean(hyperpars[,tcr.col]),sd = sd(hyperpars[,tcr.col]),
                 quant0.025=INLA::inla.qmarginal(0.025,mcfit),
                 quant0.5=INLA::inla.qmarginal(0.5,mcfit),
                 quant0.975=INLA::inla.qmarginal(0.975,mcfit),
                 samples=list(
                   TCR=hyperpars[,tcr.col],p=LL+1,sigmaf=hyperpars[,1],shift=hyperpars[,2]))
    }else{
      ret = list(mean=mean(hyperpars[,tcr.col]),sd = sd(hyperpars[,tcr.col]),
                 quant0.025=INLA::inla.qmarginal(0.025,mcfit),
                 quant0.5=INLA::inla.qmarginal(0.5,mcfit),
                 quant0.975=INLA::inla.qmarginal(0.975,mcfit),
                 samples=list(
                   TCR=hyperpars[,tcr.col],sigmaf=hyperpars[,1],shift=hyperpars[,2]))
      for(k in 1:m){
        ret$samples[[paste0("w",k)]] = ww[,k]
      }
      for(k in 1:m){
        ret$samples[[paste0("p",k)]] = LL[,k]+1
      }
    }
  }
  
  
  if(class(result) == "inla.climate"){
    if(print.progress){
      print("Exporting inla.climate object")
    }
    result$TCR = ret
    result$time$TCR = tid.mc
    result$time$Total = result$time$Total + tid.mc
    result$misc$TCR.options$nsamples = nsamples
    result$misc$TCR.options$seed = seed
    result$misc$TCR.options$Qco2 = Qco2
    return(result)
  }else{
    if(print.progress){
      print("Exporting list object")
    }
    
    ret$time = tid.mc
    return(ret)
  }
}
