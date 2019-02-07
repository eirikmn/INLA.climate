inla.climate.tcr = function(result,Qco2,nsamples=100000,seed=1234,print.progress=FALSE){
  
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
    cat("Starting Monte Carlo sampling with n = ",format(nsamples,scientific=F)," simulations..\n",sep="")
  }
  set.seed(seed)
  inla.seed = as.integer(runif(1)*.Machine$integer.max)

  #n=climate.res$misc$configs$contents$length[1]

  
  #x = inla.posterior.sample(nsamples,r,seed=inla.seed) #int.strategy=grid
  x = INLA::inla.hyperpar.sample(nsamples,climate.res)
  hyperpars = matrix(NA,nrow=nsamples,ncol=4) #c(H,sf,shift,TCR)
  zmc = 1:80

  hyperpars[,1] = 0.5+0.5/(1+exp(-x[,2]))
  hyperpars[,2] = 1/sqrt(exp(x[,3]))
  hyperpars[,3] = x[,4]

  tid.start = proc.time()[[3]]
  #if(!is.loaded('Rc_mu')){
    #print('hallo')
  #  dyn.load(file.path(.Library,"INLA.climate/libs/Rc_mu.so"))
    #dyn.load('./src/colmeansr.so')
  #}


  meansmc = numeric(80)
  for(iter in 1:nsamples){

    #zzmc = hyperpars[iter,2]*(zmc+hyperpars[iter,3])
    res = .C('Rc_mu',mumeans=as.matrix(meansmc,ncol=1),as.double(zmc),as.integer(80),
             as.double(hyperpars[iter,1]),as.double(hyperpars[iter,2]),as.double(hyperpars[iter,3]))
    #mu.eans = mu.cwrapper(as.double(zmc),as.integer(80),as.double(hyperpars[iter,1]),
    #                      as.double(hyperpars[iter,2]),as.double(hyperpars[iter,3]))
    # strukturmc = (0.5+seq(0,80-1,length.out=80))^(hyperpars[iter,1]-3/2)
    #  for(i in 1:80){
    #    meansmc[i] = rev(strukturmc[1:i])%*%zzmc[1:i]
    #  }
    Tres = Qco2/70 * res$mumeans #meansmc

    hyperpars[iter,4] = 1/20*sum(Tres[61:80])

  }
  mcfit = density(hyperpars[,4])
  tid.slutt = proc.time()[[3]]
  tid.mc=tid.slutt-tid.start

  if(print.progress){
    cat("Finished Monte Carlo sampling procedure in ",tid.mc," seconds\n",sep="")
  }

  ret = list(TCR.mean=mean(hyperpars[,4]),TCR.sd = sd(hyperpars[,4]),
             TCR.quant0.025=INLA::inla.qmarginal(0.025,mcfit),
             TCR.quant0.5=INLA::inla.qmarginal(0.5,mcfit),
             TCR.quant0.975=INLA::inla.qmarginal(0.975,mcfit),
             samples=list(
               TCR=hyperpars[,4],H=hyperpars[,1],sigmaf=hyperpars[,2],shift=hyperpars[,3]),
             time=tid.mc)
  
  return(ret)
}
