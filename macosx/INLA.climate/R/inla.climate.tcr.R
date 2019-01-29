inla.climate.tcr = function(result,Qco2,nsamples=100000,seed=1234,print.progress=FALSE){

  if(print.progress){
    cat("Starting Monte Carlo sampling with n=",nsamples," simulations..\n",sep="")
  }
  set.seed(seed)
  inla.seed = as.integer(runif(1)*.Machine$integer.max)

  #n=result$misc$configs$contents$length[1]

  #x = inla.posterior.sample(nsamples,r,seed=inla.seed) #int.strategy=grid
  x = inla.hyperpar.sample(nsamples,result)
  hyperpars = matrix(NA,nrow=nsamples,ncol=4) #c(H,sf,shift,TCR)
  zmcmc = 1:80

  hyperpars[,1] = 0.5+0.5/(1+exp(-x[,2]))
  hyperpars[,2] = 1/sqrt(exp(x[,3]))
  hyperpars[,3] = x[,4]

  tid.start = proc.time()[[3]]
  #if(!is.loaded('Rc_mu')){
    #print('hallo')
  #  dyn.load(file.path(.Library,"INLA.climate/libs/Rc_mu.so"))
    #dyn.load('./src/colmeansr.so')
  #}


  meansmcmc = numeric(80)
  for(iter in 1:nsamples){
    
    #zzmcmc = hyperpars[iter,2]*(zmcmc+hyperpars[iter,3])
    res = .C('Rc_mu',mumeans=as.matrix(meansmcmc,ncol=1),as.double(zmcmc),as.integer(80),
             as.double(hyperpars[iter,1]),as.double(hyperpars[iter,2]),as.double(hyperpars[iter,3]))
    #mu.eans = mu.cwrapper(as.double(zmcmc),as.integer(80),as.double(hyperpars[iter,1]),
    #                      as.double(hyperpars[iter,2]),as.double(hyperpars[iter,3]))
    # strukturmcmc = (0.5+seq(0,80-1,length.out=80))^(hyperpars[iter,1]-3/2)
    #  for(i in 1:80){
    #    meansmcmc[i] = rev(strukturmcmc[1:i])%*%zzmcmc[1:i]
    #  }
    Tres = Qco2/70 * res$mumeans #meansmcmc

    hyperpars[iter,4] = 1/20*sum(Tres[61:80])

  }
  mcmcfit = density(hyperpars[,4])
  tid.slutt = proc.time()[[3]]
  tid.mcmc=tid.slutt-tid.start

  if(print.progress){
    cat("Finished MCMC sampling procedure in ",tid.mcmc," seconds\n",sep="")
  }


  return(list(TCR.mean=mean(hyperpars[,4]),TCR.sd = sd(hyperpars[,4]),
              TCR.quant0.025=inla.qmarginal(0.025,mcmcfit),
              TCR.quant0.5=inla.qmarginal(0.5,mcmcfit),
              TCR.quant0.975=inla.qmarginal(0.975,mcmcfit),
              samples=list(
   TCR=hyperpars[,4],H=hyperpars[,1],sigmaf=hyperpars[,2],shift=hyperpars[,3]),
   time=tid.mcmc))
}
