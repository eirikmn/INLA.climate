#################################################################

#n_i = b*w_i
rgeneric.ar1 = function(
  cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL)
{
  
  require("INLA.climate",quietly=TRUE)
  
  tau = exp(15)
  envir = environment(sys.call()[[1]])
  
  interpret.theta = function() {
    if(!is.null(envir)){
      NN=get("N",envir)
    }
    
    kappax = exp(theta[1])
    kappaf = exp(theta[2])
    a = 3
    F0=theta[3]#F0 = -a + 2*a/(1+exp(-theta[3]))
    para=data.frame(kappax=kappax,kappaf=kappaf,F0=F0)
    denom = sum(exp(c(0,theta[3+(1:(NN-1))])))
    
    weights=numeric(NN)
    if(NN>1){
      for(i in 2:(NN)){
        weights[i] <- exp(theta[2+i])
      }
    }
    
    weights[1]=1
    weights = weights/sum(weights)
    for(i in 1:NN){
      para[[paste0("w",i)]] <- weights[i]
    }
    
    for(i in 1:NN){
      para[[paste0("p",i)]] <- 1/(1+sum(exp(-theta[ (NN+2)+1:i ])))
    }
    
    
    
    return(para)
  }
  
  mu = function() {
    
    if(!is.null(envir)){
      nn=get("n",envir)
      NN=get("N",envir)
      fforcing=get("forcing",envir)
      
      
    }
    hyperparam = interpret.theta()
    #print(hyperparam)
    
    sf = 1/sqrt(hyperparam$kappaf)
    means = numeric(nn)
    if(NN == 1){
      weights = 1
    }else{
      weights = as.numeric(hyperparam)[3+1:NN]
    }
    
    ######
    pp=as.numeric(hyperparam)[NN+3+1:NN]
    llambdas = pp-1
    
    if(!is.loaded('Rc_mu_ar1')){
      #dyn.load(file.path(.Library,"INLA.climate/libs/Rc_Q.so"))
      dyn.load(file.path("Rc_mu_ar1.so"))
    }
    res = .C('Rc_mu_ar1',mu=as.matrix(means,ncol=1),as.double(fforcing),as.integer(nn),as.integer(NN),
             as.double(weights),as.double(llambdas),as.double(sf),
             as.double(hyperparam$F0))
    #print(res$mu[1:3])
    return(c(res$mu,rep(0,NN*nn)))
  }
  
  
  graph = function()
  {
    
    if(!is.null(envir)){
      nn=get("n",envir)
      NN=get("N",envir)
      
    }else{
      nn=get("n",environment())
      NN=get("N",environment())
      
    }
    
    ii = numeric(2.5*NN*nn+nn-NN+nn*NN*NN/2)
    jj = numeric(2.5*NN*nn+nn-NN+nn*NN*NN/2)
    xx = rep(1,2.5*NN*nn+nn-NN+nn*NN*NN/2)
    
    res = .C('Rc_Q',minii=as.double(ii),minjj=as.double(jj),minxx=as.double(xx),
             as.integer(nn),as.integer(NN),as.double(rep(1/NN,NN)),as.double(rep(0.5,NN)),
             as.double(tau),as.double(1.0))
    
    G = Matrix::sparseMatrix(i=res$minii,j=res$minjj,x=res$minxx,symmetric=TRUE)
    G[G != 0] = 1
    return (G)
  }
  
  Q = function()
  {
    if(!is.null(envir)){
      nn=get("n",envir)
      NN=get("N",envir)
      
      
    }
    
    hyperparam = interpret.theta()
    
    pparam = hyperparam[NN+3+(1:NN)]
    sx = 1/sqrt(hyperparam$kappax)
    alphas = pparam[1:NN]
    if(NN > 1){
      weights = hyperparam[3+(1:NN)]
    }else{
      weights = 1
    }
    
    
    ii = numeric(2.5*NN*nn+nn-NN+nn*NN*NN/2)
    jj = numeric(2.5*NN*nn+nn-NN+nn*NN*NN/2)
    xx = numeric(2.5*NN*nn+nn-NN+nn*NN*NN/2)
    if(!is.loaded('Rc_Q')){
      #dyn.load(file.path(.Library,"INLA.climate/libs/Rc_Q.so"))
      dyn.load(file.path("Rc_Q.so"))
    }
    if(length(theta)==0){
      sx=1;weights=rep(1/NN,NN);alphas=weights
    }
    
    res = .C('Rc_Q',minii=as.double(ii),minjj=as.double(jj),minxx=as.double(xx), #skal ikke sigma inn her?
             as.integer(nn),as.integer(NN),as.double(weights),as.double(alphas),
             as.double(tau),as.double(sx))
    
    
    Q = Matrix::sparseMatrix(i=res$minii,j=res$minjj,x=res$minxx,symmetric=TRUE)
    
    return ( Q )
  }
  
  
  log.norm.const = function()
  {
    if(!is.null(envir)){
      nn=get("n",envir)
      NN=get("N",envir)
    }
    # tid.rgen.start = proc.time()[[3]]
    hyperparams = interpret.theta()
    pparam = hyperparams[NN+3+(1:NN)]
    
    #N = length(param)/2
    
    sum = nn/2*log(tau)
    for (i in 1:NN){
      sum = sum -(nn-1)/2*log(1-pparam[i]^2)
    }
    tid.rgen.slutt = proc.time()[[3]]
    
    return(sum)
  }
  
  log.prior = function()
  {
    if(!is.null(envir)){
      
      NN=get("N",envir)
      #pparam=get("params",envir)
    }
    # tid.rgen.start = proc.time()[[3]]
    #print("prior")
    params = interpret.theta()
    #H = params$H
    a = 1
    b = 0.01
    aa = 1
    bb = 0.01
    lprior = INLA::inla.pc.dprec(params$kappax, u=a, alpha=b, log=TRUE) + log(params$kappax)
    lprior = lprior + INLA::inla.pc.dprec(params$kappaf, u=aa, alpha=bb, log=TRUE) + log(params$kappaf)
    #lprior = lprior + log(0.5)+log(1+1/(1+exp(-theta[2]))) - theta[2]-2*log(1+exp(-theta[2]))
    a=3
    
    #skal ikke interne variabler brukes her? og hva så med w vektene?
    shift.a = 3
    #lprior = lprior + dnorm(theta[4],log=TRUE)
    lprior = lprior + dnorm(theta[3],log=TRUE)
    #lprior = lprior + dnorm(-shift.a+2*shift.a/(1+exp(-params$F0)),sd=0.2,log=TRUE)+log(2*shift.a)-params$F0 -2*log(1+exp(-params$F0))
    if(NN==1){
      lprior = lprior + dnorm(theta[4],log=TRUE)
      return(lprior)
    }
    for(m in 1:(NN-1)){
      lprior = lprior + dnorm(theta[3+m],log=TRUE)
      #dgamma(theta[3+m],shape=1,log=T)#dnorm(params[[3+m]],log=TRUE)
    }
    for(m in 1:(NN)){
      lprior = lprior + dnorm(theta[2+NN+m],log=TRUE)
      
    }
    
    
    return (lprior)
  }
  
  initial = function()
  {
    if(!is.null(envir)){
      NN=get("N",envir)
    }
    ini=c(-3,0.,0.)
    if(NN==1){
      ini = c(ini,0)
      return(ini)
    }else{
      for(i in 1:(NN-1)){
        ini=c(ini,0.)
      }
      for(i in 1:NN){
        ini=c(ini,0.)
      }
    }
    
    return (ini)
  }
  
  quit = function()
  {
    return ()
  }
  # if(is.null(theta)){
  #   theta = initial()
  #   #envir=.GlobalEnv
  #   
  # }
  
  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
}

# sim1 = rnorm(100000)
# sim2 = rnorm(100000)
# w0sim = exp(0)/(1+exp(sim1)+exp(sim2))
# w1sim = exp(sim1)/(1+exp(sim1)+exp(sim2))
# w2sim = exp(sim2)/(1+exp(sim1)+exp(sim2))
# hist(w1sim)
