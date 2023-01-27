

#################################################################

#n_i = b*w_i
rgeneric.lrd = function(
  cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL)
{
  
  require("INLA.climate",quietly=TRUE)
  
  tau = exp(15)
  envir = environment(sys.call()[[1]])
  

  Hmap = function(H) {

    if(!is.null(envir)){
      NN=get("N",envir)
      ffunks=get("funks",envir)
      
    }
    params = numeric(2*NN)
    for(i in 1:(2*NN)){
      params[i] = ffunks[[i]](H)
    }

    return(params)
  }

  interpret.theta = function() {

    kappa = exp(theta[1])
    H = 0.5 + 0.5 / (1 + exp(-theta[2]))
    scale = exp(theta[3])
    
    #a = 3
    #shift = -a + 2*a/(1+exp(-theta[4]))
    shift=theta[4]
    
    return(list(H = H, kappa = kappa, scale = scale, shift = shift))
  }

  mu = function() {
    
    # tid.rgen.start = proc.time()[[3]]
    if(!is.null(envir)){
      nn=get("n",envir)
      NN=get("N",envir)
      fforcing=get("forcing",envir)
    }
    hyperparam = interpret.theta()
    H = hyperparam$H
    sf = 1/sqrt(hyperparam$scale)
    shift = hyperparam$shift
    #n = 30
    means = numeric(nn)

   
    res = .C('Rc_mu',ans=as.matrix(means,ncol=1),as.double(fforcing),as.integer(nn),
             as.double(H),as.double(sf),as.double(shift), PACKAGE="INLA.climate")


    return(c(res$ans,rep(0,NN*nn)))
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
             as.double(tau),as.double(1.0), PACKAGE="INLA.climate")
    
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
    H = hyperparam$H

    kappa = hyperparam$kappa
    param = Hmap(hyperparam$H)

    sx = 1/sqrt(hyperparam$kappa)
    alphas = param[1:NN]
    weights = param[(NN+1):(2*NN)]


     ii = numeric(2.5*NN*nn+nn-NN+nn*NN*NN/2)
     jj = numeric(2.5*NN*nn+nn-NN+nn*NN*NN/2)
     xx = numeric(2.5*NN*nn+nn-NN+nn*NN*NN/2)
     #if(!is.loaded('Rc_Q')){
     #  dyn.load(file.path(.Library,"INLA.climate/libs/Rc_Q.so"))
     #}
  #print("INLA.climate"%in%.packages())
  #h.map.maker(4,1000,"fgn")
    #
    if(length(theta)==0){
      sx=1;weights=rep(1/NN,NN);alphas=weights
    }
     #res = Q.cwrapper(as.integer(n),as.double(weights),as.double(alphas),
    #                  as.integer(N),as.double(tau),as.double(sx))
    #res = .C('myQr',minii=as.double(ii),minjj=as.double(jj),minxx=as.double(xx),
    res = .C('Rc_Q',minii=as.double(ii),minjj=as.double(jj),minxx=as.double(xx),
              as.integer(nn),as.integer(NN),as.double(weights),as.double(alphas),
              as.double(tau),as.double(sx), PACKAGE="INLA.climate")



    Q = Matrix::sparseMatrix(i=res$minii,j=res$minjj,x=res$minxx,symmetric=TRUE)

    return (Q)
  }

  log.norm.const = function()
  {
    if(!is.null(envir)){
      nn=get("n",envir)
      NN=get("N",envir)
    }
    # tid.rgen.start = proc.time()[[3]]
    hyperparams = interpret.theta()
    param = Hmap(hyperparams$H)
    #N = length(param)/2

    sum = nn/2*log(tau)
    for (i in 1:NN){
      sum = sum -(nn-1)/2*log(1-param[i]^2)
    }
    tid.rgen.slutt = proc.time()[[3]]

    return(sum)
  }

  log.prior = function()
  {
    if(!is.null(envir)){
      
      llprior.fun.H=get("lprior.fun.H",envir)
    }
    # tid.rgen.start = proc.time()[[3]]
    #print("prior")
    params = interpret.theta()
    #H = params$H
    a = 1
    b = 0.01
    aa = 1
    bb = 0.01
    lprior = INLA::inla.pc.dprec(params$kappa, u=a, alpha=b, log=TRUE) + log(params$kappa)
    lprior = lprior + INLA::inla.pc.dprec(params$scale, u=aa, alpha=bb, log=TRUE) + log(params$scale)
    lprior = lprior + llprior.fun.H(theta[2])
    #lprior = lprior + log(0.5)+log(1+1/(1+exp(-theta[2]))) - theta[2]-2*log(1+exp(-theta[2]))
    a=3
    #lprior = lprior + dnorm(-a+2*a/(1+exp(-params$shift)),sd=0.2,log=TRUE)+log(2*a)-params$shift -2*log(1+exp(-params$shift))
    lprior = lprior + dnorm(theta[4],log=TRUE)
    
    return (lprior)
  }

  initial = function()
  {
    ini = c(-3,0.,-10,0.5)
    return (ini)
  }

  quit = function()
  {
    return ()
  }
  if(is.null(theta)){
    theta = initial()
    #envir=.GlobalEnv
    
  }

  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
}
