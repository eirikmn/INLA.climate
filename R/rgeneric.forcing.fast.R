

#################################################################

#n_i = b*w_i
rgeneric.forcing.fast = function(
  cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL)
{
  
  require("INLA.climate",quietly=TRUE)
  
  tau = exp(15)
  envir = environment(sys.call()[[1]])
  

  Hmap = function(H) {

    params = numeric(2*N)
    for(i in 1:(2*N)){
      params[i] = funks[[i]](H)
    }

    return(params)
  }

  interpret.theta = function() {

    kappa = exp(theta[1])
    H = 0.5 + 0.5 / (1 + exp(-theta[2]))
    scale = exp(theta[3])
    a = 3
    shift = -a + 2*a/(1+exp(-theta[4]))

    return(list(H = H, kappa = kappa, scale = scale, shift = shift))
  }

  mu = function() {
    # tid.rgen.start = proc.time()[[3]]
    hyperparam = interpret.theta()
    H = hyperparam$H
    sf = 1/sqrt(hyperparam$scale)
    shift = hyperparam$shift
    #n = 30
    means = numeric(n)

    #if(!is.loaded('Rc_mu')){
      #print('hallo')
    #  dyn.load(file.path(.Library,"INLA.climate/libs/Rc_mu.so"))
      #dyn.load('./src/colmeansr.so')
    #}
    #means = mu.cwrapper(as.double(forcing),as.integer(n),as.double(H),as.double(sf),as.double(shift))
    res = .C('Rc_mu',ans=as.matrix(means,ncol=1),as.double(forcing),as.integer(n),
             as.double(H),as.double(sf),as.double(shift))


    return(c(res$ans,rep(0,N*n)))
  }

  ar1maker = function(rho) {
    # tid.rgen.start = proc.time()[[3]]
    tauu = 1 / (1 - rho ^ 2)
    i = c(1L, n, 2L:(n - 1L), 1L:(n - 1L))
    j = c(1L, n, 2L:(n - 1L), 2L:n)
    xx = #tauu *
      c(1/(1-rho^2), 1/(1-rho^2), rep((1 + rho ^ 2)/(1-rho^2), n - 2L),
        rep(-rho/(1-rho^2), n - 1L))
    Q = Matrix::sparseMatrix(
      i = i,
      j = j,
      x = xx,
      giveCsparse = FALSE,
      symmetric = TRUE
    )

    return(Q)
  }

  graph = function()
  {
    G = Q()
    G[G != 0] = 1

    return (G)
  }

  Q = function()
  {
    hyperparam = interpret.theta()
    H = hyperparam$H

    kappa = hyperparam$kappa
    param = Hmap(hyperparam$H)

    sx = 1/sqrt(hyperparam$kappa)
    alphas = param[1:N]
    weights = param[(N+1):(2*N)]


     ii = numeric(2.5*N*n+n-N+n*N*N/2)
     jj = numeric(2.5*N*n+n-N+n*N*N/2)
     xx = numeric(2.5*N*n+n-N+n*N*N/2)
     #if(!is.loaded('Rc_Q')){
     #  dyn.load(file.path(.Library,"INLA.climate/libs/Rc_Q.so"))
     #}
  #print("INLA.climate"%in%.packages())
  #h.map.maker(4,1000,"fgn")
    #
    if(length(theta)==0){
      sx=1;weights=rep(1/N,N);alphas=weights
    }
     #res = Q.cwrapper(as.integer(n),as.double(weights),as.double(alphas),
    #                  as.integer(N),as.double(tau),as.double(sx))
    #res = .C('myQr',minii=as.double(ii),minjj=as.double(jj),minxx=as.double(xx),
    res = .C('Rc_Q',minii=as.double(ii),minjj=as.double(jj),minxx=as.double(xx),
              as.integer(n),as.integer(N),as.double(weights),as.double(alphas),
              as.double(tau),as.double(sx), PACKAGE="INLA.climate")



    Q = Matrix::sparseMatrix(i=res$minii,j=res$minjj,x=res$minxx,symmetric=TRUE)

    return (Q)
  }

  log.norm.const = function()
  {
    # tid.rgen.start = proc.time()[[3]]
    hyperparams = interpret.theta()
    param = Hmap(hyperparams$H)
    #N = length(param)/2

    sum = n/2*log(tau)
    for (i in 1:N){
      sum = sum -(n-1)/2*log(1-param[i]^2)
    }
    tid.rgen.slutt = proc.time()[[3]]

    return(sum)
  }

  log.prior = function()
  {
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
    lprior = lprior + lprior.fun.H(theta[2])
    #lprior = lprior + log(0.5)+log(1+1/(1+exp(-theta[2]))) - theta[2]-2*log(1+exp(-theta[2]))
    a=3
    lprior = lprior + dnorm(-a+2*a/(1+exp(-params$shift)),sd=0.2,log=TRUE)+log(2*a)-params$shift -2*log(1+exp(-params$shift))

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
