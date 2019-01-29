

Q.cwrapper = function(n,weights,alphas,N,tau,sx){
  means = as.matrix(numeric(n),ncol=1)
  ii = numeric(2.5*N*n+n-N+n*N*N/2)   
  jj = numeric(2.5*N*n+n-N+n*N*N/2)   
  xx = numeric(2.5*N*n+n-N+n*N*N/2)   
  #if(!is.loaded('cQ')){   
  #  dyn.load(file.path(.Library,"INLA.climate/libs/Rc_Q.so"))
  #}
  res = .C('Rc_Q',minii=as.double(ii),minjj=as.double(jj),minxx=as.double(xx),
           as.integer(n),as.integer(N),as.double(weights),as.double(alphas),
           as.double(tau),as.double(sx))
  ret = list(minii=res$minii,minjj=res$minjj,minxx=res$minxx)
  return(ret)
}
