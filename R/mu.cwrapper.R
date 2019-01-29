#   #' @useDynLib INLA.climate cmur

mu.cwrapper = function(forcing,n,H,sf,shift){
  means = as.matrix(numeric(n),ncol=1)
  #dyn.load(file.path(.Library,"INLA.climate/libs/cmur.so"))
  res = .C('Rc_mu',mumeans=means,as.double(forcing),as.integer(n),
           as.double(H),as.double(sf),as.double(shift))
  return(res$mumeans)
}

