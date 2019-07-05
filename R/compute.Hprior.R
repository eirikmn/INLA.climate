compute.Hprior = function(n, upper,  alpha, persistent=TRUE, model="fgn") {
  #cat("\n\nCOMPUTE PRIOR!!!\n\n\n")
  b=0.999  
  if(persistent){
    a=0.501
    rho.seq=seq(log(a-0.5)-log(1-a),log(b-0.5)-log(1-b),length=500)
    h=0.5+0.5/(1+exp(-rho.seq))
  }else{
    a=0.001
    rho.seq=seq(log(a)-log(1-a),log(b)-log(1-b),length=1000)
    h=1/(1+exp(-rho.seq))
  }
  
  #prior=func.priorH.new(n,h,upper,alpha,persistent)
  a.func.corr=rep(0,length(h))
  kk=1:(n-1)
  for (i in 1:length(h)){
    if(model=="fgn"){
      acf.vec=c( 1,0.5*( (kk+1)^(2*h[i])-2*kk^(2*h[i])+(kk-1)^(2*h[i]) ) )
    }else if (model=="arfima"){
      acorr = lgamma(1.5-h[i])+lgamma(kk+h[i]-0.5) - lgamma(h[i]-0.5)-lgamma(1.5-h[i]+kk) 
      acf.vec = c(1,exp(acorr))
    }else{
      print("Model is not supported...")
    }
    
    
    X = svd(toeplitz(acf.vec))
    tol = sqrt(.Machine$double.eps)
    a.func.corr[i]=sum(log(X$d[X$d > max(tol * X$d[1L], 0.0)])) 
  }
  distance=sqrt(-a.func.corr)/sqrt(n)
  
  
  h.seq2=h[h>=0.5]
  distance2=distance[h>=0.5]
  norm.const=1
  spline.dist2=splinefun(h.seq2,distance2) # Distance d2 as a function of H
  lambda=-log(2*alpha)/spline.dist2(upper)
  f.h2=norm.const*dexp(distance2,lambda)*abs(spline.dist2(h.seq2,deriv=1))
  prior=f.h2
  
  if(!persistent){
    h.seq1=h[h<0.5]
    distance1=distance[h<0.5]
    spline.dist1=splinefun(h.seq1,distance1) # Distance d1 as a function of H
    f.h1=norm.const*dexp(distance1,lambda)*abs(spline.dist1(h.seq1,deriv=1))
    prior = c(f.h1,f.h2)
  }
  if(persistent){
    fn = splinefun(rho.seq, log(prior)+log(0.5)-rho.seq-2*log(1+exp(-rho.seq)))
  }else{
    fn = splinefun(rho.seq, log(prior)-rho.seq-2*log(1+exp(-rho.seq)))
  }
  
  return (fn)
}