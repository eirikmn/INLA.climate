print.inla.climate = function(x,digits=4L,...){
  options(digits=digits)
  cat("Call:\n")
  show(x$misc$call)
  cat("\n",sep="")
  
  if(!is.null(x$time$TCR)){
    cat("Time used:\n  INLA = ",x$time$inla,", TCR = ",x$time$TCR,", Total = ",x$time$Total,"\n",sep="")
  }else{
    cat("Time used:\n  INLA = ",x$time$inla,", Total = ",x$time$Total,"\n",sep="")
  }
  
  cat("\n",x$mis$stoc," approximated using ",x$misc$m," AR(1) processes\n",sep="")
  if(!is.null(x$time$TCR)){
    cat("\nTCR obtained by ",format(x$misc$mcmcsamples,scientific=F,digits=digits)," MCMC samples and CO2 coefficient ",x$misc$C,sep="")
  }
  
  cat("\nmarginal log-likelihood:  ",round(x$log.mlikelihood[1],digits=digits),"\n",sep="")
  cat("\nDeviance Information Criterion (DIC):     ",round(x$dic$dic,digits=digits),sep="")
  cat("\nDeviance Information Criterion (saturated): ",round(x$dic$dic.sat,digits=digits),"\n",sep="")
  
  return(invisible(x))
}