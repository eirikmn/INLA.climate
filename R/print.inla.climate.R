print.inla.climate = function(x,digits=4L,...){
  
  cat("Call:\n")
  cat(deparse(x$misc$call),"\n\n",sep="")
  cat("Time used:\n",sep="")
  
  cpu = round(x$time$inla,digits=digits)
  cpu.navn="Running INLA"
  if(!is.null(x$TCR)){
    cpu=c(cpu,round(x$time$TCR,digits=digits))#fit model to observed temperatures that includes forcing contributions
    cpu.navn=c(cpu.navn,"Sampling TCR")
  }
  if(!is.null(x$mu)){
    cpu=c(cpu,round(x$time$mu,digits=digits))
    cpu.navn=c(cpu.navn,"Sampling mu")
  }
  cpu=c(cpu,round(x$time$Total,digits=digits))
  cpu.navn=c(cpu.navn,"Total")
  names(cpu)=cpu.navn
  print(cpu)
  
  if(!is.null(x$inla.result$summary.hyperpar) && length(x$inla.result$summary.hyperpar)>0L){
    cat(paste("\nModel contains ",x$inla.result$nhyper," hyperparameters:\n",sep=""))
    cat(names(x$hyperparam$means),"\n",sep=" ")
  }else{
    cat("\nThe model has no hyperparameters.\n",sep="")
  }
  if(!is.null(x$inla.result$summary.random) && length(x$inla.result$summary.random)>0L){
    cat("\nModel contains ",length(x$inla.result$summary.random)," random effects: \n",sep="")
    for(i in 1:length(x$inla.result$summary.random)){
      cat("'",names(x$inla.result$summary.random)[i],"': ", x$inla.result$model.random[[i]],"\n",sep="")
    }
  }else{
    cat("\nThe model has no random effects.\n",sep="")
  }
  
  
  cat("\nNoise model '",x$misc$model,"' was approximated using ",x$misc$m," AR(1) processes.\n",sep="")
  if(!is.null(x$time$TCR)){
    cat("\nTCR estimate obtained by ",format(x$misc$TCR.options$nsamples,scientific=F,digits=digits),
        " Monte Carlo samples with CO2 coefficient ",x$misc$TCR.options$Qco2,".\n",sep="")
  }
  if(!is.null(x$time$mu)){
    if(x$misc$mu.options$compute.mu %in% c(2,"full","complete")){
      cat("\nFull Bayesian analysis of forcing response computed from ",format(x$misc$mu.options$nsamples,scientific=F,digits=digits)," Monte Carlo samples.\n",sep="")
    }else{
      cat("\nQuick estimate of forcing response computed from ",format(x$misc$mu.options$nsamples,scientific=F,digits=digits)," Monte Carlo samples.\n",sep="")
    }
    
  }
  
  cat("\nLikelihood model: ",x$inla.result$.args$family,"\n",sep="")
  if(is.null(x$inla.result$.args$control.family$hyper$prec$fixed)){
    cat("log precision in the likelihood is fixed\n",sep="")
  }else{
    cat("log precision in the likelihood is random\n",sep="")
  }
  
  
  return(invisible(x))
}