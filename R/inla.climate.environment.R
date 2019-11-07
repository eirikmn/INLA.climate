process.inla = function(object, misc=NULL){
  if(is.null(object$climate.misc) && !is.null(misc)){
    object$climate.misc=misc
  }else if(is.null(object$climate.misc) && is.null(misc)){
    stop("Could not find inla.climate information.")
  }
  #misc: INLA.options, call, data, forcing, m, model, t0, stepLength, restart, time
  if(misc$model == "fgn"){
    var.name = "H"
    var.func = function(x)0.5+0.5/(1+exp(-x))
  }else if(misc$model == "arfima"){
    var.name = "d"
    var.func = function(x)0.5/(1+exp(-x))
  }else if(misc$model == "ar1"){
    var.name = "rho"
    var.func = function(x)2/(1+exp(-x))-1
  }
  
  
  margs = object$marginals.hyperpar
  a=3
  if(misc$model != "ar1"){
    
    mem.approx = INLA::inla.emarginal(var.func,margs$`Theta2 for idy`)
    sigmax.approx = INLA::inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta1 for idy`)
    sigmaf.approx = INLA::inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta3 for idy`)
    F0.approx = INLA::inla.emarginal(function(x) x,margs$`Theta4 for idy`)
    #F0.approx = INLA::inla.emarginal(function(x) -a + 2*a/(1+exp(-x)),margs$`Theta4 for idy`)
    
    margs.approx = object$marginals.hyperpar
    marg.mem = INLA::inla.tmarginal(var.func,margs.approx$`Theta2 for idy`)
    marg.sx = INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta1 for idy`)
    marg.sf = INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta3 for idy`)
    marg.F0 = INLA::inla.tmarginal(function(x) x,margs.approx$`Theta4 for idy`)
    #marg.F0 = INLA::inla.tmarginal(function(x) -a + 2*a/(1+exp(-x)),margs.approx$`Theta4 for idy`)
    
    zmarg.mem = INLA::inla.zmarginal(marg.mem,silent=TRUE)
    hpd.mem = INLA::inla.hpdmarginal(0.95,marg.mem)
  }else{
    sigmax.approx = INLA::inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta1 for idy`)
    sigmaf.approx = INLA::inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta2 for idy`)
    F0.approx = INLA::inla.emarginal(function(x) x,margs$`Theta3 for idy`)
    #F0.approx = INLA::inla.emarginal(function(x) -a + 2*a/(1+exp(-x)),margs$`Theta3 for idy`)
    
    margs.approx = object$marginals.hyperpar
    marg.sx = INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta1 for idy`)
    marg.sf = INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta2 for idy`)
    marg.F0 = INLA::inla.tmarginal(function(x) x,margs.approx$`Theta3 for idy`)
    #marg.F0 = INLA::inla.tmarginal(function(x) -a + 2*a/(1+exp(-x)),margs.approx$`Theta3 for idy`)
  }
  
  
  
    zmarg.sx = INLA::inla.zmarginal(marg.sx,silent=TRUE)
    zmarg.sf = INLA::inla.zmarginal(marg.sf,silent=TRUE)
    zmarg.F0 = INLA::inla.zmarginal(marg.F0,silent=TRUE)
    

    
    hpd.sx = INLA::inla.hpdmarginal(0.95,marg.sx)
    hpd.sf = INLA::inla.hpdmarginal(0.95,marg.sf)
    hpd.F0 = INLA::inla.hpdmarginal(0.95,marg.F0)
    
    #cat("Creds: ","(",hpd.H[1],",",hpd.H[2],") & (",hpd.sx[1],",",hpd.sx[2],") & (",
    #    hpd.sf[1],",",hpd.sf[2],") & (",hpd.F0[1],",",hpd.F0[2],") ",sep="")
  n=object$climate.misc$n
  T0 = object$climate.misc$T0
  results = list(inla.result=object, hyperparam = list(
      means=list(),
      sd = list(),
      quant0.025=list(),
      quant0.5=list(),
      quant0.975=list(),
      marginals=list()  ),
    latent.field = list(model.fit=list(means=object$summary.random$idy$mean[1:n]+T0,
                   sd = object$summary.random$idy$sd[1:n],
                   quant0.025=object$summary.random$idy$`0.025quant`[1:n]+T0,
                   quant0.5=object$summary.random$idy$`0.5quant`[1:n]+T0,
                   quant0.975=object$summary.random$idy$`0.975quant`[1:n]+T0) ),
    time=list(inla=object$climate.misc$time.inla) )

  if(misc$model != "ar1"){
    results$hyperparam$means[[var.name]] = mem.approx
    results$hyperparam$sd[[var.name]] = zmarg.mem$sd
    results$hyperparam$quant0.025[[var.name]] = zmarg.mem$quant0.025
    results$hyperparam$quant0.5[[var.name]] = zmarg.mem$quant0.5
    results$hyperparam$quant0.975[[var.name]] = zmarg.mem$quant0.975
    results$hyperparam$marginals[[var.name]] = marg.mem
  }

results$hyperparam$means[["sigmax"]] = sigmax.approx
results$hyperparam$means[["sigmaf"]] = sigmaf.approx
results$hyperparam$means[["F0"]] = F0.approx


results$hyperparam$sd[["sigmax"]] = zmarg.sx$sd
results$hyperparam$sd[["sigmaf"]] = zmarg.sf$sd
results$hyperparam$sd[["F0"]] = zmarg.F0$sd


results$hyperparam$quant0.025[["sigmax"]] = zmarg.sx$quant0.025
results$hyperparam$quant0.025[["sigmaf"]] = zmarg.sf$quant0.025
results$hyperparam$quant0.025[["F0"]] = zmarg.F0$quant0.025


results$hyperparam$quant0.5[["sigmax"]] = zmarg.sx$quant0.5
results$hyperparam$quant0.5[["sigmaf"]] = zmarg.sf$quant0.5
results$hyperparam$quant0.5[["F0"]] = zmarg.F0$quant0.5


results$hyperparam$quant0.975[["sigmax"]] = zmarg.sx$quant0.975
results$hyperparam$quant0.975[["sigmaf"]] = zmarg.sf$quant0.975
results$hyperparam$quant0.975[["F0"]] = zmarg.F0$quant0.975


results$hyperparam$marginals[["sigmax"]] = marg.sx
results$hyperparam$marginals[["sigmaf"]] = marg.sf
results$hyperparam$marginals[["F0"]] = marg.F0


if(misc$model %in% c("arfima","fgn","ar1")){ 
  for(comp in 1:object$climate.misc$m){
    results$latent.field[[paste0("AR.component.",comp)]] <- list(means=object$summary.random$idy$mean[1:n+n*(comp)],
                                                                 sd=object$summary.random$idy$sd[1:n+n*(comp)],
                                                                 quant0.025=object$summary.random$idy$`0.025quant`[1:n+n*(comp)],
                                                                 quant0.5=object$summary.random$idy$`0.5quant`[1:n+n*(comp)],
                                                                 quant0.975=object$summary.random$idy$`0.975quant`[1:n+n*(comp)])
  }
}  
  
  results$misc$INLA.options = object$climate.misc$INLA.options
  results$misc$TCR.options  = object$climate.misc$TCR.options
  results$misc$TCR.options$Qco2 = object$climate.misc$Qco2
  results$misc$mu.options   = object$climate.misc$mu.options
  results$misc$mu.options$compute.mu=object$climate.misc$compute.mu
  
  results$misc$call = object$climate.misc$call
  results$misc$data = object$climate.misc$data
  results$misc$forcing = object$climate.misc$forcing
  results$misc$m = object$climate.misc$m
  
  results$misc$model = object$climate.misc$model
  results$misc$T0 = object$climate.misc$T0
  
  results$misc$initialtheta = object$climate.misc$inla.options$initialtheta
  results$misc$stepLength = object$climate.misc$inla.options$stepLength
  results$misc$restart = object$climate.misc$inla.options$restart
  
  
  
  
  
  results$log.mlikelihood = object$climate.misc$mlik
  if(object$climate.misc$INLA.options$control.compute$dic && sum(is.na(object$climate.misc$data$y))==0){
    results$dic = results$inla.result$dic
  }
  
  
  
  class(results) <- "inla.climate"
  return(results)
}

process.tcr = function(object, tcr.result, misc=NULL){
  #misc 
  
  if(class(object)=="inla"){
    if(is.null(object$climate.misc) && !is.null(misc)){
      object$climate.misc = misc
    }else if(is.null(object$climate.misc) && is.null(misc)){
      stop("Could not find inla.climate information.")
    }
  
      ret = process.inla(object,misc)
  }else if(class(object) == "inla.climate"){
      ret = object
  }else{
    stop("Invalid 'object' class.")
  }
  # if("H" %in% names(tcr.result$samples)){
  #   is.lrd = TRUE
    # ret$TCR=list(mean=tcr.result$mean, sd=tcr.result$sd,
    #              quant0.025=tcr.result$quant0.025,
    #              quant0.5=tcr.result$quant0.5,
    #              quant0.975=tcr.result$quant0.975,
    #              #creds=tcr.result$creds,
    #              samples=list(
    #                TCR=tcr.result$samples$TCR, H=tcr.result$samples$H,
    #                sigmaf=tcr.result$samples$sigmaf,F0=tcr.result$samples$F0))
  
  ret$TCR = tcr.result
  ret$time$TCR = tcr.result$time
  # }else{
  #   is.lrd = FALSE
  #   m = (length(names(tcr.result$samples))-2)/2
  #   if(m==1){
  #     ret $TCR= list(mean=tcr.result$mean,sd = tcr.result$sd,
  #                quant0.025=tcr.result$quant0.025,
  #                quant0.5=tcr.result$quant0.5,
  #                quant0.975=tcr.result$quant0.975,
  #                samples=list(
  #                  tcr.result$samples
  #                  )
  #   }else{
  #     ret$TCR = list(mean=mean(hyperpars[,tcr.col]),sd = sd(hyperpars[,tcr.col]),
  #                quant0.025=INLA::inla.qmarginal(0.025,mcfit),
  #                quant0.5=INLA::inla.qmarginal(0.5,mcfit),
  #                quant0.975=INLA::inla.qmarginal(0.975,mcfit),
  #                samples=list(
  #                  TCR=hyperpars[,tcr.col],sigmaf=hyperpars[,1],shift=hyperpars[,2]))
  #     for(k in 1:m){
  #       ret$samples[[paste0("w",k)]] = ww[,k]
  #     }
  #     for(k in 1:m){
  #       ret$samples[[paste0("p",k)]] = LL[,k]+1
  #     }
  #   }
  # }
    
    # ret$misc$TCR.options$nsamples = object$climate.misc$tcr.options$nsamples
    # ret$misc$TCR.options$seed = object$climate.misc$tcr.options$seed
    
  return(ret)
}


process.ar1 = function(object, ar1.result, misc=NULL){
  #misc 
  
  if(class(object)=="inla"){
    if(is.null(object$climate.misc) && !is.null(misc)){
      object$climate.misc = misc
    }else if(is.null(object$climate.misc) && is.null(misc)){
      stop("Could not find inla.climate information.")
    }
    
    ret = process.inla(object,misc)
  }else if(class(object) == "inla.climate"){
    ret = object
  }else{
    stop("Invalid 'object' class.")
  }
  if(object$misc$m==1){
    marg = inla.tmarginal(function(x)1/(1+exp(-x)),object$inla.result$marginals.hyperpar$`Theta4 for idy`)
    zmarg = inla.zmarginal(marg,silent=TRUE)
    ret$hyperparam$means$p = zmarg$mean
    ret$hyperparam$sd$p = zmarg$sd
    ret$hyperparam$quant0.025$p = zmarg$quant0.025
    ret$hyperparam$quant0.5$p = zmarg$quant0.5
    ret$hyperparam$quant0.975$p = zmarg$quant0.975
    ret$hyperparam$marginals$p = marg
  }else{
    for(k in 1:object$misc$m){
      ret$hyperparam$means[[paste0("w",k)]] = ar1.result[[paste0("w",k)]]$mean
      ret$hyperparam$sd[[paste0("w",k)]] = ar1.result[[paste0("w",k)]]$sd
      ret$hyperparam$quant0.025[[paste0("w",k)]] = ar1.result[[paste0("w",k)]]$quant0.025
      ret$hyperparam$quant0.5[[paste0("w",k)]] = ar1.result[[paste0("w",k)]]$quant0.5
      ret$hyperparam$quant0.975[[paste0("w",k)]] = ar1.result[[paste0("w",k)]]$quant0.975
      ret$hyperparam$marginals[[paste0("w",k)]] = cbind(density(ar1.result[[paste0("w",k)]]$samples)$x,density(ar1.result[[paste0("w",k)]]$samples)$y)
    }
    for(k in 1:object$misc$m){
      ret$hyperparam$means[[paste0("p",k)]] = ar1.result[[paste0("p",k)]]$mean
      ret$hyperparam$sd[[paste0("p",k)]] = ar1.result[[paste0("p",k)]]$sd
      ret$hyperparam$quant0.025[[paste0("p",k)]] = ar1.result[[paste0("p",k)]]$quant0.025
      ret$hyperparam$quant0.5[[paste0("p",k)]] = ar1.result[[paste0("p",k)]]$quant0.5
      ret$hyperparam$quant0.975[[paste0("p",k)]] = ar1.result[[paste0("p",k)]]$quant0.975
      ret$hyperparam$marginals[[paste0("p",k)]] = cbind(density(ar1.result[[paste0("p",k)]]$samples)$x,density(ar1.result[[paste0("p",k)]]$samples)$y)
    }
  }
  
  
  return(ret)
}

process.mu = function(object, mu.result, misc=NULL){
  
  if(class(object)=="inla"){
    if(is.null(object$climate.misc) && !is.null(misc)){
      object$climate.misc = misc
    }else if(is.null(object$climate.misc) && is.null(misc)){
      stop("Could not find inla.climate information.")
    }
    ret = process.inla(object,misc)
  }else if(class(object)=="inla.climate"){
    ret = object
    
  }else{
    stop("Invalid 'object' class.")
  }
    
    
    ret$mu = list(mean=mu.result$mean, sd=mu.result$sd)
    
    if(ret$misc$mu.options$compute.mu %in% c(2,"full","complete")){
      compute.mu = 2 ###
      ret$mu$quant0.025=mu.result$quant0.025
      ret$mu$quant0.5=mu.result$quant0.5
      ret$mu$quant0.975=mu.result$quant0.975
      ret$mu$samples=list(
        mu=mu.result$samples$mu, H=mu.result$samples$H,
        sigmaf=mu.result$samples$sigmaf, F0=mu.result$samples$F0
      )
    }else{
      compute.mu = 1 ###
    }
    
    
    ret$time$mu = mu.result$time
    #ret$misc$mu.options$compute.mu = ret$climate.misc$compute.mu
    #ret$misc$mu.options$nsamples=ret$climate.misc$mu.options$nsamples
    #ret$misc$mu.options$seed = ret$climate.misc$mu.options$seed
    
    return(ret)
}


set.options = function(opt, default.opt){
  temp = default.opt
  
  if(length(opt)>0){
    for(i in 1:length(opt)){
      if(names(opt)[i] %in% names(default.opt)){
        if(!is.list(opt[[i]])){
          temp[[ names(opt)[i] ]] <- opt[[i]]
        }else{
          for(j in 1:length(opt[[i]])){
            temp[[ names(opt)[i] ]][[names(opt[[i]])[j]]] <- opt[[i]][[j]]
          }
        }
      }else{
        temp[[names(opt)[i]]] <- opt[[i]]
      }
    }
  }
  return(temp)
}




