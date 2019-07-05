process.inla = function(object, misc=NULL){
  if(is.null(object$climate.misc) && !is.null(misc)){
    object$climate.misc=misc
  }else if(is.null(object$climate.misc) && is.null(misc)){
    stop("Could not find inla.climate information.")
  }
  #misc: INLA.options, call, data, forcing, m, stoc, t0, stepLength, restart, time
  a=3
  margs = object$marginals.hyperpar
  H.approx = INLA::inla.emarginal(function(x) 0.5+0.5/(1+exp(-x)),margs$`Theta2 for idy`)
  sigmax.approx = INLA::inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta1 for idy`)
  sigmaf.approx = INLA::inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta3 for idy`)
  F0.approx = INLA::inla.emarginal(function(x) -a + 2*a/(1+exp(-x)),margs$`Theta4 for idy`)
  
  
  run.creds = TRUE
  if(run.creds){
    margs.approx = object$marginals.hyperpar
    marg.H = INLA::inla.tmarginal(function(x) 0.5+0.5/(1+exp(-x)),margs.approx$`Theta2 for idy`)
    marg.sx = INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta1 for idy`)
    marg.sf = INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta3 for idy`)
    marg.F0 = INLA::inla.tmarginal(function(x) -a + 2*a/(1+exp(-x)),margs.approx$`Theta4 for idy`)
    
    zmarg.H = INLA::inla.zmarginal(marg.H,silent=T)
    zmarg.sx = INLA::inla.zmarginal(marg.sx,silent=T)
    zmarg.sf = INLA::inla.zmarginal(marg.sf,silent=T)
    zmarg.F0 = INLA::inla.zmarginal(marg.F0,silent=T)
    
    hpd.H = INLA::inla.hpdmarginal(0.95,marg.H)
    hpd.sx = INLA::inla.hpdmarginal(0.95,marg.sx)
    hpd.sf = INLA::inla.hpdmarginal(0.95,marg.sf)
    hpd.F0 = INLA::inla.hpdmarginal(0.95,marg.F0)
    
    #cat("Creds: ","(",hpd.H[1],",",hpd.H[2],") & (",hpd.sx[1],",",hpd.sx[2],") & (",
    #    hpd.sf[1],",",hpd.sf[2],") & (",hpd.F0[1],",",hpd.F0[2],") ",sep="")
  }
  n=object$climate.misc$n
  T0 = object$climate.misc$T0
  results = list(inla.result=object,
                 hyperparam=list(
                   means=list(H=H.approx,sigmax=sigmax.approx,sigmaf=sigmaf.approx,F0=F0.approx),
                   sd = list(H=zmarg.H$sd,
                             sigmax=zmarg.sx$sd,
                             sigmaf=zmarg.sf$sd,
                             F0=zmarg.F0$sd  ),
                   quant0.025=list(H=zmarg.H$quant0.025,
                                   sigmax=zmarg.sx$quant0.025,
                                   sigmaf=zmarg.sf$quant0.025,
                                   F0=zmarg.F0$quant0.025  ),
                   quant0.5=list(H=zmarg.H$quant0.5,
                                 sigmax=zmarg.sx$quant0.5,
                                 sigmaf=zmarg.sf$quant0.5,
                                 F0=zmarg.F0$quant0.5  ),
                   quant0.975=list(H=zmarg.H$quant0.975,
                                   sigmax=zmarg.sx$quant0.975,
                                   sigmaf=zmarg.sf$quant0.975,
                                   F0=zmarg.F0$quant0.975  ),
                   marginals=list(H=marg.H, sigmax=marg.sx,sigmaf=marg.sf,F0=marg.F0) ),
                 latent.field=list(
                   model.fit=list(means=object$summary.random$idy$mean[1:n]+T0,
                                  sd = object$summary.random$idy$sd[1:n],
                                  quant0.025=object$summary.random$idy$`0.025quant`[1:n]+T0,
                                  quant0.5=object$summary.random$idy$`0.5quant`[1:n]+T0,
                                  quant0.975=object$summary.random$idy$`0.975quant`[1:n]+T0)
                 ),
                 
                 #hpd.95=list(H=hpd.H,sigmaf=hpd.sf,sigmax=hpd.sx,hpd.F0=hpd.F0),
                 time=list(inla=object$climate.misc$time.inla))
  for(comp in 1:m){
    results$latent.field[[paste0("AR.component.",comp)]] <- list(means=object$summary.random$idy$mean[1:n+n*(comp)],
                                                                 sd=object$summary.random$idy$sd[1:n+n*(comp)],
                                                                 quant0.025=object$summary.random$idy$`0.025quant`[1:n+n*(comp)],
                                                                 quant0.5=object$summary.random$idy$`0.5quant`[1:n+n*(comp)],
                                                                 quant0.975=object$summary.random$idy$`0.975quant`[1:n+n*(comp)])
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
  
  results$misc$stoc = object$climate.misc$stoc
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
    ret$TCR=list(mean=tcr.result$mean, sd=tcr.result$sd,
                     quant0.025=tcr.result$quant0.025,
                     quant0.5=tcr.result$quant0.5,
                     quant0.975=tcr.result$quant0.975,
                     #creds=tcr.result$creds,
                     samples=list(
                       TCR=tcr.result$samples$TCR, H=tcr.result$samples$H,
                       sigmaf=tcr.result$samples$sigmaf,F0=tcr.result$samples$F0))
    ret$time$TCR = tcr.result$time
    # ret$misc$TCR.options$nsamples = object$climate.misc$tcr.options$nsamples
    # ret$misc$TCR.options$seed = object$climate.misc$tcr.options$seed
    
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




