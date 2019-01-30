
inla.climate = function(data, forcing, Qco2=NULL,compute.mu=FALSE, stepLength=0.01,restart.inla=FALSE, m = 4, stoc="fgn", print.progress=FALSE,
                        inla.options = list(),
                        tcr.options = list(),
                        mu.options = list() ){

  atch = tryCatch(attachNamespace("INLA"),error=function(x){})
  if(length(find.package("INLA",quiet=TRUE))==0){
      stop("This function requires INLA. Please install at www.R-INLA.org or by calling 'install.packages(\"INLA\", repos=c(getOption(\"repos\"), INLA=\"https://inla.r-inla-download.org/R/testing\"), dep=TRUE)' from R.")
  }
  
  
  if(print.progress){
    print("Initiating inla.climate..\n")
  }
  tid.start = proc.time()[[3]]

  kall = sys.call(which=1)
  if(sum(is.null(forcing))>0) stop("Forcing contains NA values")

  old.digits=getOption("digits")
  if(old.digits < 6){
    options(digits = 6)
  }
  
  #if(is.null(inla.options))
  default.inla.options = list(
    num.threads=1,
    control.compute = list(cpo=FALSE,dic=(sum(is.na(data))>0),config=TRUE),
    control.inla = list(reordering="metis",h=stepLength[1],restart=restart.inla),
    control.family = list(hyper = list(prec = list(initial = 12, fixed=TRUE))) )
  temp = default.inla.options
  
  if(length(inla.options)>0){
    for(i in 1:length(inla.options)){
      if(names(inla.options)[i] %in% names(default.inla.options)){
        if(!is.list(inla.options[[i]])){
          temp[[ names(inla.options)[i] ]] <- inla.options[[i]]
        }else{
          for(j in 1:length(inla.options[[i]])){
            temp[[ names(inla.options)[i] ]][[names(inla.options[[i]])[j]]] <- inla.options[[i]][[j]]
          }
        }
      }else{
        temp[[names(inla.options)[i]]] <- inla.options[[i]]
      }
    }
  }
  inla.options=temp
  
  if(length(Qco2)>0){
    default.tcr.options = list(mcsamples = 100000, seed = 1234)
    temp = default.tcr.options
    
    if(length(tcr.options)>0){
      for(i in 1:length(tcr.options)){
        if(names(tcr.options)[i] %in% names(default.tcr.options)){
          if(!is.list(tcr.options[[i]])){
            temp[[ names(tcr.options)[i] ]] <- tcr.options[[i]]
          }else{
            for(j in 1:length(tcr.options[[i]])){
              temp[[ names(tcr.options)[i] ]][[names(tcr.options[[i]])[j]]] <- tcr.options[[i]][[j]]
            }
          }
        }else{
          temp[[names(tcr.options)[i]]] <- tcr.options[[i]]
        }
      }
    }
    tcr.options = temp
  }
  
  if(compute.mu){
    default.mu.options = list(full.Bayesian = FALSE, mcsamples = 100000, seed = 1234)
    temp = default.mu.options
    
    if(length(mu.options)>0){
      for(i in 1:length(mu.options)){
        if(names(mu.options)[i] %in% names(default.mu.options)){
          if(!is.list(mu.options[[i]])){
            temp[[ names(mu.options)[i] ]] <- mu.options[[i]]
          }else{
            for(j in 1:length(mu.options[[i]])){
              temp[[ names(mu.options)[i] ]][[names(mu.options[[i]])[j]]] <- mu.options[[i]][[j]]
            }
          }
        }else{
          temp[[names(mu.options)[i]]] <- mu.options[[i]]
        }
      }
    }
    mu.options = temp
  }
  
  
  
  lagmax = 1000L

  funks = h.map.maker(m,lagmax,stoc)

  n=length(forcing)
  if(class(data)=="numeric" || class(data)=="ts"){
    if(length(data)>n){
      data=c(data,rep(NA,n-length(data)))
    }
    t0 = mean(data[1:20][!is.na(data[1:20])])
    if(is.null(t0)){
      stop("There are no non-NA values among the first 20 of data, please remove these before repeating.")
    }
    y=data-t0
    df=data.frame(y=y,idy=1:n)
  }else{
    stop("'data' must be a 'numeric' object.")
  }

  if(sum(is.na(df[1:n,1]))>0){
    inla.options$control.compute$dic = FALSE
  }


  lprior.fun.H = compute.Hprior(50,0.9,0.1,persistent=T,stoc=stoc)

  model.approx = INLA::inla.rgeneric.define(rgeneric.forcing.fast,lprior.fun.H = lprior.fun.H,
                                      n=n,N=m,forcing=forcing,funks=funks)
  formula = y ~ -1+ f(idy, model=model.approx)

  if(print.progress){
    print("Starting INLA..")
  }

  tid.approx.start = proc.time()[[3]]
  
  result.approx <- tryCatch(
    do.call(INLA::inla, c(list(formula = formula,data = df,family = "gaussian"),inla.options) ),
    error=warning
    )
  
  if(is.character(result.approx)){
    feil = "\n Convergence can sometimes be improved by changing the step length h."
    stop(paste0(result.approx,feil))
  }

  tid.approx.slutt = proc.time()[[3]]
  tid.approx = tid.approx.slutt-tid.approx.start
  if(print.progress){
    cat("INLA completed in ",tid.approx," seconds\n",sep="")
  }

  if(!is.null(Qco2)){
    tcr.result = inla.climate.tcr(result.approx,Qco2,nsamples=tcr.options$mcsamples,
                                  seed=tcr.options$seed, print.progress=print.progress)
  }
  if(compute.mu){
    mu.result = inla.climate.mu(result.approx, forcing, nsamples=mu.options$mcsamples,
                                full.Bayesian=mu.options$full.Bayesian,seed=mu.options$seed,
                                print.progress=print.progress)
  }
  if(print.progress){
    cat("Finishing up...\n",sep="")
  }
  margs = result.approx$marginals.hyperpar
  H.approx = INLA::inla.emarginal(function(x) 0.5+0.5/(1+exp(-x)),margs$`Theta2 for idy`)
  sigmax.approx = INLA::inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta1 for idy`)
  sigmaf.approx = INLA::inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta3 for idy`)
  F0.approx = INLA::inla.emarginal(function(x) x,margs$`Theta4 for idy`)


  run.creds = T
  if(run.creds){ #marginals som matrix
    margs.approx = result.approx$marginals.hyperpar
    marg.H = INLA::inla.tmarginal(function(x) 0.5+0.5/(1+exp(-x)),margs.approx$`Theta2 for idy`)
    marg.sx = INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta1 for idy`)
    marg.sf = INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta3 for idy`)
    marg.F0 = margs.approx$`Theta4 for idy`
    zmarg.H = INLA::inla.zmarginal(marg.H,silent=T)
    zmarg.sx = INLA::inla.zmarginal(marg.sx,silent=T)
    zmarg.sf = INLA::inla.zmarginal(marg.sf,silent=T)
    zmarg.F0 = INLA::inla.zmarginal(margs.approx$`Theta4 for idy`,silent=T)

    hpd.H = INLA::inla.hpdmarginal(0.95,marg.H)
    hpd.sx = INLA::inla.hpdmarginal(0.95,marg.sx)
    hpd.sf = INLA::inla.hpdmarginal(0.95,marg.sf)
    hpd.F0 = INLA::inla.hpdmarginal(0.95,marg.F0)

    #cat("Creds: ","(",hpd.H[1],",",hpd.H[2],") & (",hpd.sx[1],",",hpd.sx[2],") & (",
    #    hpd.sf[1],",",hpd.sf[2],") & (",hpd.F0[1],",",hpd.F0[2],") ",sep="")
  }
  results = list(inla.result=result.approx,
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
                   latent=list(means=result.approx$summary.random$idy$mean[1:n],
                             sd = result.approx$summary.random$idy$sd[1:n],
                             quant0.025=result.approx$summary.random$idy$`0.025quant`[1:n],
                             quant0.5=result.approx$summary.random$idy$`0.5quant`[1:n],
                             quant0.975=result.approx$summary.random$idy$`0.975quant`[1:n]),
                     #hpd.95=list(H=hpd.H,sigmaf=hpd.sf,sigmax=hpd.sx,hpd.F0=hpd.F0),
                   time=list(inla=tid.approx))

  results$misc$INLA.options = inla.options
  if(length(Qco2)>0){
    results$TCR=list(mean=tcr.result$TCR.mean, sd=tcr.result$TCR.sd,
                    quant0.025=tcr.result$TCR.quant0.025,
                    quant0.5=tcr.result$TCR.quant0.5,
                    quant0.975=tcr.result$TCR.quant0.975,
                     #creds=tcr.result$creds,
                    samples=list(
                      TCR=tcr.result$samples$TCR, H=tcr.result$samples$H,
                      sigmaf=tcr.result$samples$sigmaf,F0=tcr.result$samples$F0))
        results$time$TCR.option = tcr.result$time
        results$misc$TCR.option$mcsamples = tcr.options$mcsamples
        results$misc$TCR.option$seed = tcr.options$seed
  }

  if(compute.mu){
    results$mu = list(mean=mu.result$mu.mean, sd=mu.result$mu.sd)
    if(mu.options$full.Bayesian){
      results$mu$quant0.025=mu.result$mu.quant0.025
      results$mu$quant0.5=mu.result$mu.quant0.5
      results$mu$quant0.975=mu.result$mu.quant0.975
      results$mu$samples=list(
        mu=mu.result$samples$mu, H=mu.result$samples$H,
        sigmaf=mu.result$samples$sigmaf, F0=mu.result$samples$F0
      )
    }
    results$time$mu = mu.result$time
    results$misc$mu.option$full.Bayesian = mu.options$full.Bayesian
    results$misc$mu.option$mcsamples=mu.options$mcsamples
    results$misc$mu.option$seed = mu.options$seed
  }
  
  
  if(print.progress){
    cat("Finishing up..\n",sep="")
  }
  
  results$misc$call = kall
  results$misc$data = df
  results$misc$forcing = forcing
  results$misc$Qco2 = Qco2
  results$misc$m = m
  
  results$misc$stoc = stoc
  results$misc$T0 = t0
  
  results$misc$initialtheta = inla.options$initialtheta
  results$misc$stepLength = inla.options$stepLength
  results$misc$restart = inla.options$restart
  

  results$log.mlikelihood = results$inla.result$mlik
  if(inla.options$control.compute$dic && sum(is.na(y))==0){
    results$dic = results$inla.result$dic
  }

  tid.slutt = proc.time()[[3]]
  results$time$Total = tid.slutt-tid.start


  class(results) <- "inla.climate"
  
  options(digits=old.digits)
  
  return(results)
}

