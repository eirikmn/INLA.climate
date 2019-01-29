
inla.climate = function(data, forcing, Qco2=NULL,m = 4, stoc="fgn", print.progress=FALSE,
                        inla.options = list(compute.dic=TRUE, nruns=1, initialtheta=NULL, stepLength = 0.01, verbose = FALSE),
                        tcr.options = list(mcmcsamples = 100000, seed = 1234)){

  if(print.progress){
    print("Initiating inla.climate..\n")
    print(inla.options$stepLength)
  }
  tid.start = proc.time()[[3]]

  kall = sys.call(which=1)

  lagmax = 1000L

  funks = h.map.maker(m,lagmax,stoc)
  if(class(funks)!="list"){
    if(funks==-1){
      print(paste("Model: ",stoc," is not supported...",sep=""))
      return(FALSE)
    }else if(funks==-2){
      print("Model specifications are not supported...")
      return(FALSE)
    }
  }



  n=length(forcing)
  if(class(data) == "data.frame"){
    n = length(data[,1])
    t0 = mean(data[1:20,1])
    data[,1] = data[,1]-t0
    npreds=n-length(data[,1])
    if(npreds>0){
      data[n+1:npreds,1]=NA
    }
  }else{
    n = length(data)
    t0 = mean(data[1:20])
    data = data - t0
    data=data.frame(y=data,idy=1:n)
    npreds=n-length(forcing)
    if(npreds>0){
      data=c(data,rep(NA,npreds))
    }
    #
    if(npreds>0){
      data[n+1:npreds,1]=NA
    }
  }

  if(sum(is.na(data[1:n,1]))>0){
    inla.options$compute.dic = FALSE
  }

  #yy = y - t0

  lprior.fun.H = compute.Hprior(50,0.9,0.1,persistent=T,stoc=stoc)


  model.approx = inla.rgeneric.define(rgeneric.forcing.fast,lprior.fun.H = lprior.fun.H,
                                      n=n,N=m,forcing=forcing,funks=funks)

  formula = y ~ -1+ f(idy, model=model.approx)

  if(print.progress){
    print("Starting INLA..")
  }

  ini.theta = inla.options$initialtheta
  tid.approx.start = proc.time()[[3]]
  for(run in 1:inla.options$nruns){
    result.approx = inla(formula,family="gaussian", data=data,num.threads=1,
                         control.mode=list(theta=ini.theta),
                         verbose=inla.options$verbose,control.compute=list(cpo=F,dic=inla.options$compute.dic,
                                                                           config=T),
                         control.inla=list(reordering="metis",h=inla.options$stepLength),
                         control.family = list(hyper = list(prec = list(initial = 12, fixed=TRUE))),
                         control.predictor=list(compute=T),silent=1L)
    ini.theta = result.approx$mode$theta
  }

  tid.approx.slutt = proc.time()[[3]]
  tid.approx = tid.approx.slutt-tid.approx.start
  if(print.progress){
    cat("INLA completed in ",tid.approx," seconds\n",sep="")
    if(is.null(Qco2)){
      print("Finishing up..\n")
    }
  }



  margs = result.approx$marginals.hyperpar
  H.approx = inla.emarginal(function(x) 0.5+0.5/(1+exp(-x)),margs$`Theta2 for idy`)
  sigmax.approx = inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta1 for idy`)
  sigmaf.approx = inla.emarginal(function(x) 1/sqrt(exp(x)),margs$`Theta3 for idy`)
  shift.approx = inla.emarginal(function(x) x,margs$`Theta4 for idy`)


  #cat("rgeneric approx:\n","  H:      ",H.approx,"\n  sigmax:     ",sigmax.approx,"\n  sigmaf: ",
  #    sigmaf.approx,"\n  shift: ",shift.approx,"\n     T0: ",t0,"\n",sep="")

  run.creds = T
  if(run.creds){ #marginals som matrix
    margs.approx = result.approx$marginals.hyperpar
    marg.H = inla.tmarginal(function(x) 0.5+0.5/(1+exp(-x)),margs.approx$`Theta2 for idy`)
    marg.sx = inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta1 for idy`)
    marg.sf = inla.tmarginal(function(x) 1/sqrt(exp(x)),margs.approx$`Theta3 for idy`)
    marg.shift = margs.approx$`Theta4 for idy`
    zmarg.H = inla.zmarginal(marg.H,silent=T)
    zmarg.sx = inla.zmarginal(marg.sx,silent=T)
    zmarg.sf = inla.zmarginal(marg.sf,silent=T)
    zmarg.shift = inla.zmarginal(margs.approx$`Theta4 for idy`,silent=T)

    hpd.H = inla.hpdmarginal(0.95,marg.H)
    hpd.sx = inla.hpdmarginal(0.95,marg.sx)
    hpd.sf = inla.hpdmarginal(0.95,marg.sf)
    hpd.shift = inla.hpdmarginal(0.95,marg.shift)

    #cat("Creds: ","(",hpd.H[1],",",hpd.H[2],") & (",hpd.sx[1],",",hpd.sx[2],") & (",
    #    hpd.sf[1],",",hpd.sf[2],") & (",hpd.shift[1],",",hpd.shift[2],") ",sep="")
  }
  results = list(inla.result=result.approx,
                 hyperparam=list(
                   means=list(H=H.approx,sigmax=sigmax.approx,sigmaf=sigmaf.approx,shift=shift.approx),
                   sd = list(H=zmarg.H$sd,
                             sigmax=zmarg.sx$sd,
                             sigmaf=zmarg.sf$sd,
                             shift=zmarg.shift$sd  ),
                   quant0.025=list(H=zmarg.H$quant0.025,
                                    sigmax=zmarg.sx$quant0.025,
                                    sigmaf=zmarg.sf$quant0.025,
                                    shift=zmarg.shift$quant0.025  ),
                   quant0.5=list(H=zmarg.H$quant0.5,
                                  sigmax=zmarg.sx$quant0.5,
                                  sigmaf=zmarg.sf$quant0.5,
                                  shift=zmarg.shift$quant0.5  ),
                   quant0.975=list(H=zmarg.H$quant0.975,
                                    sigmax=zmarg.sx$quant0.975,
                                    sigmaf=zmarg.sf$quant0.975,
                                    shift=zmarg.shift$quant0.975  ),
                   marginals=list(H=marg.H, sigmax=marg.sx,sigmaf=marg.sf,shift=marg.shift) ),
                   latent=list(means=result.approx$summary.random$idy$mean[1:n],
                             sd = result.approx$summary.random$idy$sd[1:n],
                             quant0.025=result.approx$summary.random$idy$`0.025quant`[1:n],
                             quant0.5=result.approx$summary.random$idy$`0.5quant`[1:n],
                             quant0.975=result.approx$summary.random$idy$`0.975quant`[1:n]),
                     #hpd.95=list(H=hpd.H,sigmaf=hpd.sf,sigmax=hpd.sx,hpd.shift=hpd.shift),
                   time=list(inla=tid.approx))

  if(length(Qco2)>0){

    tcr.result = inla.climate.tcr(result.approx,Qco2,nsamples=tcr.options$mcmcsamples,
                                  seed=tcr.options$seed, print.progress=print.progress)
    if(print.progress){
      cat("Finishing up..\n",sep="")
    }

    results$TCR=list(mean=tcr.result$TCR.mean, sd=tcr.result$TCR.sd,
                    quant0.025=tcr.result$TCR.quant0.025,
                    quant0.5=tcr.result$TCR.quant0.5,
                    quant0.975=tcr.result$TCR.quant0.975,
                     #creds=tcr.result$creds,
                    samples=list(
                      TCR=tcr.result$samples$TCR, H=tcr.result$samples$H,
                      sigmaf=tcr.result$samples$sigmaf,shift=tcr.result$samples$shift))
                    results$time$TCR = tcr.result$time
  }

  results$misc$call = kall
  results$misc$Qco2 = Qco2
  results$misc$m = m
  results$misc$mcmcsamples = tcr.options$mcmcsamples
  results$misc$stoc = stoc
  results$misc$T0 = t0
  results$misc$seed = tcr.options$seed
  results$misc$initialtheta = inla.options$initialtheta
  results$misc$stepLength = inla.options$stepLength

  results$log.mlikelihood = results$inla.result$mlik
  if(!is.null(inla.options$compute.dic)){
    results$dic = results$inla.result$dic
  }

  tid.slutt = proc.time()[[3]]
  results$time$Total = tid.slutt-tid.start


  class(results) <- "inla.climate"

  return(results)
}

#source("compute.Hprior.R")
#source("h.map.maker.R")
#source("rgeneric.forcing.fast.R")
#source("inla.climate.TCR.R")
#source("plot.inla.climate.R")
#source("summary.inla.climate.R")
#r=inla.climate(y,z,3.6)
#plot(r)
#summary(r)
