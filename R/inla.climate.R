
inla.climate = function(data, forcing, Qco2=NULL,compute.mu=NULL, stepLength=0.01,restart.inla=FALSE,
                        m = 4,model="fgn", print.progress=FALSE,
                        inla.options = list(),
                        tcr.options = list(),
                        mu.options = list(),
                        ar1.options = list()){

  catch = tryCatch(attachNamespace("INLA"),error=function(x){})
  if(length(find.package("INLA",quiet=TRUE))==0){
      stop("This function requires INLA. Please install at www.R-INLA.org or by calling 'install.packages(\"INLA\", repos=c(getOption(\"repos\"), INLA=\"https://inla.r-inla-download.org/R/testing\"), dep=TRUE)' from R.")
  }
  
  
  if(print.progress){
    print("Initiating inla.climate..")
  }
  tid.start = proc.time()[[3]]

  inla.climate.call = sys.call(which=1)
  if(sum(is.null(forcing))>0) stop("Forcing contains NA values")

  old.digits=getOption("digits")
  if(old.digits < 7){
    options(digits = 7)
  }
  
  
  default.inla.options = list(
    num.threads=1,
    control.compute = list(cpo=FALSE,dic=(sum(is.na(data))==0),config=TRUE),
    control.inla = list(reordering="metis",h=stepLength[1],restart=restart.inla),
    control.family = list(hyper = list(prec = list(initial = 12, fixed=TRUE))) )
  
  inla.options=set.options(inla.options,default.inla.options)
  
  
  if(length(Qco2)>0){
    default.tcr.options = list(nsamples = 100000, seed = 1234)
    # temp = default.tcr.options
    
    tcr.options=set.options(tcr.options,default.tcr.options)
    
  }
  if(length(compute.mu)==0){
    compute.mu="No"
  }
  if(compute.mu %in% c(1,2,"full","complete","quick","fast")){
    default.mu.options = list(nsamples = 100000, seed = 1234)
    # temp = default.mu.options
    
    mu.options=set.options(mu.options,default.mu.options)
    
  }
  
  if(model == "ar1"){
    default.ar1.options = list(nsamples = 100000, seed = 1234)
    # temp = default.mu.options
    
    ar1.options=set.options(ar1.options,default.ar1.options)
    
  }

  
  if(class(data)=="numeric" || class(data)=="ts"){
    n = length(data)
    y0 = data
  }else if (class(data)=="data.frame"){
    n = length(data[,1])
    
    y0 = data[,1]
    
  }else{
    stop("'data' input not recognized. Only objects of class 'numeric', 'ts' and 'data.frame' are valid.")
  }
  
  if(class(forcing)=="numeric" || class(forcing)=="ts"){
    if(length(forcing)>n){
      y0 = c(y0,rep(NA,length(forcing)-n))
      n = length(forcing)
      
    }
    T0 = mean(data[1:20][!is.na(data[1:20])])
    if(is.null(T0)){
      stop("There are no non-NA values among the first 20 of data, please remove these before repeating.")
    }
    y=y0-T0
    df=data.frame(y=y,idy=1:n)
  }else{
    stop("'forcing' must be a 'numeric' object.")
  }

  if(sum(is.na(df[1:n,1]))>0){
    inla.options$control.compute$dic = FALSE
  }

  
  is.lrd = model %in% c("arfima","fgn")
  
  if(is.lrd){
    lagmax = 1000L
    funks = h.map.maker(m,lagmax,model)
    lprior.fun.H = compute.Hprior(50,0.9,0.1,persistent=TRUE,model=model)
    model.approx = INLA::inla.rgeneric.define(rgeneric.lrd,lprior.fun.H = lprior.fun.H,
                                              n=n,N=m,forcing=forcing,funks=funks)
  }else if(model == "ar1"){
    #m = 1 #only one component is available so far
    model.approx = INLA::inla.rgeneric.define(rgeneric.ar1,n=n,N=m,forcing=forcing)
    #model.approx = INLA::inla.rgeneric.define(rgeneric.forcing.3AR1.free,n=n,N=m,forcing=forcing)
  }
  
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
  
  if(model == "ar1"){
    ar1.result = inla.climate.ar1(result.approx, m=m, nsamples=ar1.options$nsamples,
                                  seed=ar1.options$seed, print.progress=print.progress)
  }
  
  if(!is.null(Qco2)){
    tcr.result = inla.climate.tcr(result.approx,Qco2,nsamples=tcr.options$nsamples,
                                  seed=tcr.options$seed, print.progress=print.progress,model=model)
  }
  if(compute.mu %in% c(1,2,"full","complete","quick","fast") ){
    if(compute.mu %in% c(2,"full","complete")){
      mu.quick = FALSE
    }else{
      mu.quick=TRUE
    }
    mu.result = inla.climate.mu(result.approx, forcing, quick=mu.quick, T0.corr = T0, nsamples=mu.options$nsamples,
                                seed=mu.options$seed, print.progress=print.progress,model=model)
  }
  
  
  if(print.progress){
    cat("Finishing up...\n",sep="")
  }
  
  misc = list(call=inla.climate.call, INLA.options=inla.options, TCR.options=tcr.options, mu.options=mu.options,
              data = data.frame(y=(df$y+T0),idy=(df$idy)), forcing=forcing, n=n, m=m, model=model, T0=T0,
              stepLength=stepLength, restart.inla=restart.inla, Qco2=Qco2, compute.mu=compute.mu,
              time.inla = tid.approx)
  
  results = process.inla(result.approx, misc)
  
  
  if(length(Qco2) > 0){
    results = process.tcr(results,tcr.result)
  }

  if(compute.mu %in% c(1,2,"full","complete","quick","fast")){
    results = process.mu(results,mu.result)
  }
  if(model == "ar1"){
    misc$ar1.options = ar1.options
    results = process.ar1(results,ar1.result)
  }
  
  if(print.progress){
    cat("Finishing up..\n",sep="")
  }
  
  tid.slutt = proc.time()[[3]]
  results$time$Total = tid.slutt-tid.start


  class(results) <- "inla.climate"
  
  options(digits=old.digits)
  
  return(results)
}

