plot.inla.climate = function(x,
                             plot.random = TRUE,
                             plot.hyperpars = TRUE,
                             plot.tcr.samples = TRUE,
                             plot.mu = TRUE,
                             postscript = FALSE,
                             pdf = FALSE,
                             prefix= "inla.climate.plots/figure-",
                             ...){
  n = x$inla.result$misc$configs$contents$length[1]
  figure.count=1L
  
  
  if(!postscript && !pdf){
    dev=getOption("device")
    if(!is.character(dev) || dev != "RStudioGD"){
      dev.new(...)
    }
  }else{
    dir = dirname(prefix)
    if (!file.exists(dir) && nchar(dir) > 0L) {
      dir.create(dir, recursive=TRUE)
    } else {
      stopifnot(file.info(dir)$isdir)
    }
  }
  
  
  if(plot.random){
    
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    plot(x$model.fit$quant0.5,lwd=1.5,type="l",ylim=range(x$model.fit$quant0.025,x$model.fit$quant0.975),
         xlab="Posterior mean: 2.5%, 50%, 97.5%",ylab="", main="Model fit")
    lines(x$model.fit$quant0.025[1:n],lty=2)
    lines(x$model.fit$quant0.975[1:n],lty=2)
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }
  if(plot.hyperpars){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(2,2),mar=(c(2.6,2.2,2.3,1.4)+0.1))
    plot(x$hyperparam$marginals$H,xlab="",ylab="",
         main=expression(paste("Posterior density [H]")),type="l")
    abline(v=x$hyperparam$means$H,lwd=0.8,col="black")
    abline(v=x$hyperparam$quant0.025$H,lwd=0.8,col="gray")
    abline(v=x$hyperparam$quant0.975$H,lwd=0.8,col="gray")
    plot(x$hyperparam$marginals$sigmax,xlab="",ylab="",
         main=expression(paste("Posterior density [",sigma[epsilon],"]")),type="l")
    abline(v=x$hyperparam$means$sigmax,lwd=0.8,col="black")
    abline(v=x$hyperparam$quant0.025$sigmax,lwd=0.8,col="gray")
    abline(v=x$hyperparam$quant0.975$sigmax,lwd=0.8,col="gray")
    plot(x$hyperparam$marginals$sigmaf,xlab="",ylab="",
         main=expression(paste("Posterior density [",sigma[f],"]")),type="l")
    abline(v=x$hyperparam$means$sigmaf,lwd=0.8,col="black")
    abline(v=x$hyperparam$quant0.025$sigmaf,lwd=0.8,col="gray")
    abline(v=x$hyperparam$quant0.975$sigmaf,lwd=0.8,col="gray")
    plot(x$hyperparam$marginals$F0,xlab="",ylab="",
         main=expression(paste("Posterior density [",F[0],"]")),type="l")
    abline(v=x$hyperparam$means$F0,lwd=0.8,col="black")
    abline(v=x$hyperparam$quant0.025$F0,lwd=0.8,col="gray")
    abline(v=x$hyperparam$quant0.975$F0,lwd=0.8,col="gray")
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
    
  }

  if(plot.tcr.samples && !is.null(x$TCR)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    hist(x$TCR$samples$TCR,col="orange",xlab="TCR",ylab="Density",main="",
         probability=T,breaks=min(80,round(length(x$TCR$samples$TCR)/250)))
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }
  if(plot.mu && !is.null(x$mu)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    
    if(x$misc$mu.option$compute.mu %in% c(2,"full","complete")){
      plot(x$misc$data$y,pch=19,cex=0.8,col="gray",ylab="", main="Forcing response vs data",
           xlab="Posterior mean: 2.5%, 50%, 97.5%",
           ylim=range(x$misc$data$y,x$mu$quant0.025,x$mu$quant0.975))
      lines(x$mu$quant0.5[1:n],lwd=1.5)
      lines(x$mu$quant0.025[1:n],lty=2)
      lines(x$mu$quant0.975[1:n],lty=2)
    }else{
      plot(x$misc$data$y,pch=19,cex=0.8,col="gray",ylab="",xlab="Posterior mean",
           ylim=range(x$misc$data$y,x$mu$mean),main="Forcing response vs data")
      lines(x$mu$mean[1:n],lwd=1.5)
    }
    
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }
  


  return(invisible(x))
}

new.plot = function(postscript,pdf,prefix,figure.count,...)
{

  #dev = getOption("device")
  if(postscript && pdf){
    stop("Multiple file types have been seleced.")
  }
  if(postscript) {
    ending=".eps"
  }else if(pdf){
    ending=".pdf"
  }
  if(!postscript && !pdf){
    dev=getOption("device")
    if(!is.character(dev) || dev != "RStudioGD"){
      dev.new(...)
    }
  }else{
    
    file.found=FALSE
    while(!file.found){
      filename=paste(prefix,figure.count,ending,sep="")
      
      if(file.exists(filename)){
        figure.count <- figure.count +1L
      }else{
        file.found=TRUE
      }
    }
    if(postscript){
      postscript(file=filename,...)
    }else if(pdf){
      pdf(file=filename,...)
    }
  }
  return (invisible(figure.count))
}

