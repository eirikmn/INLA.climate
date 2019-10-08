plot.inla.climate = function(x,
                             plot.random = TRUE,
                             plot.ar = TRUE,
                             plot.hyperpars = TRUE,
                             plot.tcr.samples = TRUE,
                             plot.mu = TRUE,
                             postscript = FALSE,
                             pdf = FALSE,
                             prefix= "inla.climate.plots/figure-",
                             ...){
  n = x$inla.result$misc$configs$contents$length[1]
  figure.count=1L
  
  if(x$misc$model == "fgn"){
    var.name = "H"
  }else if(x$misc$model == "arfima"){
    var.name = "d"
  }else if(x$misc$model == "ar1"){
    var.name = c()
  }
  
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
    plot(x$latent.field$model.fit$quant0.5,lwd=1.5,type="l",ylim=range(x$latent.field$model.fit$quant0.025,x$latent.field$model.fit$quant0.975),
         xlab="Posterior mean: 2.5%, 50%, 97.5%",ylab="", main="Model fit")
    lines(x$latent.field$model.fit$quant0.025[1:n],lty=2)
    lines(x$latent.field$model.fit$quant0.975[1:n],lty=2)
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }
  if(plot.ar && x$misc$m > 1){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    skip.this = FALSE
    if(x$misc$m==1){
      par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    }else if(x$misc$m==2){
      par(mfrow=c(2,1),mar=(c(4,3,1,1.5)))
    }else if(x$misc$m==3){
      par(mfrow=c(2,2),mar=(c(3.5,2.5,0.8,0.8)))
    }else if(x$misc$m==4){
      par(mfrow=c(2,2),mar=(c(3.5,2.5,0.8,0.8)))
    }else if(x$misc$m==5){
      par(mfrow=c(3,2),mar=(c(3.5,2.5,0.8,0.8)))
    }else if(x$misc$m==6){
      par(mfrow=c(3,2),mar=(c(3.5,2.5,0.8,0.8)))
    }else{
      warning(paste0("The built-in plot function for AR(1) components is not compatible with m=",x$misc$m," and will be skipped."))
      skip.this = TRUE
    }
    if(!skip.this){
      rekke=numeric(0)
      for(iter in 1:x$misc$m){
        rekke=range(rekke,x$latent.field[[paste0("AR.component.",iter)]])
      }
      for(iter in 1:x$misc$m){
        plot(x$latent.field[[paste0("AR.component.",iter)]]$quant0.5,lwd=1.5,type="l",ylim=rekke,
             xlab="",ylab="", main="",cex.axis=0.8)
        lines(x$latent.field[[paste0("AR.component.",iter)]]$quant0.025[1:n],lwd=0.7,lty=1,col="gray")
        lines(x$latent.field[[paste0("AR.component.",iter)]]$quant0.975[1:n],lwd=0.7,lty=1,col="gray")
        title(xlab=paste0("AR(1) component #",iter),ylab="",line=2.,cex.lab=0.8)
      }
      
      if(postscript || pdf){
        if (names(dev.cur()) != "null device") {
          dev.off()
        }
      }
    }
    
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
  }
  if(plot.hyperpars){
    if(x$misc$model != "ar1"){
      figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
      par(mfrow=c(2,2),mar=(c(3.8,2.6,0.8,1)))
      plot(x$hyperparam$marginals[[var.name]],xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means[[var.name]],lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025[[var.name]],lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975[[var.name]],lwd=0.8,col="gray")
      plot(x$hyperparam$marginals$sigmax,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$sigmax,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$sigmax,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$sigmax,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",sigma[epsilon],")")),ylab="",line=2.5,cex.lab=1)
      
      plot(x$hyperparam$marginals$sigmaf,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$sigmaf,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$sigmaf,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$sigmaf,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",sigma[f],")")),ylab="",line=2.5,cex.lab=1)
      
      plot(x$hyperparam$marginals$F0,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$F0,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$F0,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$F0,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",F[0],")")),ylab="",line=2.5,cex.lab=1)
      
      
    }else if(x$misc$model == "ar1" && x$misc$m ==1){
      figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
      par(mfrow=c(2,2),mar=(c(3.8,2.6,0.8,1)))
      
      plot(x$hyperparam$marginals$p,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$p,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$p,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$p,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",phi,")")),ylab="",line=2.5,cex.lab=1)
      plot(x$hyperparam$marginals$sigmax,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$sigmax,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$sigmax,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$sigmax,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",sigma[epsilon],")")),ylab="",line=2.5,cex.lab=1)
      
      plot(x$hyperparam$marginals$sigmaf,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$sigmaf,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$sigmaf,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$sigmaf,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",sigma[f],")")),ylab="",line=2.5,cex.lab=1)
      
      plot(x$hyperparam$marginals$F0,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$F0,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$F0,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$F0,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",F[0],")")),ylab="",line=2.5,cex.lab=1)
    }else if(x$misc$model == "ar1" && x$misc$m <=3 && x$misc$m >1){
      figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
      par(mfrow=c(3,3),mar=(c(3.8,2.6,0.8,1)))
      
      for(k in 1:x$misc$m){
        plot(x$hyperparam$marginals[[paste0("w",k)]],xlab="",ylab="",
             main="",type="l")
        abline(v=x$hyperparam$means[[paste0("w",k)]],lwd=0.8,col="black")
        abline(v=x$hyperparam$quant0.025[[paste0("w",k)]],lwd=0.8,col="gray")
        abline(v=x$hyperparam$quant0.975[[paste0("w",k)]],lwd=0.8,col="gray")
        title(xlab=bquote(paste("Posterior density (",w[.(k)],")")),ylab="",line=2.5,cex.lab=1)
      }
      for(k in 1:x$misc$m){
        plot(x$hyperparam$marginals[[paste0("p",k)]],xlab="",ylab="",
             main="",type="l")
        abline(v=x$hyperparam$means[[paste0("p",k)]],lwd=0.8,col="black")
        abline(v=x$hyperparam$quant0.025[[paste0("p",k)]],lwd=0.8,col="gray")
        abline(v=x$hyperparam$quant0.975[[paste0("p",k)]],lwd=0.8,col="gray")
        title(xlab=bquote(paste("Posterior density (",phi[.(k)],")")),ylab="",line=2.5,cex.lab=1)
      }
      
      
      plot(x$hyperparam$marginals$sigmax,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$sigmax,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$sigmax,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$sigmax,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",sigma[epsilon],")")),ylab="",line=2.5,cex.lab=1)
      
      plot(x$hyperparam$marginals$sigmaf,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$sigmaf,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$sigmaf,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$sigmaf,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",sigma[f],")")),ylab="",line=2.5,cex.lab=1)
      
      plot(x$hyperparam$marginals$F0,xlab="",ylab="",
           main="",type="l")
      abline(v=x$hyperparam$means$F0,lwd=0.8,col="black")
      abline(v=x$hyperparam$quant0.025$F0,lwd=0.8,col="gray")
      abline(v=x$hyperparam$quant0.975$F0,lwd=0.8,col="gray")
      title(xlab=expression(paste("Posterior density (",F[0],")")),ylab="",line=2.5,cex.lab=1)
      
      
      
      par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
      if(postscript || pdf){
        if (names(dev.cur()) != "null device") {
          dev.off()
        }
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
    
    if(x$misc$mu.options$compute.mu %in% c(2,"full","complete")){
      plot(x$misc$data$y,pch=19,cex=0.8,col="gray",ylab="", main="Temperature response to forcing",
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

