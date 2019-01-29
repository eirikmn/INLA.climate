plot.inla.climate = function(x,
                             plot.latent = TRUE,
                             plot.hyperpars = TRUE,
                             plot.samples.tcr = TRUE,...){
  n = x$inla.result$misc$configs$contents$length[1]
  if(plot.latent){
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    plot(x$latent$quant0.5,lwd=1.5,type="l",ylim=range(x$latent$quant0.025,x$latent$quant0.975),
         xlab="Posterior Mean: 2.5%, 50%, 9.75%",ylab="")
    lines(x$latent$quant0.025[1:n],lty=2)
    lines(x$latent$quant0.975[1:n],lty=2)
  }
  if(plot.hyperpars){
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
    plot(x$hyperparam$marginals$shift,xlab="",ylab="",
         main=expression(paste("Posterior density [",F[0],"]")),type="l")
    abline(v=x$hyperparam$means$shift,lwd=0.8,col="black")
    abline(v=x$hyperparam$quant0.025$shift,lwd=0.8,col="gray")
    abline(v=x$hyperparam$quant0.975$shift,lwd=0.8,col="gray")
  }

  if(plot.samples.tcr && !is.na(x$TCR)){
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    hist(x$TCR$samples$TCR,col="orange",xlab="TCR",ylab="Density",main="",
         probability=T,breaks=min(80,round(length(x$TCR$samples$TCR)/250)))
  }

  return(invisible(x))
}

