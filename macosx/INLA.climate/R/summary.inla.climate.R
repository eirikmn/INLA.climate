summary.inla.climate = function(object,digits=4L,...){
  is.tcr=!is.null(object$TCR)
  #format(, nsmall = 4)
  #options(digits=6)
  timeString = paste("Time used:\n"," Computing INLA",sep="")
  if(is.tcr){
    timeString=paste(timeString,"Sampling TCR",sep="\t")
    timeString=paste(timeString,"\tOther\t","Total\n",sep="")
    timeString=paste(timeString," \t ",round(object$time$inla,digits=digits),"\t       ",sep="")
  }else{
    timeString=paste(timeString,"\tOther\t","Total\n",sep="")
    timeString=paste(timeString," \t ",round(object$time$inla,digits=digits),"\t ",sep="")
  }
  
  
  tid=object$time$inla
  if(is.tcr){
    timeString=paste(timeString,round(object$time$TCR,digits=digits),"\t",sep="")
    tid=tid+object$time$TCR
  }
  timeString = paste(timeString,round(object$time$Total-tid,digits=digits),"\t",
                     round(object$time$Total,digits=digits),sep="")
  
  
  hyperString=paste("\n\nModel hyperparameters:\n","\t","mean\t","sd\t","0.025q\t",
                    "0.5q\t","0.975q\n",sep="")
  hyperString=paste(hyperString,"H\t",round(object$hyperparam$means$H,digits=digits),"\t",
                    round(object$hyperparam$sd$H,digits=digits),"\t",
                    round(object$hyperparam$quant0.025$H,digits=digits),"\t",
                    round(object$hyperparam$quant0.5$H,digits=digits),"\t",
                    round(object$hyperparam$quant0.975$H,digits=digits),
                    "\n",sep="")
  hyperString=paste(hyperString,"Sigmax\t",round(object$hyperparam$means$sigmax,digits=digits),"\t",
                    round(object$hyperparam$sd$sigmax,digits=digits),"\t",
                    round(object$hyperparam$quant0.025$sigmax,digits=digits),"\t",
                    round(object$hyperparam$quant0.5$sigmax,digits=digits),"\t",
                    round(object$hyperparam$quant0.975$sigmax,digits=digits),
                    "\n",sep="")
  hyperString=paste(hyperString,"Sigmaf\t",round(object$hyperparam$means$sigmaf,digits=digits),"\t",
                    round(object$hyperparam$sd$sigmaf,digits=digits),"\t",
                    round(object$hyperparam$quant0.025$sigmaf,digits=digits),"\t",
                    round(object$hyperparam$quant0.5$sigmaf,digits=digits),"\t",
                    round(object$hyperparam$quant0.975$sigmaf,digits=digits),
                    "\n",sep="")
  hyperString=paste(hyperString,"Shift\t",round(object$hyperparam$means$shift,digits=digits),"\t",
                    round(object$hyperparam$sd$shift,digits=digits),"\t",
                    round(object$hyperparam$quant0.025$shift,digits=digits),"\t",
                    round(object$hyperparam$quant0.5$shift,digits=digits),"\t",
                    round(object$hyperparam$quant0.975$shift,digits=digits),
                    "\n",sep="")
  #cat(hyperString)
  TCRstring=""
  if(is.tcr){
    
    TCRstring=paste(TCRstring,"\nTransient climate response (",format(object$misc$mcmcsamples,scientific=F,digits=digits)," samples):",
                    "\n\tmean\tsd\t0.025q\t0.5q\t0.975q\n",sep="")
    
    TCRstring=paste(TCRstring,"TCR\t",round(object$TCR$mean,digits=digits),"\t",
                      round(object$TCR$sd,digits=digits),"\t",
                      round(object$TCR$quant0.025,digits=digits),"\t",
                      round(object$TCR$quant0.5,digits=digits),"\t",
                      round(object$TCR$quant0.975,digits=digits),
                      "\n",sep="")
    #cat(TCRstring)
  }
  
  approxstring = paste("\n",object$mis$stoc," approximated using ",object$misc$m," AR(1) processes",sep="")
  
  likelihoodstring=paste("\nlog marginal log-likelihood:  ",round(object$inla.result$mlik[1],digits=digits),"\n",sep="")
  dicstring1=paste("\nDeviance Information Criterion (DIC):     ",round(object$dic$dic,digits=digits),sep="")
  dicstring2=paste("\nDeviance Information Criterion (saturated): ",round(object$dic$dic.sat,digits=digits),"\n",sep="")
  utString=paste(timeString,hyperString,TCRstring,likelihoodstring,dicstring1,dicstring2,sep="")
  cat("Call:\n")
  show(object$misc$call)
  cat("\n",utString)
}
