h.map.maker = function(m=4,lagmax=1000,stoc="fgn"){
  
  if(!(stoc %in% c("fgn","arfima"))){
    stop(paste0("'",stoc,"' is not supported."))
  }
  
  if(length(list.files(path=file.path(.Library,"/INLA.climate/Hmapping/",fsep="")))==0)
    stop(paste0("Directory ",file.path(.Library,"/INLA.climate/Hmapping/",fsep="")," is empty. Please make sure that INLA.climate package is properly installed."))
  
  filename = paste(stoc,"Lag",lagmax,"m",m,".txt",sep="") 
  
  if(!( filename %in% list.files(path=file.path(.Library,"/INLA.climate/Hmapping/",fsep="")))){
    stop("Model specifications are not supported.")
  }
  filepath = file.path(.Library,"/INLA.climate/Hmapping/",filename,fsep="")
  
  data = read.table(filepath)
  
  funks = numeric(0)
  for(i in 1:(2*m)){
    fun = splinefun(data[,1],data[,i+1],method="natural")
    funks = c(funks,fun)
  }
  return(funks)
}





