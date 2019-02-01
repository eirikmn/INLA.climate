h.map.maker = function(m=4,lagmax=1000,stoc="fgn"){
  
  if(!(stoc %in% c("fgn","arfima"))){
    stop(paste0("'",stoc,"' is not supported."))
  }
  
#  load(file=system.file("inst/Hmapping/",package="INLA.climate"))
  
  folder.path=c()
  for(i in 1:length(.libPaths())){
    if("INLA.climate" %in% list.files(.libPaths()[i])){
      folder.path=paste0(.libPaths()[i],"/INLA.climate/Hmapping/")
    }
  }
  if(length(folder.path)==0){
    stop("Could not find package directory, please make sure that INLA.climate is installed within one of the libraries displayed by '.libPaths()'.")
  }
  

  
  #if(length(list.files(path=file.path(.Library,"/INLA.climate/Hmapping/",fsep="")))==0)
  if(length(list.files(path=folder.path))==0){
    stop(paste0("Directory ",folder.path," is empty. Please make sure that INLA.climate package is properly installed."))
    
  }
    
  
  filename = paste(stoc,"Lag",lagmax,"m",m,".txt",sep="") 
  
  if(!( filename %in% list.files(path=folder.path))){
    stop("Model specifications are not supported.")
  }
  filepath = file.path(folder.path,filename,fsep="")
  
  data = read.table(filepath)
  
  funks = numeric(0)
  for(i in 1:(2*m)){
    fun = splinefun(data[,1],data[,i+1],method="natural")
    funks = c(funks,fun)
  }
  return(funks)
}





