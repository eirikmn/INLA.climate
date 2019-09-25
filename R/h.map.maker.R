h.map.maker = function(m=4,lagmax=1000,model="fgn"){
  
  if(!(model %in% c("fgn","arfima"))){
    stop(paste0("'",model,"' is not supported."))
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
    
  
  filename = paste(model,"Lag",lagmax,"m",m,".txt",sep="") 
  
  if(!( filename %in% list.files(path=folder.path))){
    stop("Model specifications are not supported.")
  }
  filepath = file.path(folder.path,filename,fsep="")
  
  coofdata = read.table(filepath)
  
  warn.old = getOption("warn") ## grid spacing in coofdatais too narrow for splinefun to distinguish
  options(warn=-1)              ## them properly. It produces a warning which is here suppressed.
  funks = numeric(0)
  for(i in 1:(2*m)){
    fun = splinefun(coofdata[,1],coofdata[,i+1],method="natural")
    funks = c(funks,fun)
  }
  options(warn=warn.old)
  
  return(funks)
}





