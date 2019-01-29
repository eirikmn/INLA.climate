h.map.maker = function(m=4,lagmax=1000,stoc="fgn"){
  
  if(!(stoc %in% c("fgn","arfima"))){
    return(-1)
  }
  
  filename = paste(stoc,"Lag",lagmax,"m",m,".txt",sep="") 
  if(!( filename %in% list.files(path=file.path(.Library,"/INLA.climate/Hmapping/",fsep="")))){
    return(-2)
  }
  # if(stoc=="fgn"){
  #   if( ( lagmax==1000 && !(m%in%c(3,4,5,6)) ) || ( lagmax==4296 &&!(m%in%c(3,4,5)) ) ){
  #     return(-2)
  #   }
  # }else if(stoc=="arfima"){
  #   if(  lagmax==1000 && !(m%in%c(4))){
  #     return(-2)
  #   }
  # }else{
  #   return(-1)
  # }
  
  filepath = file.path(.Library,"/INLA.climate/Hmapping/",filename,fsep="")
  
  data = read.table(filepath)
  
  funks = numeric(0)
  for(i in 1:(2*m)){
    fun = splinefun(data[,1],data[,i+1],method="natural")
    funks = c(funks,fun)
  }
  return(funks)
}





