#
complete.measurements <- function(data= NULL, data.file= NULL, flag= TRUE, flag.name= 'complete', measure.id= NULL, rep.id= NULL, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
  # Load libraries
  require(data.table, quietly=TRUE)  # data table is very efficeint way to work with large tables. so we used it.
  
  # Check and create directories and files
  if(is.null(data)){
    if(!check.directory(subDir= data.file)==1){
      stop(paste0('Error: ', data.file, ' does not exist!' ))
    }
  }
  
  if(write.file){
    if(is.null(output.file)){
      stop('Error: set output file.')
    }
    if(!create.directory(subDir= substr(output.file, 1, regexpr('/', output.file, fixed= TRUE)))==0){
      stop(paste0('Error: failed to create ', substr(output.file, 1, regexpr('/', output.file, fixed= TRUE))))
    }
  }
  
  # Read data
  if(is.null(data)){
    data <- fread(input= data.file, sep= '\t', header= TRUE, showProgress= FALSE)
  } 
  setkeyv(data, c(measure.id, rep.id))
  rm(data.file)
  
  # Keep only fragments detected en each run
  nreps <- ulength(unlist(data[, get(rep.id)], recursive= TRUE, use.names= FALSE))
  data[, times_detected:= ulength(get(rep.id)), by= measure.id]
  
  if(!flag){
    data <- data[times_detected==nreps]
  } else if(!is.null(flag.name)){
    data[, eval(flag.name)] <- data[, times_detected]==nreps
  } else {
    stop('Error flag name needed')
  } 
  
  data[, times_detected:= NULL]
  #rm(list= c('fragments', 'fragment.stat', 'nruns'))
  
  if(write.file){
    write.table(x= data, file= output.file, quote= FALSE, sep= '\t', row.names= FALSE)
  }
  
  if(return.table){
    return(data)
  } else {
    if(write.file){
      return(output.file)
    } else {
      #return(0)
    }
  }
  
}