#
min.measurements <- function(data= NULL, data.file= NULL, flag= FALSE, flag.name= 'minimum', min= 3, measure.id= NULL, group.id= NULL, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
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
  } else {
    data <- copy(data)
  }
  setkeyv(data, c(measure.id, group.id))
  rm(data.file)
  
  # Keep only precursors with more or equal then 'min' fragments
  data[, n_measurements:= ulength(get(measure.id)), by= group.id]
  keep <- unique(data[n_measurements>=min, get(group.id)])
  
  if(!flag){
    data <- copy(data)[get(group.id) %in% keep]
  } else if(!is.null(flag.name)){
    data[, eval(flag.name)] <- (data[, get(group.id)] %in% unlist(keep, recursive= TRUE, use.names= FALSE))
  } else {
    stop('Error flag name needed')
  }
  
  data[, n_measurements:= NULL]
  #rm(list= c('precursors', 'precursor.stat', 'min'))
    
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