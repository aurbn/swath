#
threshold.measurements <- function(data= NULL, data.file= NULL, flag= TRUE, flag.name= 'threshold', measure.id= NULL, value.var= NULL, threshold= 0, direction= NULL, include= TRUE, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
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
  setkeyv(data, measure.id)
  rm(data.file)
  
  # Remove fragments below threshold
  if(direction=='lower'&include==TRUE){
    remove <- unique(data[get(value.var)<=threshold, measure.id, with= FALSE])
  } else if(direction=='lower'&include==FALSE){
    remove <- unique(data[get(value.var)<threshold, measure.id, with= FALSE])
  } else if(direction=='higher'&include==TRUE){
    remove <- unique(data[get(value.var)>=threshold, measure.id, with= FALSE])
  } else if(direction=='higher'&include==FALSE){
    remove <- unique(data[get(value.var)>threshold, measure.id, with= FALSE])
  } else if(direction=='equal'){
    remove <- unique(data[get(value.var)==threshold, measure.id, with= FALSE])
  } else {
    stop('Wrong direction: lower, higher or equal')
  }
    
  rm(threshold)
  
  if(!flag){
    data <- data[!(get(measure.id) %in% unlist(remove, recursive= TRUE, use.names= FALSE))]
  } else if(!is.null(flag.name)){
    data[, eval(flag.name)] <- !(data[, get(measure.id)] %in% unlist(remove, recursive= TRUE, use.names= FALSE))
  } else {
    stop('Error flag name needed')
  }
  rm(remove)
  
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