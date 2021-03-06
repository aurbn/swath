#
score <- function(data= NULL, data.file= NULL, measure.id= NULL, value.var= NULL, score.function= NULL, score.name= NULL, rep.id= NULL, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
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
  
  setkeyv(data, c(measure.id, rep.id))
  rm(data.file)
  
  if(is.function(score.function)){
    score.function <- match.fun(score.function)
  } else {
    score.function <- function(x){'Not a function!'}
  }
  
  data[, eval(score.name):= score.function(get(value.var)), by=c(measure.id, rep.id)]
 
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