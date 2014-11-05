#
produce.stat <- function(data= NULL, data.file= NULL, measure.id= NULL, value.var= NULL, stat.name= NULL, group.id= NULL, stat.function= NULL, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
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
  
  stat.name <- c(stat.name)
  value.var <- rep(value.var, length.out=length(stat.name))
  if(!is.list(stat.function)){
    stat.function <- list(stat.function)
  }
  
  if(length(stat.name)==length(stat.function)){
    for(i in 1:length(stat.name)){
      if(is.function(stat.function[[i]])){
        stat.function.i <- match.fun(stat.function[[i]])
      } else {
        stat.function.i <- function(x){'Not a function!'}
      }
      data[, eval(stat.name[i]):= stat.function.i(get(value.var[[i]])), by=c(measure.id, group.id)] # SLOW! Does grouping again for each function
    }
  }
  
    
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