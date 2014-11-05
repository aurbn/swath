#
normalize.ms <- function(data= NULL, measure.id= 'fragment_id', value.var= 'intensity', coef= NULL, by_group= TRUE, data.file= NULL, coef.file= NULL, output.file= NULL, return.table= TRUE, return.coef= FALSE, write.file= !is.null(output.file), output.coef= NULL, write.coef= !is.null(output.coef)){
  
  # Load libraries
  require(data.table, quietly=TRUE)  # data table is very efficeint way to work with large tables. so we used it.
  
  # Check and create directories and files
  if(is.null(data)){
    if(!check.directory(subDir= data.file)==1){
      stop(paste0('Error: ', data.file, ' does not exist!' ))
    }
  }
  
  if(!is.null(coef.file)){
    if(!check.directory(subDir= coef.file)==1){
      stop(paste0('Error: ', coef.file, ' does not exist!' ))
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
  
  if(write.coef){
    if(is.null(output.coef)){
      stop('Error: set output coef file.')
    }
    if(!create.directory(subDir= substr(output.coef, 1, regexpr('/', output.coef, fixed= TRUE)))==0){
      stop(paste0('Error: failed to create ', substr(output.coef, 1, regexpr('/', output.coef, fixed= TRUE))))
    }
  }
  
  # Read data
  if(is.null(data)){
    data <- fread(input= data.file, sep= '\t', header= TRUE, showProgress= FALSE)
  }
  setkey(data, fragment_id, run_id)
  rm(data.file)
  
  if(is.null(coef)&is.null(coef.file)){
    if(by_group){
      coef <- scaleNorm(dt= data, measure.id= measure.id, value= value.var, rep.id= 'run_id', group.id= 'tech_id')
    } else {
      coef <- scaleNorm(dt= data, measure.id= measure.id, value= value.var, rep.id= 'run_id', group.id= NULL)
      setnames(coef, 'reps', 'run_id')
    }
  } else if(is.null(coef)){
    coef <- fread(input= coef.file, sep= '\t', header= TRUE, showProgress= FALSE) 
  }
  
  if(return.table){
    setkey(data, run_id)
    setkey(coef, run_id)
    data<-data[coef]
    setnames(data, measure.id, '__modify__')
    data[, __modify__:= __modify__*coef]
    setnames(data, '__modify__', measure.id)
    data[, coef:= NULL]
  }
    
  if(write.file){
    write.table(x= data, file= output.file, quote= FALSE, sep= '\t', row.names= FALSE)
  }
  
  if(write.coef){
    write.table(x= coef, file= output.coef, quote= FALSE, sep= '\t', row.names= FALSE)
  }
  
  if(return.table&!return.coef){
    return(data)
  } else if(!return.table&return.coef){
    return(coef)
  } else if(return.table&return.coef){    
    return(list(data= data, coef= coef))
  } else if(write.file&!write.coef){
    return(output.file)
  } else if(!write.file&write.coef){
    return(output.coef)
  } else if(write.file&write.coef){
    return(list(data= output.file, coef= output.coef))
  } else {
    #return(0)
  }
  
}

normalize.ms.coef <- function(data= NULL, data.file= NULL, by_group= TRUE, output.coef= NULL, write.coef= !is.null(output.coef)){
  # Alias for special default filename arguments
  normalize.ms(data= data, coef= NULL, data.file= data.file, by_group= by_group, coef.file= NULL, output.file= NULL, return.table= FALSE, return.coef= TRUE, write.file= FALSE, write.coef= write.coef, output.coef= output.coef)  
}