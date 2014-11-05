#
normalize <- function(data= NULL, measure.id= NULL, value.var= NULL, rep.id= NULL, group.id= NULL, coef= NULL, data.file= NULL, coef.file= NULL, output.file= NULL, return.table= TRUE, return.coef= FALSE, write.file= !is.null(output.file), output.coef= NULL, write.coef= !is.null(output.coef)){
  
  # Load libraries
  require(data.table, quietly=TRUE)  # data table is very efficeint way to work with large tables. so we used it.
  
  # Check and create directories and files
  if(is.null(data)){
    if(!check.directory(subDir= data.file)==1){
      stop(paste0('Error: data file', data.file, ' does not exist!' ))
    }
  }
  
  if(!is.null(coef.file)){
    if(!check.directory(subDir= coef.file)==1){
      stop(paste0('Error: coef file ', coef.file, ' does not exist!' ))
    }
  }
  
  if(write.file){
    if(is.null(output.file)){
      stop('Error: set output file.')
    }
    if(!create.directory(subDir= substr(output.file, 1, regexpr('/', output.file, fixed= TRUE)))==0){
      stop(paste0('Error: failed to create output directory ', substr(output.file, 1, regexpr('/', output.file, fixed= TRUE))))
    }
  }
  
  if(write.coef){
    if(is.null(output.coef)){
      stop('Error: set output coef file.')
    }
    if(!create.directory(subDir= substr(output.coef, 1, regexpr('/', output.coef, fixed= TRUE)))==0){
      stop(paste0('Error: failed to create output directory', substr(output.coef, 1, regexpr('/', output.coef, fixed= TRUE))))
    }
  }
  
  # Read data
  if(is.null(data)){
    data <- fread(input= data.file, sep= '\t', header= TRUE, showProgress= FALSE)
  }
  setkeyv(data, c(measure.id, rep.id, group.id))
  rm(data.file)
  
  if(is.null(coef)&is.null(coef.file)){
    coef <- scaleNorm(dt= data, measure.id= measure.id, value= value.var, rep.id= rep.id, group.id= group.id)    
  } else if(is.null(coef)){
    coef <- fread(input= coef.file, sep= '\t', header= TRUE, showProgress= FALSE)
  }
  
  if(return.table|write.file){
    setkeyv(data, c(rep.id, group.id))
    setkeyv(coef, c(rep.id, group.id))
    data<-data[coef]
    data[, eval(value.var):= get(value.var)*coef, by= c(measure.id, rep.id)]
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
    return(0)
  }
  
}