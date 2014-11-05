#
average.tech <- function(data= NULL, data.file= NULL, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
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
  
  setkey(data, fragment_id, run_id)
  rm(data.file)
  
  # Average tech repeats
  
  data <- average.ms(data)
  
  data[, c('mean','se','cv'):= list(mean(intensity), se_corr(intensity), se_corr(intensity)/mean(intensity)), by=c('fragment_id', 'bio_sample')]
  data <- unique(data[, c('intensity', 'tech_id', 'tech', 'sample_code'):= NULL])
  setnames(data, 'mean', 'intensity')
    
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