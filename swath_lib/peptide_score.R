#
peptide.score <- function(data= NULL, data.file= NULL, measure_id= 'bio_sample', output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
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
  rm(data.file)
  
  
  keys <- c('peptide_id', measure_id)
  setkeyv(data, keys)
  
  data[, peptide_score:= sum(intensity), by= keys]
  data[, prec_num:= ulength(precursor_id), by= peptide_id]
  data[, frg_num:= ulength(fragment_id), by= peptide_id]
  data <- unique(data[, c('protein_id', 'peptide_id', 'prec_num', 'frg_num', measure_id, 'peptide_score'), with= FALSE])
  
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