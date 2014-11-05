#
rank.peptides <- function(data= NULL, data.file= NULL, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
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
  
  setkey(data, peptide_id)
  
  data[, frgsum:= sum(intensity), by= peptide_id]
  rankUs <- unique(data[, list(protein_id, peptide_id, frgsum)])
  setnames(rankUs, 'frgsum', 'rankby')
  setkey(rankUs, protein_id, peptide_id)
  data[, frgsum:= NULL]
   
  prots <- unique(rankUs$protein_id)
  byprot <- sapply(prots, FUN= function(x)rankUs[protein_id==x], simplify= FALSE)
  byprot.sort <- sapply(byprot, simplify= FALSE,FUN= function(x)x[order(-x$rankby), ])
  byprot.sort <- sapply(byprot.sort, simplify= FALSE, FUN= function(x)x[, rank:= .I])
  rankUs <- rbindlist(byprot.sort)[, rankby:= NULL]
  setkey(rankUs, protein_id, peptide_id)
  setkey(data , protein_id, peptide_id)
  data <- data[rankUs]
  
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