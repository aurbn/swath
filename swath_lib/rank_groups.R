#
rank.groups <- function(data= NULL, data.file= NULL, measure.id= NULL, stat.function= NULL, group.id= NULL, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
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
    
  if(is.function(stat.function)){
    stat.function <- match.fun(stat.function)
  } else {
    stat.function <- function(x){'Not a function!'}
  }
  
  setkeyv(data, c(measure.id, group.id))
  data[, rankby:= stat.function(.SD), by= list(get(measure.id), get(group.id))]
  rankUs <- unique(data[, list(get(measure.id), get(group.id), rankby)])
  setnames(rankUs, c('V1', 'V2'), c(measure.id, group.id))
  setkeyv(rankUs, c(measure.id, group.id))
  data[, rankby:= NULL]
   
  groups <- unique(rankUs[, get(group.id)])
  byprot <- sapply(groups, FUN= function(x)rankUs[get(group.id)==x], simplify= FALSE)
  byprot.sort <- sapply(byprot, simplify= FALSE, FUN= function(x)x[order(-x$rankby), ])
  byprot.sort <- sapply(byprot.sort, simplify= FALSE, FUN= function(x)x[, rank:= .I])
  rankUs <- rbindlist(byprot.sort)[, rankby:= NULL]
  setkeyv(rankUs, c(measure.id, group.id))
  setkeyv(data , c(measure.id,group.id))
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