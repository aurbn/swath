#
cluster.measurements <- function(data= NULL, flag= TRUE, flag.name= 'cluster', measure.id= NULL, group.id= NULL, rep.id= NULL, value.var= NULL, data.file= NULL, dbscn.eps.init= 5*pi/180, dbscn.eps.limit= 0.001, dbscn.MinPts= 3, dbscn.step= 0.98, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
  # Load libraries
  require(data.table, quietly=TRUE)  # data table is very efficeint way to work with large tables. so we used it.
  require(fpc, quietly=TRUE)
  
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
  setkeyv(data, c(measure.id, rep.id))
  rm(data.file)
  
  
  reps <- unique(data[, get(rep.id)]) # vector with sample names
  
  
  clusterizeMe <- data.table(dcast(data = data, sprintf('%s + %s ~ %s', eval(measure.id), eval(group.id), eval(rep.id)),
                                     value.var = eval(value.var)))
  
  setkeyv(clusterizeMe, c(group.id, measure.id))
  
  clusterizeMe$cluster <- FALSE # set all cluster flags to FALSE
  
  for (group in unique(clusterizeMe[, get(group.id)])) {
    # get value matrix
    prec.mtrx <- as.matrix(clusterizeMe[group, reps, with = FALSE])
    
    # calculate angle distance
    prec.dist <- apply(prec.mtrx, 1, function(s1) {
      apply(prec.mtrx, 1, function(s2) { angleDist(s1, s2)} )
    })
    
    
    dbscn.eps <- dbscn.eps.init ;  #  maximum distance
    # number of dbscan clusters. set to 2 to start cycle
    prec.clusters <- 2
    
    # run cycle
    while(prec.clusters > 1 & dbscn.eps >= dbscn.eps.limit) {
      # dbscan
      prec.dbscn <- dbscan(prec.dist, eps = dbscn.eps, MinPts = dbscn.MinPts, method = 'dist')
      # calc number of clusters
      prec.clusters <- ulength(prec.dbscn$cluster[prec.dbscn$cluster > 0])
      # change eps
      dbscn.eps <- dbscn.eps * dbscn.step
    }
    
    # change cluster flag for 1-clustered transitions
    # think about this (may be faster)
    clusterizeMe[which(clusterizeMe[, get(group.id)]  == group)[which(prec.dbscn$cluster == 1)],
                   cluster := TRUE]
  }
  
  # clear memory
  rm(list = c('prec.mtrx', 'prec.dist', 'prec.dbscn', 'prec.clusters', 'dbscn.eps',
              'dbscn.step', 'dbscn.MinPts', 'reps'))
  
  clusterizeMe <- clusterizeMe[,c(measure.id, 'cluster'), with= FALSE]
  # ok. now we can add flag data to our swath.mwdata.all
  setkeyv(clusterizeMe, measure.id)
  setkeyv(data, measure.id)
    
  if(!flag){
    data <- data[get(measure.id) %in% unique(unlist(clusterizeMe[cluster==TRUE, get(measure.id)], recursive= TRUE, use.names= FALSE))]
  } else if(!is.null(flag.name)){
    setnames(clusterizeMe, 'cluster', eval(flag.name))
    data <- data[clusterizeMe]
    rm(clusterizeMe)
  } else {
    stop('Error flag name needed')
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