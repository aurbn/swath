#
cluster.fragments <- function(data= NULL, data.file= NULL, dbscn.eps.init= 5*pi/180, dbscn.MinPts= 3, dbscn.step= 0.98, output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){
  
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
  setkey(data, fragment_id, tech_id)
  rm(data.file)
  
  clusterizeMe <- data
  
  techs <- unique(clusterizeMe$tech_id) # vector with sample names
  
  
  clusterizeMe <- data.table(dcast(data = clusterizeMe, fragment_id + precursor_id ~ tech_id,
                                     value.var = 'intensity'))
  
  setkey(clusterizeMe, precursor_id, fragment_id)
  
  clusterizeMe$cluster <- FALSE # set all cluster flags to FALSE
  
  for (precursor in unique(clusterizeMe$precursor_id)) {
    # get value matrix
    prec.mtrx <- as.matrix(clusterizeMe[precursor, techs, with = FALSE])
    
    # calculate angle distance
    prec.dist <- apply(prec.mtrx, 1, function(s1) {
      apply(prec.mtrx, 1, function(s2) { angleDist(s1, s2)} )
    })
    
    
    dbscn.eps <- dbscn.eps.init ;  #  maximum distance
    # number of dbscan clusters. set to 2 to start cycle
    prec.clusters <- 2
    
    # run cycle
    while(prec.clusters > 1 & dbscn.eps >= 0.001) {
      # dbscan
      prec.dbscn <- dbscan(prec.dist, eps = dbscn.eps, MinPts = dbscn.MinPts, method = 'dist')
      # calc number of clusters
      prec.clusters <- ulength(prec.dbscn$cluster[prec.dbscn$cluster > 0])
      # change eps
      dbscn.eps <- dbscn.eps * dbscn.step
    }
    
    # change cluster flag for 1-clustered transitions
    # think about this (may be faster)
    clusterizeMe[which(clusterizeMe$precursor == precursor)[which(prec.dbscn$cluster == 1)],
                   cluster := TRUE]
  }
  
  # clear memory
  rm(list = c('prec.mtrx', 'prec.dist', 'prec.dbscn', 'prec.clusters', 'dbscn.eps',
              'dbscn.step', 'precursor', 'dbscn.MinPts', 'techs'))
  
  clusterizeMe <- clusterizeMe[,c('fragment_id', 'cluster'), with= FALSE]
  # ok. now we can add flag data to our swath.mwdata.all
  setkey(clusterizeMe, fragment_id)
  setkey(data, fragment_id)
  data <- data[clusterizeMe]
  rm(clusterizeMe)
  
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