#
cluster.measurements <- function(data,
                                 measure.id, 
                                 group.id, 
                                 rep.id, 
                                 value.var, 
                                 flag.name= 'cluster',
                                 dbscn.eps.init= 5*pi/180, 
                                 dbscn.eps.limit= 0.001, 
                                 dbscn.MinPts= 3,
                                 dbscn.step= 0.98)
{
    loginfo("Started clustering")
    # Check and create directories and files
    
    reps <- unique(data[, get(rep.id)]) # vector with sample names
    
    
    clusterizeMe <- data.table(dcast(data = data, sprintf('%s + %s ~ %s', 
                                                          eval(measure.id), 
                                                          eval(group.id), eval(rep.id)),
                                     value.var = eval(value.var)))
    
    ### For empty runs --- Check it
    clusterizeMe = na.omit(clusterizeMe) 
    
    setkeyv(clusterizeMe, c(group.id, measure.id))
    
    clusterizeMe$cluster <- FALSE # set all cluster flags to FALSE
    
    groups = unique(clusterizeMe[, get(group.id)])
    loginfo("in %i groups", length(groups))
    
    for (group in groups) {
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
            prec.dbscn <- dbscan(prec.dist, eps = dbscn.eps, 
                                 MinPts = dbscn.MinPts, method = 'dist')
            
            # calc number of clusters
            prec.clusters <- ulength(prec.dbscn$cluster[prec.dbscn$cluster > 0])
            # change eps
            dbscn.eps <- dbscn.eps * dbscn.step
        }
        
        # change cluster flag for 1-clustered transitions
        # think about this (may be faster)
        clusterizeMe[which(clusterizeMe[, get(group.id)]  == group)[which(prec.dbscn$cluster == 1)]
                     ,cluster := TRUE]
    }
    
    # clear memory
    rm(list = c('prec.mtrx', 'prec.dist', 'prec.dbscn', 'prec.clusters', 'dbscn.eps',
                'dbscn.step', 'dbscn.MinPts', 'reps'))
    
    clusterizeMe <- clusterizeMe[,c(measure.id, 'cluster'), with= FALSE]
    # ok. now we can add flag data to our swath.mwdata.all
    setkeyv(clusterizeMe, measure.id)
    setkeyv(data, measure.id)
    
    #data <- data[get(measure.id) %in% unique(unlist(clusterizeMe[cluster==TRUE, 
    #             get(measure.id)], recursive= TRUE, use.names= FALSE))]
    setnames(clusterizeMe, 'cluster', eval(flag.name))
    data <- data[clusterizeMe]
    
    setkeyv(data, group.id)

    loginfo("%i measurements of %i are clustered into %i clusters", 
                 sum(data[,get(flag.name) == TRUE]), data[,.N], 
                length(unique(data[get(flag.name) == TRUE, get(group.id)])))

    rm(clusterizeMe)

    return(data)
}