#' Leave only fragments detected in all runs
#' 
#' @param data data table
#' @param flag
complete.measurements <- function(data,
                                  measure.id= NULL,
                                  rep.id= NULL,
                                  ndetect=0,
                                  flag= TRUE,
                                  flag.name= 'complete')
{
    
    # Load libraries
    require(data.table, quietly=TRUE) 
    
    setkeyv(data, c(measure.id, rep.id))
    
    # Keep only fragments detected en each run
    nreps <- ulength(unlist(data[, get(rep.id)], recursive= TRUE, use.names= FALSE))
    data[, times_detected:= ulength(get(rep.id)), by= measure.id]
    
    if(!flag){
        data <- data[times_detected==nreps]
    } else if(!is.null(flag.name)){
        data[, eval(flag.name)] <- data[, times_detected]==nreps
    } else {
        stop('Error flag name needed')
    } 
    
    data[, times_detected:= NULL]
    #rm(list= c('fragments', 'fragment.stat', 'nruns'))
    
    return(data)
}