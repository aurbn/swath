#' Leave only fragments detected in all runs or in % of rins 
#' 
#' If this function is called with \code{detect=0} 
#' it leaves only fragments detected in all runs. 
#' Or it leaves only fragments detected at least \code{detect} times or detected 
#' in \code{detect} percent of cases, according to \code{percent} flage.
#' 
#' @param data data table
#' @param flag
#' @param measure.id
#' @param rep.id
#' @param detect
#' @param percent
complete.measurements <- function(data,
                                  measure.id= NULL,
                                  rep.id= NULL,
                                  detect=0,
                                  percent=FALSE,
                                  flag.name= 'complete')
{
    # Load libraries
    require(data.table, quietly=TRUE) 
        
    # count each fragment occurencies
    # http://stackoverflow.com/questions/19869145/counting-in-r-data-table
    setkeyv(data, c(measure.id))
    data[ , `:=`( times_detected = .N ) , by =  get(measure.id)]

    #TODO(urban): Probably slow, check & rewrite
    setkeyv(data, c(rep.id))
    nreps = length(unique(data[, get(rep.id)]))

    if (detect > 0)
    {
        if (percent)
        {
            data[,get(flag.name) := times_detected / nreps >= detect/100.0]    
        }else
        {
            data[,get(flag.name) := times_detected >= detect]    
        }  
    }else
    {
        data[,get(flag.name) := times_detected == nreps]    
    }
    
    data[, times_detected:= NULL]
        
    return(data)
}