#' Leave only fragments detected in all runs or in % of rins 
#' 
#' If this function is called with \code{detect=0} 
#' it leaves only fragments detected in all runs. 
#' Or it leaves only fragments detected at least \code{detect} times or detected 
#' in \code{detect} percent of cases, according to \code{percent} flage.
#' 
#' @param data data table
#' @param measure.id identified fragment id
#' @param rep.id experimnt id
#' @param detect amount aof detected events
#' @param percent is \code{detect} percent or not
#' @param flag.name name of the column with results
complete.measurements <- function(data,
                                  measure.id= NULL,
                                  rep.id= NULL,
                                  detect=0,
                                  percent=FALSE,
                                  flag.name= 'complete')
{
    loginfo("Checking data completeness with detect=%i %s", as.integer(detect),
            ifelse(percent, "percent.", "measurements."))
    
    # count each fragment occurencies
    setkeyv(data, rep.id)
    
    #data[ , times_detected := .N  , by =  eval(measure.id)] #FAST
    data[ , times_detected := length(unique(get(rep.id))), by =  eval(measure.id)] #Reusable

    #TODO(urban): Probably slow, check & rewrite
    setkeyv(data, c(rep.id))
    nreps = length(unique(data[, get(rep.id)]))

    if (detect > 0)
    {
        if (percent)
        {
            data[,eval(flag.name) := times_detected / nreps >= detect/100.0]    
        }else
        {
            data[,eval(flag.name) := times_detected >= detect]    
        }  
    }else
    {
        data[,eval(flag.name) := times_detected == nreps]    
    }
    
    data[, times_detected:= NULL]
    
    loginfo("Of %i values %i are complete and %i are incomplete.", data[, .N],
            data[get(flag.name)==TRUE, .N], data[get(flag.name)==FALSE, .N])
        
    return(data)
}