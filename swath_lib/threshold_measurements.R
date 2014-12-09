#' Compare \code{value.var} with \code{threshhold} by \code{operator} and 
#' set the \code{flag}
#' 
#' @param data data.table
#' @param measure.id 
#' @param value.var value column to be compared
#' @param threshhold value to compare with
#' @param operator comparison operator
#' @param flag.name flag to be set
threshold.measurements <- function(data, 
                                   measure.id, 
                                   value.var, 
                                   threshold = 0,
                                   operator = '>',
                                   flag.name = 'threshold')
{
    loginfo("Thresholding with measure=%s, value=%s%s%i", measure.id, value.var,
            operator, threshold)
    
    setkeyv(data, measure.id)
        
    # Remove fragments below threshold
    ok <- unique( data[ !(get(operator)( get( value.var), threshold)),
                            get(measure.id)])
    
    data[, eval(flag.name)] <- !(data[, get(measure.id)] %in% unlist(
        ok, recursive = TRUE, use.names = FALSE))
    rm(ok)
    
    loginfo("%i under threshold and %i above", sum(data[,get(flag.name) == FALSE]),
            sum(data[,get(flag.name) == TRUE]))
    
    return(data)
}