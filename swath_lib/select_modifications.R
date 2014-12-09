#' Select modifications and set a flag
#' 
#' @param data data table
#' @param sel.modifications vector of modification to be selected
#' @param flag.name flag will be added to the data

select.modifications <- function(data, sel.modifications = c("","C[CAM]"), 
                                 flag.name= 'selected_modifications')
{
    data[, eval(flag.name):= (modifications %in% sel.modifications)]
    
    return(data)
}