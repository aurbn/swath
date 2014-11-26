write_table <- function(filename)
{
    loginfo("Writing table: %s", filename)
    if (!create.directory(dirname(output.file)==0)
    {
        logerror('Error: failed to create ', dirname(output.file))
        stop()
    }
    write.table(x=data, file=output.file, quote= FALSE, sep='\t', row.names=FALSE)
}

read_table <- function(filename)
{
    data <- fread(input= data.file, sep= '\t', header= TRUE, showProgress= FALSE)
} 

#' Filter data.table rows according to logical columns
#' 
#' Construct a new data.table with rows from \code{data}
#' where all (\code{logic="&"}) or any (\code{logic="|"}) columns 
#' from \code{...} are true. Columns must be logical. 
#' Also it is optiannaly removes these columns from the data.
#' @param data data.table
#' @param ... one or more column names
#' @param logic & or |
#' @param remove.columns remove selected columns or not
filter.measurements <- function(data, ..., logic="&", remove.columns=TRUE)
{
    r <- NULL 
    k <- c()
    for ( x in list(...))
    {
        k <- c(k, x)
        t <- paste0(x,"==TRUE")
        if (is.null(r))
        {
            r <- t
        }else
        {
            r <- paste(sep=logic, r, t)
        }
    }
    
    loginfo("Filtering data with vars: %s, logic=%s and remove.columns=%i", 
            paste(k), logic, remove.columns)
    n_before = data[,.N]
            
    setkeyv(data, k)
    data_ <- data[eval(parse(text=r)),]
    
    if (remove.columns)
    {
        for ( x in k)
        {
            data_[,eval(parse(text=x)) := NULL]
        }
    }
    
    n_after = data_[,.N]
    loginfo("%i values of %i was filtered out", n_before-n_after, n_before)
    
    data_
}
