write_table <- function(filename)
{
    loginfo("Writing table: %s", filename)
    if (!create.directory(dirname(output.file)==0))
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

build.check.expr <- function(flags, logic="&")
{
    r <- NULL 
    k <- c()
    for ( x in flags)
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

    return(parse(text=r))
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
    loginfo("Filtering data with vars: %s, logic=%s and remove.columns=%i", 
            paste(list(...)), logic, remove.columns)
    
    n_before = data[,.N]
    
    setkeyv(data, unlist(list(...)))
    data_ <- data[eval(build.check.expr(...,logic=logic)),]
    
    if (remove.columns)
    {
        for ( x in list(...))
        {
            data_[,eval(parse(text=x)) := NULL]
        }
    }
    
    n_after = data_[,.N]
    loginfo("%i values of %i was filtered out", n_before-n_after, n_before)
    
    data_
}

filter.apply.combine <- function(data, func, columns, ...)
{
    browser()
    fn <- match.fun(func)
    tmp <- unique(do.call(filter.measurements, c(data, columns))[,list(fragment_id, run_id)])
    setkey(tmp, fragment_id, run_id)
    tmp <- fn(data=tmp, ...)
    setkey(tmp, fragment_id, run_id)
    setkey(data, fragment_id, run_id)
    data = tmp[data]
    return(data)
}

splittmp <- function(data, flags, columns)
{
    setkeyv(data, columns)
    return(unique(data[eval(build.check.expr(flags)), columns, with = FALSE]))
}

combinetmp.n <- function(data, tmp)
{
    columns <- intersect(names(data), names(tmp))
    setkeyv(data, columns)
    setkeyv(tmp, columns)
    r = tmp[data]
    fl = setdiff(names(tmp), names(data))
    r[fl==NA, fl:=FALSE, with=FALSE]
    r
}

combinetmp <- function(data, tmp, group)
{
    fl = setdiff(names(tmp), names(data))
    ok <- unique(tmp[eval(build.check.expr(fl)), group, with=FALSE])
    ok[,eval(fl):=TRUE]
    setkeyv(ok, group)
    setkeyv(data, group)
    data = ok[data]
    data[is.na(get(fl)), (fl) := FALSE]
}

