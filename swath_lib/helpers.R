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