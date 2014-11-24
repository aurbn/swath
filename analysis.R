source(paste(sep="/", script_dir, "swath_lib/swath_functions.R"))

#Aggregate all SWATH data
data <- collect_data(data.directory=data.path,
                     sample.description.file=samples.path, 
                     unique.file=unique.path, 
                     tryptic.file=tryptic.path, 
                     not_from_mv=not_from_mv)

data <- complete.measurements(data= tmp,
                              ndetect=0
                              flag= FALSE,
                              flag.name= NULL,
                              measure.id= 'fragment_id', 
                              rep.id= 'run_id',
                              return.table= TRUE)