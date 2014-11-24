source(paste(sep="/", script_dir, "swath_lib/swath_functions.R"))

# Aggregate all SWATH data
data <- collect.data(data.directory = data.path,
                     sample.description.file = samples.path, 
                     unique.file = unique.path, 
                     tryptic.file = tryptic.path, 
                     not_from_mv = not_from_mv)

# Leave only fragments detected in all runs
data <- complete.measurements(data,
                              measure.id= "fragment_id",
                              rep.id= "run_id")


