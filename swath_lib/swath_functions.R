#
wd <- getwd()

setwd(paste(sep="/", script_dir, "swath_lib"))

require(data.table, quietly=TRUE)
require(reshape2, quietly=TRUE)
require(fpc, quietly=TRUE)

source('measure_and_rep_levels.R')
source('suppl_functions.R')
source('collect_data.R')
source('threshold_measurements.R')
source('complete_measurements.R')
source('min_measurements.R')
source('cluster_measurements.R')
source('normalize.R')
source('produce_stat.R')
source('select_modifications.R')
source('rank_groups.R')
source('score.R')
source('clustering_results.R')

setwd(wd)
rm(wd)
