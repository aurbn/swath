#
home <- '.'
data.directory <- 'data/'
data.files <- 'swathdata/'
libs.directory <- 'libs/'
sample.desc <- 'samples/'
sample.description.file <- 'samples.txt'
unique.file <- 'unique.txt'
tryptic.file <- 'tryptic.txt'
resuls.directory <- 'results/'
aggregated.data.file <- 'aggregated_data.txt'
treated.data.file <- 'treated_data.txt'
precursor.sample.score.file <- 'precursor_score_by_sample.txt'
protein.biosample.score.file <- 'protein_score_by_biosample.txt'
preclust.ms.coef.file <- 'preclust_ms_coef.txt'
preclust.mean.ms.file <- 'preclust_mean_ms.txt'
ms.coef.file <- 'ms_coef.txt'
tech.coef.file <- 'tech_coef.txt'
biosample.coef.file <- 'biosample_coef.txt'

not_from_mv <- TRUE

data.path <- paste0(data.directory, data.files)
samples.path <- paste0(data.directory, sample.desc, sample.description.file)
unique.path <- paste0(data.directory, libs.directory, unique.file)
tryptic.path <- paste0(data.directory, libs.directory, tryptic.file)
aggregated.path <- paste0(resuls.directory, aggregated.data.file)
treated.path <- paste0(resuls.directory, treated.data.file)
precursor.sample.score.path <- paste0(resuls.directory, precursor.sample.score.file)
protein.sample.score.path <- paste0(resuls.directory, protein.biosample.score.file)
preclust.ms.coef.path <- paste0(resuls.directory, preclust.ms.coef.file)
preclust.mean.ms.path <- paste0(resuls.directory, preclust.mean.ms.file)
ms.coef.path <- paste0(resuls.directory, ms.coef.file)
tech.coef.path <- paste0(resuls.directory, tech.coef.file)
biosample.coef.path <- paste0(resuls.directory, biosample.coef.file)

rm(list= c('aggregated.data.file', 'preclust.ms.coef.file', 'preclust.mean.ms.file', 'ms.coef.file', 'tech.coef.file', 'biosample.coef.file', 'data.directory', 'data.files', 'libs.directory', 'resuls.directory', 'sample.desc', 'sample.description.file', 'tryptic.file', 'treated.data.file', 'unique.file', 'precursor.sample.score.file', 'protein.biosample.score.file'))