#!/usr/bin/env Rscript

# Check if all packages are present and install them otherwise
list.of.packages <- c("data.table", "fpc", "reshape2", "argparser", "logging")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) 
{
    install.packages(new.packages)
}

require(logging)
basicConfig()
addHandler(writeToFile, file="log.txt")

logdebug("Checking required packages...")
logdebug("Packages are OK.")

if(commandArgs(trailingOnly=FALSE)[1] != "RStudio")
{   
    #Try to determine current script location
    initial.options <- commandArgs(trailingOnly = FALSE)
    script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
    script_dir <- dirname(script.name)
        
    require(argparser, quietly=TRUE)
    parser <- arg.parser("SWATH data analyzer.")
    parser <- add.argument(parser, "--data_dir", help="data directory",
                           default="data/swathdata/")
    parser <- add.argument(parser, "--samples", help="samples description file",
                           default="data/samples/samples.txt")
    parser <- add.argument(parser, "--unique", help="unique file",
                           default="data/libs/unique.txt")
    parser <- add.argument(parser, "--tryptic", help="tryptic file",
                           default="data/libs/tryptic.txt")
    parser <- add.argument(parser, "--results_dir", help="base directpry fro output files",
                           default="results/")
    parser <- add.argument(parser, "--mv", help="flag for old file format",
                           flag=TRUE, default=FALSE)
    parser <- add.argument(parser, "--pdf", help="create pdf with results of clustering",
                           flag=TRUE, default=FALSE)
    parser <- add.argument(parser, "--hungry", help="hungry mode:consumes lots of memory",
                           flag=TRUE, default=FALSE)
    parser <- parse.args(parser, argv = commandArgs(trailingOnly = TRUE))
    
    data.files <- parser$data_dir
    sample.description.file <- parser$samples
    unique.file <- parser$unique
    tryptic.file <- parser$tryptic
    resuls.directory <- parser$results_dir
    cluster_pdf <- parser$pdf
    not_from_mv <- !parser$mv
    
    #TODO(urban): make it adjustable
    aggregated.data.file <- 'aggregated_data.txt'
    treated.data.file <- 'treated_data.txt'
    precursor.sample.score.file <- 'precursor_score_by_sample.txt'
    protein.biosample.score.file <- 'protein_score_by_biosample.txt'
    preclust.ms.coef.file <- 'preclust_ms_coef.txt'
    preclust.mean.ms.file <- 'preclust_mean_ms.txt'
    ms.coef.file <- 'ms_coef.txt'
    tech.coef.file <- 'tech_coef.txt'
    biosample.coef.file <- 'biosample_coef.txt'
    
    data.path <- data.files
    samples.path <- sample.description.file
    unique.path <- unique.file
    tryptic.path <- tryptic.file
    aggregated.path <- paste0(resuls.directory, aggregated.data.file)
    treated.path <- paste0(resuls.directory, treated.data.file)
    precursor.sample.score.path <- paste0(resuls.directory, precursor.sample.score.file)
    protein.sample.score.path <- paste0(resuls.directory, protein.biosample.score.file)
    preclust.ms.coef.path <- paste0(resuls.directory, preclust.ms.coef.file)
    preclust.mean.ms.path <- paste0(resuls.directory, preclust.mean.ms.file)
    ms.coef.path <- paste0(resuls.directory, ms.coef.file)
    tech.coef.path <- paste0(resuls.directory, tech.coef.file)
    biosample.coef.path <- paste0(resuls.directory, biosample.coef.file)
} else # Run from RStudio
{
    rm(list= ls())
    script_dir <- "."
    
    #Default paramen
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
    
    cluster_pdf = TRUE
    not_from_mv = FALSE
}

source(paste(sep="/", script_dir, "analysis.R"), verbose = T)
