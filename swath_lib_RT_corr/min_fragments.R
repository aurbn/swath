rm(list=ls())
setwd('~/swath lib RT corr/')
library(package= 'data.table')
lib <- fread(input='~/swath lib RT corr/Sergiev_ecoli_merged_lib_june_january_RT_june.txt', sep='\t', header= T)

lib[, precursor_id:= paste0(protein_name, '_', modification_sequence, '_', prec_z)]
lib[, fragment_id:= paste0(protein_name, '_', modification_sequence, '_', prec_z, '_', frg_type, '_', frg_nr, '_', frg_z)]

lib[, frg_by_prec:= length(unique(fragment_id)), by= precursor_id]

summary <- unique(lib[, list(precursor_id, frg_by_prec)])
hist(summary$frg_by_prec,breaks=max(summary$frg_by_prec))
nrow(summary[frg_by_prec==60])
