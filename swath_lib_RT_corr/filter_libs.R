rm(list=ls())
setwd('~/swath lib RT corr/myc merged lib october january - june/')
library(package= 'data.table')

filter <- fread(input='LIBRARY_myc_merged_october_january_march_may_april_may_june_I_june_II_no_mopd_or_ccam_unique_tryptic.txt', sep='\t', header= T)

filter[, fragment_id:= paste0(protein_name, '_', modification_sequence, '_', prec_z, '_', frg_type, '_', frg_nr, '_', frg_z)]
setkey(filter, fragment_id)
filter <- unique(filter[, fragment_id])

lib <- fread(input='~/swath lib RT corr/myc june II/LIBRARY_myc_merged_october_january_march_may_april_may_june_I_june_II_RT_CORRECTED_BY_myc_swath_lib_june_II.txt', sep='\t', header= T)
lib[, fragment_id:= paste0(protein_name, '_', modification_sequence, '_', prec_z, '_', frg_type, '_', frg_nr, '_', frg_z)]
setkey(lib, fragment_id)
lib <- lib[fragment_id %in% filter]
lib[, fragment_id:= NULL]
write.table(lib, file='~/swath lib RT corr/myc june II/LIBRARY_myc_merged_october_january_march_may_april_may_june_I_june_II_RT_CORRECTED_BY_myc_swath_lib_june_II_no_mod_or_ccam_unique_tryptic.txt', sep='\t', col.names= TRUE, row.names= FALSE,na="", quote= FALSE)
