rm(list=ls())
setwd('~/swath lib RT corr')
library(package= 'data.table')

data <- fread(input='Sergiev_ecoli_merged_lib_june_january_RT_june_no_mod_or_ccam.txt', sep='\t', header= T)
data[, clean_peptide_id:= paste0(uniprot_id, '_', stripped_sequence)]
setkey(data, clean_peptide_id)

uniq <-  fread(input='Sergiev_ecoli_merged_lib_june_january_RT_june_no_mod_or_ccam_for_sort-Unique.txt', sep='\t', header= T)
tryptic <-  fread(input='Sergiev_ecoli_merged_lib_june_january_RT_june_no_mod_or_ccam_for_sort-tryptic.txt', sep='\t', header= T)

uniq <- unique(uniq[, clean_peptide_id])
tryptic <- unique(tryptic[, clean_peptide_id])

data_ut<-data[(clean_peptide_id %in% uniq)&(clean_peptide_id %in% tryptic)]
data_ut[, c('clean_peptide_id', 'peptide_id'):= NULL]

write.table(data_ut, file='LIBRARY_Sergiev_ecoli_merged_lib_june_january_RT_june_no_mod_or_ccam_unique_tryptic.txt', sep='\t', col.names= TRUE, row.names= FALSE,na="", quote= FALSE)
rm(list=ls())
