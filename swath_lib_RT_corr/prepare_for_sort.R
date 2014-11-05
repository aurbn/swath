rm(list=ls())
setwd('~/swath lib RT corr/')
library(package= 'data.table')

data <- fread(input='Sergiev_ecoli_merged_lib_june_january_RT_june_no_mod_or_ccam.txt', sep='\t', header= T)
data[, clean_peptide_id:= paste0(uniprot_id, '_', stripped_sequence)]
setkey(data, clean_peptide_id)

data_for_sort <- unique(data[, list(uniprot_id, stripped_sequence, clean_peptide_id)])

write.table(data_for_sort, file='Sergiev_ecoli_merged_lib_june_january_RT_june_no_mod_or_ccam_for_sort.txt', sep='\t', col.names= TRUE, row.names= FALSE, quote= FALSE)
rm(list=ls())