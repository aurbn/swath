#
rm(list=ls())
setwd('~/swath lib RT corr/myc merged lib october january - june/')
library(package= 'data.table')

data <- fread(input='LIBRARY_myc_merged_october_january_march_may_april_may_june_I_june_II.txt', sep='\t', header= T)
data[, clean_peptide_id:= paste0(protein_name, '_', stripped_sequence)]
setkey(data, clean_peptide_id)

data_for_sort <- unique(data[, list(protein_name, stripped_sequence, clean_peptide_id)])

write.table(data_for_sort, file='LIBRARY_myc_merged_october_january_march_may_april_may_june_I_june_II_for_sort.txt', sep='\t', col.names= TRUE, row.names= FALSE, quote= FALSE)
rm(list=ls())