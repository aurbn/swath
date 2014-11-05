rm(list=ls())
setwd('~/swath lib RT corr/')
library(package= 'data.table')

data <- fread(input='Sergiev_ecoli_merged_lib_june_january_RT_june.txt', sep='\t', header= T)
data[, peptide_id:= paste0(protein_name, '_', modification_sequence)]
setkey(data, peptide_id)
data_clean <- unique(data[,list(protein_name, peptide_id, modification_sequence, stripped_sequence)])

# select mods
data_clean$modifications <- ''  # create empty column with mods
# substr mods string
data_clean$modifications[grep('\\[', data_clean$modification_sequence)] <- gsub(pattern= '.*?(.{0,1}\\[.*\\]).*', x= data_clean$modification_sequence[grep('\\[', data_clean$modification_sequence)], replacement= '\\1')
# replace between mods seq with ;
data_clean$modifications <- gsub(x= data_clean$modifications, pattern= '(\\]).*?(.{1}\\[)', replacement= '\\1;\\2')
data_clean$modifications_rm <- gsub('C[CAM]','',data_clean$modifications, fixed= TRUE)
data_clean$modifications_rm <- gsub(';','',data_clean$modifications_rm, fixed= TRUE)
data_clean$nomodorccam <- data_clean$modifications_rm == ''

select <- unique(data_clean[nomodorccam == TRUE, peptide_id])

data_selected <- data[peptide_id %in% select]

write.table(data_selected, file='Sergiev_ecoli_merged_lib_june_january_RT_june_no_mod_or_ccam.txt', sep='\t',na="", col.names= TRUE, row.names= FALSE, quote= FALSE)
rm(list=ls())