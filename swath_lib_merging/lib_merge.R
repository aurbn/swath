rm(list=ls())
setwd('~/swath lib merging')
library(package= 'data.table')

files<-dir(pattern= '\\.txt$')
data.list<-list() # create list before cycle

for (file in files) {
  tmp <- fread(input= file, header= TRUE,sep= '\t', colClasses= c(Q1= 'numeric', Q3= 'numeric', uniprot_id= 'character', modification_sequence= 'character', RT_detected= 'numeric', mods= 'character', nterm= 'character', cterm= 'character'), showProgress= FALSE)
  tmp[, fragment_id:= paste0(protein_name, '_', modification_sequence, '_', prec_z, '_', frg_type, '_', frg_nr, '_', frg_z)]
  tmp[, precursor_id:= paste0(protein_name, '_', modification_sequence, '_', prec_z)]
  tmp[, peptide_id:= paste0(protein_name, '_', modification_sequence)]
  setkey(tmp)
  data.list[[file]] <- tmp
}
rm(list= c('tmp', 'file', 'files'))

data <- data.list[[1]]
peaks <- unique(data[, fragment_id])
for(i in 2:length(data.list)){
  data <- merge(data, data.list[[i]][!(fragment_id %in% peaks)], all= TRUE)
  peaks <- unique(data[, fragment_id])
}

data[, RT_detected:= mean(RT_detected), by= peptide_id]
data[, iRT:= mean(iRT), by= peptide_id]
data[, Q1:= mean(Q1), by= precursor_id]
data[, Q3:= mean(Q3), by= fragment_id]

data[, unq1:= length(unique(Q1)), by= precursor_id]
data[, unq3:= length(unique(Q3)), by= fragment_id]
data[, unirt:= length(unique(iRT)), by= peptide_id]
data[, unrt:= length(unique(RT_detected)), by= peptide_id]

if(!unique(data[,unirt])==1) stop('non-unique irt')
if(!unique(data[,unrt])==1) stop('non-unique rt')
if(!unique(data[,unq1])==1) stop('non-unique q1')
if(!unique(data[,unq3])==1) stop('non-unique q3')


data[, c('fragment_id', 'peptide_id', 'precursor_id', 'unirt', 'unrt', 'unq1', 'unq3'):= NULL]

data <- unique(data)

write.table(data, file= 'merged_lib.txt', quote=FALSE, sep= '\t', row.names= FALSE, na= '')

rm(list=ls())