#
collect_data <- function(data.directory= 'data/swathdata/', not_from_mv= FALSE, sample.description.file= 'data/samples/samples.txt', unique.file= 'data/libs/unique.txt', tryptic.file= 'data/libs/tryptic.txt', output.file= NULL, return.table= TRUE, write.file= !is.null(output.file)){

# Function for merging SWATH data for further analysis
# Version: 0.1
# Date: 18.06.2014
# Authors: Ivan Butenko, Dima Ischenko

# Check and create directories and files

if(!check.directory(subDir= data.directory)==0){
  stop(paste0('Error: ', data.directory, ' does not exist!' ))
}

if(!check.directory(subDir= sample.description.file)==1){
  stop(paste0('Error: ', sample.description.file, ' does not exist!' ))
}

files <- list.files(path= data.directory, pattern= '.txt')

if(length(files)==0){
  stop(paste0('Error: ', data.directory, ' has no data files!' ))
}

if(write.file){
  if(is.null(output.file)){
    stop('Error: set output file.')
  }
  if(!create.directory(subDir= substr(output.file, 1, regexpr('/', output.file, fixed= TRUE)))==0){
    stop(paste0('Error: failed to create ', substr(output.file, 1, regexpr('/', output.file, fixed= TRUE))))
  }
}

if(!is.null(unique.file)){
  if(!check.directory(subDir= unique.file)==1){
    stop(paste0('Error: ', unique.file, ' does not exist!' ))
  }
}
  
if(!is.null(tryptic.file)){
  if(!check.directory(subDir= tryptic.file)==1){
    stop(paste0('Error: ', tryptic.file, ' does not exist!' ))
  }
}

# Load libraries
require(data.table, quietly=TRUE)  # data table is very efficeint way to work with large tables. so we used it.
require(reshape2, quietly=TRUE)

# get all ms data files and merge its to one table
# it is good practice to read all tables in list
#  and then rbind its to a one data table with data.table::rbindlist
data.list<-list() # create list before cycle
for (file in files) {
  tmp <- fread(input= sprintf(paste0(data.directory,'%s'), file), header= TRUE,sep= '\t', showProgress= FALSE)
  if(not_from_mv){
    setnames(tmp, c('Precursor Charge', 'Ion Type', 'Fragment Charge'), c('prec_z', 'frg_type', 'frg_z'))
    tmp[, c('Peak Name'):= paste0(Protein, '.', Peptide, '.+', prec_z, frg_type, Residue, '+', frg_z), ]
    tmp[, c('Row', 'Index'):= 1:nrow(tmp)]
    tmp$Use <- TRUE
    setnames(tmp, c('Fragment MZ', 'RT', 'Protein'), c('m/z', 'Ret. Time', 'Group'))
    tmp[, c('Peptide', 'Precursor MZ', 'prec_z', 'frg_z', 'frg_type', 'Residue'):= NULL]
    setcolorder(tmp, c('Row', 'Index', 'Peak Name', 'm/z', 'Ret. Time', 'Group', 'Use', names(tmp)[!(names(tmp) %in% c('Row', 'Index', 'Peak Name', 'm/z', 'Ret. Time', 'Group', 'Use'))]))
    tmp[is.na(tmp)] <- 0
    tmp[, ]
    
  }
  data.list[[file]] <- melt(tmp[, c(3, 8:ncol(tmp)), with= FALSE], id= 'Peak Name')
}
# rbind all data.tables to one
data <- rbindlist(data.list)
# clear memory
rm(list= c('data.list', 'tmp', 'file', 'files', 'data.directory'))

# colnames renaming
setnames(data, c('PeakName', 'run_code', 'intensity'))
# add key Sample to data for fast and efficient processing
setkey(data, run_code, PeakName)


# create data.table with unique transitions from our data
fragments <- data.table(PeakName= unique(data$PeakName))

# parse peak IDs
fragments$peptide_sequence <- gsub(x= fragments$PeakName, pattern= '(.*)\\.(.*)\\..*', replacement= '\\2')
# delete modifications from seq
fragments$clean_peptide_sequence <- gsub(x= fragments$peptide_sequence, pattern= '\\[.*?\\]', replacement= '')
fragments$clean_peptide_sequence <- gsub(x= fragments$clean_peptide_sequence, pattern= '-', replacement= '', fixed= TRUE)
# select mods
fragments$modifications <- ''  # create empty column with mods
# substr mods string
fragments$modifications[grep('\\[', fragments$peptide_sequence)] <- gsub(pattern= '.*?(.{0,1}\\[.*\\]).*', x= fragments$peptide_sequence[grep('\\[', fragments$peptide_sequence)], replacement= '\\1')
# replace between mods seq with ;
fragments$modifications <- gsub(x= fragments$modifications, pattern= '(\\]).*?(.{1}\\[)', replacement= '\\1;\\2')
# Get protein IDs
fragments$protein_id <- gsub(x= fragments$PeakName, pattern= '(.*)\\.(.*)\\..*', replacement= '\\1')
# get precursor charge
fragments$precursor_charge <- gsub(x= fragments$PeakName, pattern= '.*\\.([\\+\\-]\\d+)([y|b]\\d+)[\\+|\\-].*', replacement= '\\1')
# Get fragment properties
fragments$fragment <- gsub(x= fragments$PeakName, pattern = '.*\\.[\\+\\-]\\d+([y|b]\\d+)[\\+|\\-].*', replacement= '\\1')
fragments$fragment_series <- gsub(fragments$fragment, pattern= '[[:digit:]]', replacement= '')
fragments$fragment_number <- gsub(fragments$fragment, pattern= '[[:alpha:]]', replacement= '')
fragments$fragment <- NULL
fragments$fragment_charge <- gsub(x = fragments$PeakName, pattern= '.*\\.([\\+\\-]\\d+)([y|b]\\d+)([\\+|\\-].*)', replacement= '\\3')

setkey(fragments, PeakName)


# split Sample into sample and ms
loaded.samples <- data.table(run_code = unique(data$run_code))
loaded.samples$sample_code <- gsub(x= loaded.samples$run_code, pattern= '(.*)\\_.*', replacement= '\\1')
loaded.samples$ms <- gsub(x= loaded.samples$run_code, pattern= '(.*)\\_(.*)', replacement= '\\2')
setkey(loaded.samples, sample_code)
sample.codes <- unique(loaded.samples$sample_code) # vector with sample names

# read samples data
sample.description <- fread(input= sample.description.file, header= TRUE, sep= '\t')
setkey(sample.description, sample) #change description files and code later
setnames(sample.description, 'sample', 'sample_code')#change description files and code later
rm(sample.description.file)

# merge our samples with samples from file
loaded.samples <-loaded.samples[sample.description, nomatch= 0]

loaded.samples[,c('run_id', 'tech_id', 'bio_sample','nms'):= list(paste0(bio, '_', state,'_', tech, '_', ms), paste0(bio, '_', state, '_', tech), paste0(bio, '_', state), NULL)]

rm(list= c('sample.codes', 'sample.description'))

fragments[,peptide_id:= paste0(protein_id, '_', peptide_sequence)]
fragments[,clean_peptide_id:= paste0(protein_id, '_', clean_peptide_sequence)]
fragments[,precursor_id:= paste0(peptide_id, '_', precursor_charge)]
fragments[,fragment_id:= paste0(precursor_id, '_', fragment_series, '_', fragment_number, '_', fragment_charge)]

setkey(fragments, PeakName)
setkey(loaded.samples, run_code)

setkey(data, run_code)
data<-data[loaded.samples]

setkey(data, PeakName)
data<-data[fragments]

setkey(data, fragment_id, run_id)

rm(list= c('fragments', 'loaded.samples'))
if(!is.null(unique.file)){
  uniq <- fread(input= unique.file, sep= '\t', header= TRUE, showProgress= FALSE)
  if (!not_from_mv){
     setnames(uniq, c('ProteinName', 'NoModSeq'), c('protein_id', 'clean_peptide_sequence'))
     uniq[, clean_peptide_id:= paste0(protein_id, '_', clean_peptide_sequence)]
  }
  setkey(uniq)
  uniq <- unique(uniq[, clean_peptide_id])
  setkey(data, clean_peptide_id)
  data$uniq <- data[, clean_peptide_id] %in% uniq
}

if(!is.null(tryptic.file)){
  tryptic <- fread(input= tryptic.file, sep= '\t', header= TRUE, showProgress= FALSE)
  if(!not_from_mv){
    setnames(tryptic, c('ProteinName', 'NoModSeq'), c('protein_id', 'clean_peptide_sequence'))
    
    tryptic[, clean_peptide_id:= paste0(protein_id, '_', clean_peptide_sequence)]
  }
  setkey(tryptic)
  tryptic <- unique(tryptic[, clean_peptide_id])
  setkey(data, clean_peptide_id)
  data$tryptic <- data[, clean_peptide_id] %in% tryptic
}


if(write.file){
  write.table(x= data, file= output.file, quote= FALSE, sep= '\t', row.names= FALSE)
}

if(return.table){
  return(data)
} else {
    if(write.file){
      return(output.file)
    } else {
      #return(0)
    }
}

}