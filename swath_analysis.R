#
rm(list= ls())
#setwd("~/swath/mycoplasma_revamped/")

## Write pdf file with clustering results or not?
CLUSTER_PDF=F

## Old file format
not_from_mv <- F


source('./swath_lib/swath_functions.R')

source('./default_settings.R')
#setwd(home)

data <- collect_data(data.directory= data.path,
                     sample.description.file= samples.path, 
                     unique.file= unique.path, 
                     tryptic.file= tryptic.path, 
                     not_from_mv= not_from_mv)
#data <- dt.read(aggregated.path)

setkey(data, fragment_id, run_id)

## Keep only fragments intensity > 0 (??) 
data <- threshold.measurements(data= data,
                               flag.name= 'above_zero',
                               measure.id= 'fragment_id', 
                               value.var= 'intensity', 
                               threshold= 0, 
                               direction= 'lower') 
setkey(data, fragment_id, run_id)

##??
tmp <- unique(data[above_zero==TRUE, list(fragment_id, run_id)])
setkey(tmp, fragment_id, run_id)

## run_id tech_id

## Fragments detected in all runs
complete <- complete.measurements(data= tmp, 
                                  flag= FALSE,
                                  flag.name= NULL,
                                  measure.id= 'fragment_id', 
                                  rep.id= 'run_id',
                                  return.table= TRUE)
setkey(complete, fragment_id)
complete <- unique(complete[, fragment_id])

setkey(data, fragment_id)
data$complete <- data[, fragment_id] %in% complete #slow  
rm(complete)
rm(tmp)

setkey(data, fragment_id, precursor_id)
tmp <- unique(data[above_zero==TRUE&complete==TRUE, list(fragment_id, precursor_id)])
setkey(tmp, fragment_id, precursor_id)

## Keep only precursors with more or equal then 'min' fragments
min_frg <- min.measurements(data= tmp,
                            min= 3, 
                            measure.id= 'fragment_id',
                            group.id= 'precursor_id')

setkey(min_frg, precursor_id)
min_frg <- unique(min_frg[, precursor_id])
data$min_frg <- data[, precursor_id] %in% min_frg #slow
rm(min_frg)
rm(tmp)

setkey(data, peptide_id)
data <- select.modifications(data)

setkey(data, fragment_id, run_id, tech_id)
tmp <- unique(data[complete==TRUE&min_frg==TRUE&above_zero==TRUE,
                   list(fragment_id, run_id, tech_id, intensity)])


##?? Normalization (Good)
setkey(tmp, fragment_id, run_id, tech_id)
coef_run_preclust <- normalize(data= tmp, 
                               measure.id= 'fragment_id', 
                               value.var= 'intensity', 
                               rep.id= 'run_id',
                               group.id= 'tech_id', 
                               return.table= FALSE, 
                               return.coef= TRUE, 
                               output.coef= preclust.ms.coef.path)

##?? Normalization (All)
setkey(coef_run_preclust, tech_id, run_id)
preclust_data <- normalize(data= data,
                           measure.id= 'fragment_id',
                           value.var= 'intensity', 
                           rep.id= 'run_id',
                           group.id= 'tech_id',
                           coef= coef_run_preclust)
rm(tmp)

setkey(preclust_data, precursor_id, fragment_id, tech_id, run_id)
tmp <- unique(preclust_data[, list(fragment_id, precursor_id, tech_id, run_id, intensity,
                                   complete, min_frg, above_zero)])
rm(preclust_data)

## Calculate Mean, StdErr, StdDev of intensities over each tech id(??)
#super-slow 1 - divide for cv separately  2 - produce.stat
#subsets data table for each function separately
for_selection <-  produce.stat(data= tmp,
                               measure.id= 'fragment_id',
                               value.var= 'intensity',
                               stat.name= c('mean', 'se_corr', 'cv'),
                               group.id= 'tech_id',
                               stat.function= list(mean, se_corr, function(x){se_corr(x)/mean(x)}))

for_selection[, c('intensity', 'mean'):= list(mean, NULL)]
for_selection <- unique(for_selection[, names(for_selection)[names(for_selection) %in% ms.rep]:= NULL])
write.file(data= for_selection, path= preclust.mean.ms.path)

setkey(data, precursor_id, fragment_id)
setkey(for_selection, fragment_id, precursor_id, tech_id)
tmp <- unique(for_selection[complete==TRUE&min_frg==TRUE&above_zero==TRUE, list(fragment_id, precursor_id, tech_id, intensity)])
setkey(tmp, fragment_id, precursor_id, tech_id)

## Clustering fragments of each precursor 
clustered <- cluster.measurements(data= tmp, 
                                  flag.name= 'cluster',
                                  measure.id= 'fragment_id',
                                  group.id= 'precursor_id',
                                  rep.id= 'tech_id',
                                  value.var= 'intensity',
                                  dbscn.eps.init= 5*pi/180, 
                                  dbscn.eps.limit= 0.001, 
                                  dbscn.MinPts= 3, 
                                  dbscn.step= 0.98)

setkey(clustered, fragment_id)
clustered <- unique(clustered[cluster==TRUE, fragment_id])
data$cluster <- data[, fragment_id] %in% clustered
rm(clustered)

## Precursors which is clustered and have not less than 3 fragments
setkey(data, fragment_id, precursor_id)
tmp <- unique(data[complete==TRUE&min_frg==TRUE&above_zero==TRUE&cluster==TRUE, list(fragment_id, precursor_id)])
setkey(tmp, fragment_id, precursor_id)
min_frg_clust <- min.measurements(data= tmp, 
                                  min= 3, 
                                  measure.id= 'fragment_id',
                                  group.id= 'precursor_id')

setkey(min_frg_clust, precursor_id)
min_frg_clust <- unique(min_frg_clust[, precursor_id])
data$min_frg_clust<- data[, precursor_id] %in% min_frg_clust
rm(min_frg_clust)
rm(tmp)

## Normalize one more time ????
setkey(data, fragment_id, run_id, tech_id)
tmp <- unique(data[complete==TRUE&min_frg==TRUE&above_zero==TRUE&cluster==TRUE&min_frg_clust==TRUE, 
                   list(fragment_id, run_id, tech_id, intensity)])

setkey(tmp, fragment_id, run_id, tech_id)
coef_run <- normalize(data= tmp, 
                      measure.id= 'fragment_id',
                      value.var= 'intensity',
                      rep.id= 'run_id',
                      group.id= 'tech_id',
                      return.table= FALSE, 
                      return.coef= TRUE, 
                      output.coef= ms.coef.path)

setkey(coef_run, tech_id, run_id)
data <- normalize(data= data, 
                  measure.id= 'fragment_id',
                  value.var= 'intensity',
                  rep.id= 'run_id',
                  group.id= 'tech_id',
                  coef= coef_run)
rm(tmp)

## And again ????
setkey(data, fragment_id, run_id, tech_id, bio_sample)
tmp <- unique(data[complete==TRUE&min_frg==TRUE&above_zero==TRUE&cluster==TRUE&min_frg_clust==TRUE, 
                   list(fragment_id, run_id, tech_id, bio_sample, intensity)])

setkey(tmp, fragment_id, run_id, tech_id, bio_sample)
tmp <- produce.stat(data= tmp, 
                    measure.id= 'fragment_id',
                    value.var= 'intensity',
                    stat.name= c('mean', 'se_corr', 'cv'), 
                    group.id= 'tech_id',
                    stat.function= list(mean, se_corr, function(x){se_corr(x)/mean(x)}))

tmp[, c('intensity', 'mean'):= list(mean, NULL)]
tmp <- unique(tmp[, names(tmp)[names(tmp) %in% ms.rep]:= NULL])

setkey(tmp, fragment_id, tech_id, bio_sample)
coef_tech <- normalize(data= tmp, 
                       measure.id= 'fragment_id',
                       value.var= 'intensity',
                       rep.id= 'tech_id',
                       group.id= 'bio_sample',
                       return.table= FALSE, 
                       return.coef= TRUE, 
                       output.coef= tech.coef.path)

setkey(coef_tech, tech_id, bio_sample)
data <- normalize(data= data, 
                  measure.id= 'fragment_id',
                  value.var= 'intensity',
                  rep.id= 'tech_id',
                  group.id= 'bio_sample',
                  coef= coef_tech)
rm(tmp)

setkey(data, fragment_id, run_id, tech_id, bio_sample)
tmp <- unique(data[complete==TRUE&min_frg==TRUE&above_zero==TRUE&cluster==TRUE&min_frg_clust==TRUE, 
                   list(fragment_id, run_id, tech_id, bio_sample, intensity)])

setkey(tmp, fragment_id, run_id, tech_id, bio_sample)
tmp <- produce.stat(data= tmp, 
                    measure.id= 'fragment_id',
                    value.var= 'intensity',
                    stat.name= c('mean', 'se_corr', 'cv'), 
                    group.id= 'bio_sample',
                    stat.function= list(mean, se_corr, function(x){se_corr(x)/mean(x)}))

setkey(tmp, fragment_id, tech_id, bio_sample)
tmp[, c('intensity', 'mean'):= list(mean, NULL)]
tmp <- unique(tmp[, names(tmp)[names(tmp) %in% tech.rep]:= NULL])
coef_biosample <- normalize(data= tmp, 
                            measure.id= 'fragment_id',
                            value.var= 'intensity',
                            rep.id= 'bio_sample',
                            group.id= NULL, 
                            return.table= FALSE, 
                            return.coef= TRUE, 
                            output.coef= biosample.coef.path, 
                            write.coef= TRUE)
setnames(coef_biosample, 'reps', 'bio_sample')

setkey(coef_biosample, bio_sample)
data <- normalize(data= data, 
                  measure.id= 'fragment_id',
                  value.var= 'intensity',
                  rep.id= 'bio_sample',
                  group.id= NULL, 
                  output.file= treated.path, 
                  coef= coef_biosample)
rm(tmp)

#super-slow 1 - divide for cv separately  2 - produce.stat subsets data table for each function separately
setkey(data, fragment_id, run_id, precursor_id)
clustering <- produce.stat(data= data, 
                           measure.id= 'fragment_id',
                           value.var= 'intensity',
                           stat.name= c('mean', 'se_corr', 'cv'), 
                           group.id= 'tech_id',
                           stat.function= list(mean, se_corr, function(x){se_corr(x)/mean(x)}))

clustering[, c('intensity', 'mean'):= list(mean, NULL)]
setkey(clustering, fragment_id, tech_id)
clustering <- unique(clustering[, names(clustering)[names(clustering) %in% ms.rep]:= NULL])
write.file(data= clustering, path= 'results/clustering.txt')

select_outcome <- function(x){
  if(x$complete==FALSE){
    return('not_complete')
  } else if(x$cluster==FALSE) {
    return('not_clustered')
  } else {
    return('ok')
  }
}

clustering[, frg_outcome:= select_outcome(unique(.SD)), by= fragment_id, .SDcols=c('complete', 'cluster')]

if(CLUSTER_PDF)
{
    pdf('results/clustering_results.pdf')
    for(precursor in unique(clustering[, precursor_id])){
        print(clustering_results(data= clustering,
                                 measure.id= 'fragment_id',
                                 group.id= 'precursor_id',
                                 value.var= 'intensity',
                                 rep.id= 'tech_id',
                                 flag= 'frg_outcome',
                                 target= precursor, 
                                 normalize= TRUE, 
                                 log= TRUE))
    }
    rm(precursor)
    dev.off()
}
    
rank_fun_sum_int <- function(data){
  if(!is.data.table(data)){
    stop('data not a DT')
  }
  sum(data[, intensity])
}

use <- data[uniq==TRUE&tryptic==TRUE&above_zero==TRUE&complete==TRUE&min_frg==TRUE&
                selected_modifications==TRUE&cluster==TRUE&min_frg_clust==TRUE]

use <- rank.groups(data= use, 
                   data.file= NULL, 
                   measure.id= 'peptide_id',
                   group.id= 'protein_id',
                   stat.function= rank_fun_sum_int)

#Precursor score = sum of fragments intensities
precursor_score <- score(data= use, 
                         data.file= NULL, 
                         measure.id= 'precursor_id',
                         value.var= 'intensity', 
                         score.function= function(x)sum(x),
                         score.name= 'precursor_score', 
                         rep.id= 'tech_id')

setkey(precursor_score, tech_id, precursor_id)
precursor_score[, intensity:= NULL]
precursor_score[, names(precursor_score)[names(precursor_score) %in% ms.rep]:= NULL]
precursor_score[, names(precursor_score)[names(precursor_score) %in% frg.lev]:= NULL]
precursor_score <- unique(precursor_score)

write.file(precursor_score, precursor.sample.score.path)

topN <- 3

#Protein score = sum of precursors scores
protein_score <- score(data= use[rank<=topN],
                       data.file= NULL, 
                       measure.id= 'protein_id',
                       value.var= 'intensity',
                       score.function= function(x)sum(x),
                       score.name= 'protein_score',
                       rep.id= 'bio_sample')

setkey(protein_score, bio_sample, protein_id)
protein_score[, intensity:= NULL]
protein_score[, names(protein_score)[names(protein_score) %in% bio.rep]:= NULL]
protein_score[, names(protein_score)[names(protein_score) %in% pep.lev]:= NULL]
protein_score <- unique(protein_score)
protein_score[, c('above_zero', 'complete', 'min_frg', 'min_frg_clust', 'rank'):= NULL]
setkey(protein_score, protein_id, bio_sample)
protein_score <- unique(protein_score)

#???
protein_score.w <- data.table::dcast.data.table(data= protein_score, formula= protein_id~bio_sample,
                                                value.var= 'protein_score')

pepnum <- use[rank<=topN, list(protein_id, peptide_id)]
pepnum[, pep_num:= ulength(peptide_id), by= protein_id]
pepnum[, peptide_id:= NULL]
pepnum <- unique(pepnum)

setkey(protein_score.w, protein_id)
setkey(protein_score, protein_id)
protein_score <- protein_score[pepnum]
protein_score.w <- protein_score.w[pepnum]

write.file(protein_score, protein.sample.score.path)
write.file(data= protein_score.w, path= 'results/protein_score_wide.txt')

protein_score[, state:= substr(bio_sample, 1+regexpr('_', bio_sample, fixed= TRUE), nchar(bio_sample)), by=list(protein_id, bio_sample)]
protein_score[, bio:= substr(bio_sample, 1, -1+regexpr('_', bio_sample, fixed= TRUE)), by=list(protein_id, bio_sample)]
protein_score[, bio_sample:= NULL]

protein_score_by_state <- data.table::dcast.data.table(data= protein_score, formula= protein_id+bio+pep_num~state, value.var= 'protein_score')
