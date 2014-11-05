#
clustering_results <- function(data= NULL, measure.id= NULL, group.id= NULL, value.var= NULL, rep.id= NULL, flag=NULL, target= NULL, normalize= TRUE, log= TRUE, style= NULL){
  library(ggplot2, quietly= TRUE)
  if(!all(c(measure.id, group.id, rep.id) %in% key(data))){
    setkeyv(data, c(measure.id, group.id, rep.id))
  }
  plotMe <- unique(data[get(group.id)==target, list(get(measure.id), get(group.id), get(rep.id), get(value.var), get(flag))])
  setnames(plotMe, c('V1', 'V2', 'V3', 'V4', 'V5'), c(measure.id, group.id, rep.id, value.var, flag))
  if(normalize==TRUE){
    plotMe[, mean:= mean(get(value.var)), by= measure.id]
    if(log==TRUE){
      plotMe[, eval(value.var):= log(get(value.var)/mean, base= 2), by=list(get(measure.id), get(rep.id))]
    } else {
      plotMe[, eval(value.var):= get(value.var)/mean, by=list(get(measure.id), get(rep.id))]
    }
  }
  if(normalize==TRUE){
    if(log==TRUE){
      return(ggplot(data= plotMe, aes_string(x= rep.id, y= value.var, group= measure.id, colour= flag)) + geom_line() + ggtitle(target) + coord_cartesian(ylim = c(-3, 3)) + theme(axis.text.x= element_blank(), axis.title.x= element_blank(), axis.title.x= element_blank()) + style)
    } else {
      return(ggplot(data= plotMe, aes_string(x= rep.id, y= value.var, group= measure.id, colour= flag)) + geom_line() + ggtitle(target) + coord_cartesian(ylim = c(0, 5)) + theme(axis.text.x= element_blank(), axis.title.x= element_blank())) + style
    }
  } else {
    return(ggplot(data= plotMe, aes_string(x= rep.id, y= value.var, group= measure.id, colour= flag)) + geom_line() + ggtitle(target) + theme(axis.text.x= element_blank(), axis.title.x= element_blank())) + style
  }
}