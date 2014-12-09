#
# Supllementary functions for analysis of SWATH data 
# Version: 0.1
# Date: 18.06.2014
# Authors: Ivan Butenko, Dima Ischenko

findMaxDens <- function(x) {
  # simple function for finding maximum value from density of the data
  #  Args:
  #   x: vector with data
  #  Return: value (something like mean value) that has maximum probability from
  #    density
  
  x.d <- density(x)
  return(x.d$x[which.max(x.d$y)])
}

#Standar error of mean corrected for small samples
se_corr<-function(x,na.rm=TRUE){
  if(na.rm==TRUE){x<-x[!is.na(x)]}
  return(sd(x)*sqrt((length(x)-1)/((length(x)-1.5)*length(x))))
}

# Calculate angle between two vectors
angleDist <- function(x, y) {
  #  Args:
  #   x: first vector
  #   y: second vector
  #  Return: angle between vector
  
  acos(sum(x * y) / sqrt(sum(x ^ 2) * sum(y ^ 2)));
}


# length of unique elements in vecotor
ulength <- function(x) { length(unique(x)) }


scaleNorm<-function(dt, #data.table - data in 'long' format
                    measure.id, #name of column with measured feature id (fragment id)
                    value, #value assigned to feature (fragment AUC)
                    rep.id, #name of column with id's of sets of repeated measures (run id)
                     
                    group.id=NULL #name of column of groups to be scaled separately
)
{
  require(data.table, quietly=TRUE)
  require(reshape2, quietly=TRUE)

  if(!is.null(group.id)){
    setkeyv(dt,group.id)
    
    groups<-as.character(unique(unlist(dt[,group.id,with=FALSE])))
    if(length(groups)==0) stop("No groups")
    l<-lapply(groups,FUN=function(x)
      data.table(group=x,
                 scaleNorm(dt[x],measure.id,value,rep.id))
    )
    
    return(setnames(rbindlist(l),c(group.id,rep.id,"coef")))         
  } else {
    
    setkeyv(dt,measure.id)
    
    reps<-as.character(unique(unlist(dt[,rep.id,with=FALSE])))
    if(length(reps)==0) stop("No repeats")
    if(length(reps)==1)return(data.table(reps=reps[1],coef=1))
    
    dt.wide<-dcast.data.table(dt,formula(sprintf("%s ~ %s", measure.id, rep.id)),value.var=value)
    
    l<-lapply(reps[2:length(reps)],FUN=function(x)
      data.table(reps=x,coef=
                   10^(-findMaxDens(log(unlist(
                     dt.wide[!(get(x)==0|get(reps[1])==0),x,with=FALSE]/dt.wide[!(get(x)==0|get(reps[1])==0),reps[1],with=FALSE]
                     ,use.names=FALSE),base=10)
                   ))
      ))
 
    l[[length(l)+1]]<-data.table(reps=reps[1],coef=1)
    return(rbindlist(l))
  }
}

# Test data for scaleNorm
#  test<-data.table(group=c(rep("a",times=20),rep("b",times=20)),rep=rep(c("a1","a2","b1","b2"),each=10),meas=rep(c("z","x","c","v","g","n","m","s","d","f"),times=4),value=c(2000+rnorm(10),1000+rnorm(10),5000+rnorm(10),1000+rnorm(10)))
#  test2<-test[1:20]
#  {dt=test2;measure.id="meas";value="value";rep.id="rep"}
#  group.id=NULL



#Check and create directories

check.directory <- function(mainDir= getwd(), subDir= 'results', silent= TRUE){
  
  dir <- paste(mainDir, subDir, sep = "/", collapse = "/")
  
  if (file.exists(paste(mainDir, subDir, "/", sep = "/", collapse = "/"))) {
    if(!silent){cat(paste0(dir, " exists and is a directory.\n"))}
    return(0)
  } else if (file.exists(paste(mainDir, subDir, sep = "/", collapse = "/"))) {
    if(!silent){cat(paste0(dir, " exists but is a file.\n"))}
    return(1)
    # you will probably want to handle this separately
  } else {
    if(!silent){cat(paste0(dir, " does not exist.\n"))}
    return(2)
  }
}

create.directory <- function(mainDir= getwd(), subDir= 'results', silent= TRUE){

check <- check.directory(mainDir, subDir, silent)

if(check==0|check==1){
  return(check)
}

dir <- paste(mainDir, subDir, sep = "/", collapse = "/")

if(!silent){cat(paste0('Creating ', dir, '\n'))}
dir.create(file.path(mainDir, subDir))

check <- check.directory(mainDir, subDir, silent=TRUE)

if(check==0){
  if(!silent){cat(paste0(dir, ' created.\n'))}
  return(check)
} else {
  if(!silent){cat(paste0('Failed to create ', dir, '\n'))}
  return(3)
}
}

dt.read <- function(path= NULL, ...){
  return(fread(input= path, sep= '\t', header= TRUE, showProgress= FALSE, ...))
}

write.file <- function(data= NULL, path= NULL, ...){
  write.table(x= data, file= path, append= FALSE, quote= FALSE, sep= '\t', eol= '\r\n', na= 'NA', dec= '.', row.names= FALSE, col.names= TRUE, qmethod= 'escape')
}