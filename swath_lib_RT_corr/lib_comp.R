rm(list=ls())
setwd("~/swath lib RT corr")
library(package="data.table")

files<-dir(pattern="\\.txt$")
if(!length(files)==2)stop("Only two files in folder allowed: LIBRARY_*.txt and CORR_*.txt")
libfile<-files[sapply(files,FUN=function(x){substring(x,1,8)=="LIBRARY_"})]
corrfile<-files[sapply(files,FUN=function(x){substring(x,1,5)=="CORR_"})]
if(length(libfile)==0)stop("Library missing! Library name should be LIBRARY_*.txt")
if(length(corrfile)==0)stop("Correction file missing! Correction file name should be CORR_*.txt")
rm(files)
name<-paste0(substring(libfile,1,nchar(libfile)-4),"_RT_CORRECTED_BY_",substring(corrfile,6,nchar(corrfile)))

rtcorr<-fread(corrfile, sep="\t",header=TRUE,colClasses=c(uniprot_id="character",modification_sequence="character",RT_detected="numeric"),showProgress=FALSE,select=c("uniprot_id","modification_sequence","RT_detected"))
rm(corrfile)

setkey(rtcorr)
rtcorr<-unique(unique(rtcorr)[,RT_detected:=mean(RT_detected),by=list(uniprot_id,modification_sequence)])
setkey(rtcorr,uniprot_id,modification_sequence)
setnames(rtcorr,"RT_detected","RT_detected.corr")

colClasses<-c("numeric","numeric","numeric","character","boolean","numeric","character","character","numeric","character","numeric","numeric","numeric","character","numeric","boolean","numeric","numeric","boolean","numeric","numeric","character","character","character")

libinit<-fread(libfile, sep="\t",header=TRUE,colClasses=colClasses,showProgress=FALSE)
rm(list=c("colClasses","libfile"))

setkey(libinit,uniprot_id,modification_sequence)
rtinit<-unique(unique(libinit[,c("uniprot_id","modification_sequence","RT_detected"),with=FALSE],by=NULL)[,RT_detected:=mean(RT_detected),by=list(uniprot_id,modification_sequence)])
setnames(rtinit,"RT_detected","RT_detected.init")


rt<-rtinit[rtcorr,nomatch=0]
rm(rtcorr)


lmcorr<-lm(RT_detected.corr~RT_detected.init+1,data=rt)
rm(rt)


k<-lmcorr$coefficients[2]
b<-lmcorr$coefficients[1]
rm(lmcorr)

colOrder<-names(as.data.frame(libinit))
corrected<-libinit[,c("RT_detected","iRT"):=NULL][rtinit[,c("RT_detected","iRT"):=RT_detected.init*k+b][,RT_detected.init:=NULL]]
setcolorder(corrected,colOrder)

write.table(corrected,file=name,quote=FALSE,sep="\t",row.names=FALSE,na="")

rm(list=ls())