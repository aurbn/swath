rm(list=ls())
setwd("~/swath lib RT corr/")
library("data.table")
fasta<-fread("Uniprot taxonomy Escherichia coli str. K-12 substr. MG1655 keyword complete proteome 4069prot.fasta",header=FALSE)
fasta[,order:=1:nrow(fasta)]
setkey(fasta,V1)
fasta[,start:=(grepl(">",V1,fixed=TRUE))]
setkey(fasta,order)
start<-fasta$start
start2<-vector(mode="numeric",length=length(start))
start2[1]<-1
for(i in 2:length(start2)){start2[i]<-start2[i-1]+start[i]}
fasta<-cbind(fasta,belongTo=start2)
rm(list=c("i","start","start2"))
fasta[,storder:=order-min(order),by=belongTo]
names<-fasta[start==TRUE][,c("V1","start","name"):=list(NULL,NULL,substr(paste(V1," "),2,regexpr(" ",V1,fixed=TRUE)-1))]
setkey(fasta,storder)
seqs<-fasta[start==FALSE][,paste(V1,collapse=""),by=belongTo]
setkey(names,belongTo)
setkey(seqs,belongTo)
setnames(names,"name","ProteinName")
setnames(seqs,"V1","NoModSeq")
libforsort<-names[seqs][,c('belongTo','order','storder'):=NULL]
write.table(x=libforsort,file="Uniprot taxonomy Escherichia coli str. K-12 substr. MG1655 keyword complete proteome 4069prot.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
