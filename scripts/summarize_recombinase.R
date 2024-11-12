#!/usr/bin/Rscript
library(data.table)
Args<-commandArgs(trailingOnly=T)
output<-Args[1]
int_file<-Args[2]
mge_files<-Args[3:length(Args)]
int<-fread(int_file,header = F,stringsAsFactors = F,quote = '"',sep = " ")
names(int)<-c("seqid","start","end","strand","seq_len","accession","query_name","hmm_accession","E-value","score","bias","E-value2","score2","bias2","exp","reg","clu","ov","env","dom","rep","inc","description of target","seq","protein")
int<-int[,c("seqid","start","end","strand","seq_len","query_name","E-value","score","bias","protein")]
my_trim<-function(x,sep){return(paste(unlist(strsplit(x,sep))[-length(unlist(strsplit(x,sep)))],collapse=sep))}
my_trim<-Vectorize(my_trim,vectorize.args = "x")
int$seqid<-my_trim(int$seqid,"_")

my_split<-function(x){
    seq<-paste(unlist(strsplit(x,":"))[-length(unlist(strsplit(x,":")))],collapse=":")
    poss<-unlist(strsplit(unlist(strsplit(x,":"))[length(unlist(strsplit(x,":")))],"-"))
    return(c(seq,poss))
}
my_split<-Vectorize(my_split)
res<-data.frame()
for(i in mge_files){
    mges<-read.csv(i,sep="\t",header = T,stringsAsFactors=F,quote = '"')
    if (nrow(mges) > 0) {
        tmp<-my_split(mges$loc)
        mges<-data.frame(mges[,1:3],seqid=as.character(tmp[1,]),start_ins=as.integer(tmp[2,]),end_ins=as.integer(tmp[3,]),mges[,5:6])
        res<-rbind(res,merge(as.data.frame(int),mges,by="seqid"))
    }
}
write.table(res,output,sep = "\t",row.names = F,quote = F)