
library(IntAssoPlot)
###load gwas file
gwa=read.csv(file.choose(),header=T,sep="\t")
head(gwa)
gwa=gwa[,c(2,1,3,13)]
gwa[which(gwa$p==min(gwa$p,na.rm=T)),]
gwa=gwa[which(!is.na(gwa$pC)),]
names(gwa)=c("Marker","Locus","Site","p")
gwa5=gwa[which(gwa$Locus==4),]
splitfun=function(a){
  as.integer(strsplit(as.character(a[1]),":")[[1]][2])
}
gwa5$Site =apply(gwa5,1,splitfun)  

gwa5[which(gwa5$p==min(gwa5$p,na.rm=T)),]
marker=gwa5[which(gwa5$p<2e-7),c(1,2,3)]
names(marker)=c("rs","chrom","pos")
###load gtf file
gtf1=read.csv(file.choose(),header = F,sep="\t")
head(gtf)
gtf5=gtf1[which(gtf1$V1==4),]
###load hapmap
hmp=read.csv(file.choose(),header = T,sep="\t",stringsAsFactors = F)
head(hmp)
names(hmp)[1]="rs"
hmp5=hmp[which(hmp$chrom==4),]
splitfun2=function(a){
  as.integer(strsplit(as.character(a[1]),":")[[1]][2])
}
hmp5$pos=apply(hmp5,1,splitfun2)  

###compulate LD and load
#draw
IntRegionalPlot(chr=4,left=(15998883-10000),right=(15998883+10000),gtf=gtf5,association=gwa5,hapmap=hmp5,hapmap_ld=hmp5,threshold=5,leadsnp_size=2,label_gene_name = TRUE,marker2label = marker,marker2label_angle = 0,marker2label_size = 3)



