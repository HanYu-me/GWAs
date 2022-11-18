library(GenABEL)
library(parallel)
###data save path
#run data
path1=""
#draw figure
path2=""
#software path
gcta=""
plink=""
gec=""

###add phe data to GenABEL data
#load the 1001 arabidopsis data
load("data.RData")
#load the phe data
phe=as.data.frame(read.csv(""))
fl.data.merge=add.phdata(fl.data,Add_phe)
phdata(fl.data.merge)[1:10,]
nids(fl.data.merge)

###Quality Control
#QC to all individual
qc1 <- check.marker(fl.data.merge, p.level=0)
fl.data.merge=fl.data.merge[qc1$idok,qc1$snpok]
ibs_m=ibs(fl.data.merge,weight="freq")
dis_m=as.dist(0.5-ibs_m)
cmd=cmdscale(dis_m)
#choose best k
wss <- (nrow(cmd)-1)*sum(apply(cmd,2,var))
for(i in 2:10){
  wss[i]=kmeans(cmd,centers=i)$withinss  
}
plot(1:10, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
i=2
km=kmeans(cmd,centers=i)
plot(cmd[,1],cmd[,2],main=paste("k=",as.character(i)))
points(cmd[which(km$cluster==1),1],cmd[which(km$cluster==1),2],col="red")
points(cmd[which(km$cluster==2),1],cmd[which(km$cluster==2),2],col="blue")
i=3
km=kmeans(cmd,centers=i)
plot(cmd[,1],cmd[,2],main=paste("k=",as.character(i)))
points(cmd[which(km$cluster==1),1],cmd[which(km$cluster==1),2],col="red")
points(cmd[which(km$cluster==2),1],cmd[which(km$cluster==2),2],col="blue")
points(cmd[which(km$cluster==3),1],cmd[which(km$cluster==3),2],col="orange")
i=4
km=kmeans(cmd,centers=i)
plot(cmd[,1],cmd[,2],main=paste("k=",as.character(i)))
points(cmd[which(km$cluster==1),1],cmd[which(km$cluster==1),2],col="red")
points(cmd[which(km$cluster==2),1],cmd[which(km$cluster==2),2],col="blue")
points(cmd[which(km$cluster==3),1],cmd[which(km$cluster==3),2],col="orange")


#choose the smallest group as outer group
sub_id2=names(km$cluster[which(km$cluster==1)])
sub_id1=names(km$cluster[which(km$cluster!=1)])
#choose main group
fl_1=fl.data.merge[which(fl.data.merge@phdata$id%in% sub_id1),]
qc2 <- check.marker(fl_1, p.level=0)
fl_1=fl_1[qc2$idok,qc2$snpok]
nids(fl_1)
nids(fl.data.merge)

###heritability and GWAs
#fun for relationship matrix
GCTA_estKin <- function(gcta,output,input){
  cmd_k <- paste(gcta,"--bfile",input, "--make-grm" ,"--make-grm-alg 1","--out",output,sep=" ")
  system(cmd_k)
}
#fun for heritability
GCTA_REML <- function(t1,output,gcta,phenotype,grm){
  cmd_biGREML <- paste(gcta,"--pheno",phenotype,"--grm",grm,"--reml","--mpheno",t1, "--out",output,sep=" ")
  cat(cmd_biGREML,"\n")
  system(cmd_biGREML)
  #print(paste(output,"is finished"))
}
#fun for mlm
GCTA_mlm <- function(gcta,bfile,grm,pheno,mlm,mark=0){
  print(mark)
  if(mark==1){
    #gcta_v1.94.0Beta_windows_x86_64 --bfile CD_B --make-grm --sparse-cutoff 0.05 --thread-num 10 --out sp_grm
    cmd_grm=paste(gcta,"--bfile",bfile,"--make-grm --sparse-cutoff 0.05 --thread-num 10 --out",grm,sep=" ")
    cat(cmd_grm,"\n")
    system(cmd_grm)  
  }
  #	gcta_v1.94.0Beta_windows_x86_64 --bfile CD_B --grm-sparse sp_grm --fastGWA-mlm --pheno CD_B_phe.txt --thread-num 10 --out geno_assoc
  cmd_mlm=paste(gcta,"--bfile",bfile,"--grm-sparse",grm,"--fastGWA-mlm --pheno",pheno,"--thread-num 10 --out",mlm)
  cat(cmd_mlm,"\n")
  system(cmd_mlm)
}
export.plink(fl_1)
#convert to binary file
cmd_plink=paste(plink,"--tfile C:\\Users\\MSI\\OneDrive\\Documents\\plink --make-bed --out","",sep=" ")
system(cmd_plink)
#Kinship
GCTA_estKin(gcta=gcta,output = paste(path1,"oak",sep=""),input=paste(path1,"out_","",sep=""))
#write phenotype data
phe_file=phdata(fl_2)[,c(1,2,3)]
phe_file[,2]=as.numeric(phe_file[,1])
phe_file[,1]=1:nrow(phe_file)
phe_file[,1]=as.numeric(phe_file[,1])
file_path=paste(path1,names(phdata(fl_2))[i],"_phe.txt",sep="")
write.table(phe_file,file = file_path,col.names = F,row.names = F)
#heritability
GCTA_REML(t1=1,output=paste(path1,names(phdata(fl_1))[3],sep=""),gcta=gcta,phenotype=paste(path1,names(phdata(fl_1))[3],"_phe.txt",sep=""),grm=paste(path1,"oak",sep=""))
#GWAs
GCTA_mlm(gcta=gcta,bfile=paste(path1,"out_",phname[i],sep=""),grm=paste(path1,"mlm\\",names(phdata(fl_1))[3],"sp_oak",sep=""),pheno=paste(path1,names(phdata(fl_1))[3],"_phe.txt",sep=""),mlm=paste(path1,"mlm//",names(phdata(fl_1))[3],"_mlm",sep=""),mark=1)

###use qqman to visulization
data=read.csv(file=paste(path1,"mlm\\",names(phdata(fl_1))[3],"_mlm.fastGWA",sep=""),sep="\t",header=T)
data=data[!is.na(data$P),]
data=data[,c(1,2,3,10)]
names(data)=c("CHR","SNP","BP","P")
lam <- estlambda(data$P,plot = F) 
chis <- qchisq(data$P,df=1,lower.tail = FALSE)
data$P <- pchisq(chis/lam$estimate,df=1, lower.tail = FALSE)
png(paste(path2,"raw_data\\qq_",phname[i],".png",sep=""))
lam <- estlambda(data$P,plot = T) 
text(5, 15, paste("lambda =",substr(as.character(lam$estimate),1,6)))
dev.off()
png(paste(path2,"raw_data\\",phname[i],".png",sep=""))
qqman::manhattan(data,col = c("blue4", "orange3"))
dev.off()
