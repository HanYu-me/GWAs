library(GenABEL)
library(readxl)
library(parallel)
###data save path
#run data
path1="C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\工作数据\\2022.09.27\\"
#draw figure
path2="C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.09.27\\"
#software path
gcta <- "C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\gcta_v1.94.0Beta_windows_x86_64"
plink <- "C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\plink_win64_20220402\\plink"
gec="java -jar C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.09.05\\gec\\gec.jar"

###processing raw phenotype data
#load the 1001 arabidopsis data
load("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\工作数据\\2022.07.11\\fl.data.RData")
#load the phe data
phe=as.data.frame(read_excel("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\工作数据\\2022.07.27\\成都和红原同质园生态型拟南芥开花物候数据.xlsx",sheet=1,col_names=T))
for(i in 12:ncol(phe)){
  phe[,i]=as.numeric(phe[,i])
}
table(phe$Commongarden)
CD_phe=phe[which(phe$Commongarden=="Chengdu"),]
HY_phe=phe[which(phe$Commongarden=="Hongyuan"),]
#load id
arab_name=as.data.frame(read_excel("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\工作数据\\2022.07.27\\生态型拟南芥种源地信息.xlsx",sheet=1,col_names=T))
#caculate CV(coeffcient variation)
CD_CV=data.frame(matrix(nrow=nrow(CD_phe)/12,ncol=8))
names(CD_CV)=names(phe)[12:19]
HY_CV=data.frame(matrix(nrow=nrow(HY_phe)/8,ncol=8))
names(HY_CV)=names(phe)[12:19]
#CD CV
for(i in seq(from=1,to=nrow(CD_phe),by=12)){
  for(j in 12:19){
    B=CD_phe[i:(i+11),j]  
    B=sd(B,na.rm = TRUE)/mean(B,na.rm = TRUE)
    CD_CV[(i%/%12)+1,(j-11)]=B
  }
}
#HY CV
for(i in seq(from=1,to=nrow(HY_phe),by=8)){
  for(j in 12:19){
    B=HY_phe[i:(i+7),j]
    B=sd(B,na.rm = TRUE)/mean(B,na.rm = TRUE)
    HY_CV[(i%/%8)+1,(j-11)]=B
  }
}
#save raw CV
HY_CV2=HY_CV
CD_CV2=CD_CV
#caculate the mean ph data (if cv >0.05 or cv>95%CV, remove the maxmal and minmal 2 indivadual)
Add_phe=data.frame(matrix(nrow=nrow(HY_phe)/8,ncol=17))
names(Add_phe)=c("id",c(paste(rep("CD",8),names(phe)[12:19],sep="_")),c(paste(rep("HY",8),names(phe)[12:19],sep="_")))
sort(table(arab_name$id_1001))
for(i in 1:240){
  Add_phe[i,1]=as.numeric(as.character(arab_name[which(arab_name$label==i),9]))
  for(j in 12:19){
    CD=CD_phe[((i-1)*12+1):(i*12),j]
    if(!is.na(sd(CD,na.rm = TRUE)/mean(CD,na.rm = TRUE))){
      if(mean(CD_CV2[,(j-11)],na.rm=T)<0.05){
        if(sd(CD,na.rm = TRUE)/mean(CD,na.rm = TRUE)>quantile(sort(CD_CV2[,(j-11)]),0.95)){
          CD=sort(CD)
          CD=CD[3:(length(CD)-2)]
          CD_CV[i,(j-11)]=sd(CD,na.rm=T)/mean(CD,na.rm = T)
          CD=mean(CD,na.rm=T)
        }else{
          CD=mean(CD,na.rm=T)
        }
      }else if(mean(CD_CV2[,(j-11)],na.rm=T)>=0.05){
        if(sd(CD,na.rm = TRUE)/mean(CD,na.rm = TRUE)>0.05){
          CD=sort(CD)
          CD=CD[3:(length(CD)-2)]
          CD_CV[i,(j-11)]=sd(CD,na.rm=T)/mean(CD,na.rm = T)
          CD=mean(CD,na.rm=T)
        }else{
          CD=mean(CD,na.rm=T)
        }
      }
      
    }else{
      CD=NA
    }
    Add_phe[i,(j-10)]=CD
  }
  
  for(j in 12:19){
    CD=HY_phe[((i-1)*8+1):(i*8),j]
    if(!is.na(sd(CD,na.rm = TRUE)/mean(CD,na.rm = TRUE))){
      if(mean(HY_CV2[,(j-11)],na.rm=T)<0.05){
        if(sd(CD,na.rm = TRUE)/mean(CD,na.rm = TRUE)>quantile(sort(HY_CV2[,(j-11)]),0.95)){
          CD=sort(CD)
          CD=CD[3:(length(CD)-2)]
          HY_CV[i,(j-11)]=sd(CD,na.rm=T)/mean(CD,na.rm = T)
          CD=mean(CD,na.rm=T)
        }else{
          CD=mean(CD,na.rm=T)
        }
      }else if(mean(HY_CV2[,(j-11)],na.rm=T)>=0.05){
        if(sd(CD,na.rm = TRUE)/mean(CD,na.rm = TRUE)>0.05){
          CD=sort(CD)
          CD=CD[3:(length(CD)-2)]
          HY_CV[i,(j-11)]=sd(CD,na.rm=T)/mean(CD,na.rm = T)
          CD=mean(CD,na.rm=T)
        }else{
          CD=mean(CD,na.rm=T)
        }
      }
      
    }else{
      CD=NA
    }
    Add_phe[i,(j-2)]=CD
  }
  
}
Add_phe=Add_phe[which(!is.na(Add_phe$id)),]
                
###CV visualization
CD_CV$id=0
CD_CV2$id=0
HY_CV$id=0
HY_CV2$id=0
for(i in 1:240){
  CD_CV$id[i]=as.numeric(as.character(arab_name[which(arab_name$label==i),9]))
  CD_CV2$id[i]=as.numeric(as.character(arab_name[which(arab_name$label==i),9]))
  HY_CV$id[i]=as.numeric(as.character(arab_name[which(arab_name$label==i),9]))
  HY_CV2$id[i]=as.numeric(as.character(arab_name[which(arab_name$label==i),9]))
}
CD_CV_main=CD_CV[which(CD_CV$id%in%phdata(fl_1)[,1]),]
CD_CV2_main=CD_CV2[which(CD_CV2$id%in%phdata(fl_1)[,1]),]
CD_CV_all=CD_CV[which(CD_CV$id%in%phdata(fl.data.merge)[,1]),]
CD_CV2_all=CD_CV2[which(CD_CV2$id%in%phdata(fl.data.merge)[,1]),]

HY_CV_main=HY_CV[which(HY_CV$id%in%phdata(fl_1)[,1]),]
HY_CV2_main=HY_CV2[which(HY_CV2$id%in%phdata(fl_1)[,1]),]
HY_CV_all=HY_CV[which(HY_CV$id%in%phdata(fl.data.merge)[,1]),]
HY_CV2_all=HY_CV2[which(HY_CV2$id%in%phdata(fl.data.merge)[,1]),]

'
for(i in 1:8){
  png(paste(path,"CD_",names(CD_CV)[i],"_CV.png",sep=""))
  par(mfrow=c(2,2))
  hist(CD_CV2_main[,i],breaks=90,main=paste("Main group raw CV of CD",names(CD_CV)[i]),xlab = "value")
  hist(CD_CV2_all[,i],breaks=90,main=paste("All group raw CV of CD",names(CD_CV)[i]),xlab = "value")
  hist(CD_CV_main[,i],breaks=90,main=paste("Main group adjusted CV of CD",names(CD_CV2)[i]),xlab = "value")  
  hist(CD_CV_all[,i],breaks=90,main=paste("All group adjusted CV of CD",names(CD_CV2)[i]),xlab = "value")  
  dev.off()
}
for(i in 1:8){
  png(paste(path,"HY_",names(HY_CV)[i],"_CV.png",sep=""))
  par(mfrow=c(2,2))
  hist(HY_CV2_main[,i],breaks=90,main=paste("Main group raw CV of HY",names(HY_CV)[i]),xlab = "value")
  hist(HY_CV2_all[,i],breaks=90,main=paste("All group raw CV of HY",names(HY_CV)[i]),xlab = "value")
  hist(HY_CV_main[,i],breaks=90,main=paste("Main group adjusted CV of HY",names(HY_CV2)[i]),xlab = "value")  
  hist(HY_CV_all[,i],breaks=90,main=paste("All group adjusted CV of HY",names(HY_CV2)[i]),xlab = "value")  
  dev.off()
}
'


#drop duplicated data
Add_phe_dup <- Add_phe[which(Add_phe$id %in% Add_phe[duplicated(Add_phe$id),1]),]
Add_phe=Add_phe[which(!Add_phe$id%in%Add_phe_dup$id),]
sort(table(Add_phe_dup$id))
sort(table(Add_phe$id))
#input location data to compare with id data
data2=read.csv("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\工作数据\\2022.07.28\\1001genomes-accessions.csv",header=T,sep=";")
data2=data2[,c(1,3,6,7)]
data=arab_name[which(arab_name$id_1001 %in% Add_phe_dup$id),c(9,2,4,5)]
names(data)[1]="id"
data3=merge.data.frame(data,data2,by="id")
data3=data3[,c(1,2,3,6,4,7,5)]
k=data3[which((data3$latitude.x!=data3$latitude.y)|(data3$longitude.x!=data3$longitude.y)),]
unique(k$id)
#droup all different data
Add_phe_dup=Add_phe_dup[which(!Add_phe_dup$id %in% unique(k$id)),]
name=unique(Add_phe_dup$id)
for( i in 1:length(name)){
  a=Add_phe_dup[which(Add_phe_dup$id==name[i]),]
  print(a)
  for(j in 1:ncol(a)){
    a[1,j]=median(a[,j],na.rm=T)
  }
  Add_phe=rbind.data.frame(Add_phe,a[1,])
}
##check data 
Add_phe$CD_f_lfdays=Add_phe$CD_lf_doy-Add_phe$CD_ff_doy
Add_phe$CD_b_ffdays=Add_phe$CD_ff_doy-Add_phe$CD_bol_doy
Add_phe$CD_b_lfdays=Add_phe$CD_lf_doy-Add_phe$CD_bol_doy
Add_phe$CD_s_bdays=Add_phe$CD_bol_doy-25
Add_phe$CD_s_ffdays=Add_phe$CD_ff_doy-25
Add_phe$CD_s_lfdays=Add_phe$CD_lf_doy-25

Add_phe$HY_f_lfdays=Add_phe$HY_lf_doy-Add_phe$HY_ff_doy
Add_phe$HY_b_ffdays=Add_phe$HY_ff_doy-Add_phe$HY_bol_doy
Add_phe$HY_b_lfdays=Add_phe$HY_lf_doy-Add_phe$HY_bol_doy
Add_phe$HY_s_bdays=Add_phe$HY_bol_doy-145
Add_phe$HY_s_ffdays=Add_phe$HY_ff_doy-145
Add_phe$HY_s_lfdays=Add_phe$HY_lf_doy-145

#add the subtract trait 
names(Add_phe)
Add_phe=Add_phe[,c(1:9,18,10:17,19)]
names(Add_phe)
for(i in 1:9){
  Add_phe=cbind.data.frame(Add_phe,(Add_phe[,(i+10)]-Add_phe[,(i+1)]))
  names(Add_phe)[19+i]=paste(names(Add_phe)[i+10],"_",names(Add_phe)[i+1],sep="")
}

##add phe data to GenABEL data
fl.data.merge=add.phdata(fl.data,Add_phe)
phdata(fl.data.merge)[1:10,]


#del empty phe individual
names(phdata(fl.data.merge)[1,])
fl.data.merge=del.phdata(fl.data.merge,c("FT10_mean","FT16_mean","FT10_sd","FT16_sd"))
names(phdata(fl.data.merge)[1,])
fl.data.merge=fl.data.merge[which(fl.data.merge@phdata$id%in%Add_phe$id),]
nids(fl.data.merge)

###Quality Control
#QC to all individual
qc1 <- check.marker(fl.data.merge, p.level=0)
fl.data.merge=fl.data.merge[qc1$idok,qc1$snpok]
#export to caculate admixture
#export.plink(fl.data.merge)
#population structure and QC to main group
par(mfrow=c(2,2))
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

#use LEA to run analysis of population structure
#install.packages(c("fields","RColorBrewer","mapplots"))
#BiocManager::install("LEA")
library(LEA)
# genabel to bed to ped
#plink.exe --tfile plink --make-bed --out out_all
#plink.exe --bfile ut_all --recode --tab --out out_all
memory.limit(1000000)
output = ped2lfmm("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\plink_win64_20220402\\out_all.ped","C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\工作数据\\2022.09.26\\out_all.lfmm")
project = snmf(paste(path1,"out_all.lfmm",sep=""),K = 1:10,entropy = TRUE,repetitions = 1,project = "new")
project = load.snmfProject(paste(path1,"out_all.snmfProject",sep=""))
plot(project, col = "blue", pch = 19, cex = 1.2)
#choose the lowest entropy K
obj.snmf=snmf(paste(path1,"out_all.lfmm",sep=""), K = 5,entropy = TRUE, repetitions = 1,project = "new")
qmatrix = Q(obj.snmf, K = 5)
qmatrix=as.data.frame(qmatrix)
head(qmatrix)
#visulise
qmatrix=qmatrix[order(qmatrix$V1,qmatrix$V2),]
barplot(t(qmatrix), col = rainbow(3), border = NA, space = 0,xlab = "Individuals", ylab = "Admixture coefficients")



#choose the smallest group as outer group
sub_id2=names(km$cluster[which(km$cluster==1)])
sub_id1=names(km$cluster[which(km$cluster!=1)])
#choose main group
fl_1=fl.data.merge[which(fl.data.merge@phdata$id%in% sub_id1),]
qc2 <- check.marker(fl_1, p.level=0)
fl_1=fl_1[qc2$idok,qc2$snpok]
save(fl_1,file=paste(path1,"fl_1.RData",sep=""))
#load(file=paste(path1,"fl_1.RData",sep=""))
save(fl.data.merge,file=paste(path1,"fl.data.merge.RData",sep=""))
#load(file=paste(path1,"fl.data.merge.RData",sep=""))
nids(fl_1)
nids(fl.data.merge)

###heritability and GWAs

CD_CV$id=1:240
CD_CV$id=as.numeric(as.character(arab_name[which(arab_name$label==CD_CV$id),9]))
CD_CV=CD_CV[,c(9,1:8)]
All_CV=cbind.data.frame(CD_CV,HY_CV)
All_CV=cbind.data.frame(cbind.data.frame(All_CV[,c(1:9)],c(rep(0,nrow(All_CV)))),All_CV[,c(10:18)])

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

for(i in 3:29){
  if((i<19)&(i!=11)){
    print(paste(names(phdata(fl_1))[i],names(All_CV)[i-1]))
  }
}
#caculating main group ...
for(i in 3:29){
  print(i)
  if((i<19)&(i!=11)){
    fl_2=fl_1[which(!fl_1@phdata$id%in%All_CV[which(All_CV[,(i-1)]>quantile(All_CV[,(i-1)],0.95,na.rm = T)),1]),]  
    export.plink(fl_2)
    cmd_plink=paste(plink,"--tfile C:\\Users\\MSI\\OneDrive\\Documents\\plink --make-bed --out",paste(path1,"out_",names(phdata(fl_1))[i],sep=""),sep=" ")
    system(cmd_plink)
  }else{
    export.plink(fl_1)
    fl_2=fl_1
    cmd_plink=paste(plink,"--tfile C:\\Users\\MSI\\OneDrive\\Documents\\plink --make-bed --out",paste(path1,"out_",names(phdata(fl_1))[i],sep=""),sep=" ")
    system(cmd_plink)
  }
  GCTA_estKin(gcta=gcta,output = paste(path1,"oak",sep=""),input=paste(path1,"out_",names(phdata(fl_1))[i],sep=""))
  phe_file=phdata(fl_2)[,c(1,2,i)]
  phe_file[,2]=as.numeric(phe_file[,1])
  phe_file[,1]=1:nrow(phe_file)
  phe_file[,1]=as.numeric(phe_file[,1])
  file_path=paste(path1,names(phdata(fl_2))[i],"_phe.txt",sep="")
  write.table(phe_file,file = file_path,col.names = F,row.names = F)
  GCTA_REML(t1=1,output=paste(path1,"heritability\\",names(phdata(fl_1))[i],sep=""),gcta=gcta,phenotype=paste(path1,names(phdata(fl_1))[i],"_phe.txt",sep=""),grm=paste(path1,"oak",sep=""))
}

#caculating all group
for(i in 3:29){
  print(i)
  if((i<19)&(i!=11)){
    fl_2=fl.data.merge[which(!fl.data.merge@phdata$id%in%All_CV[which(All_CV[,(i-1)]>quantile(All_CV[,(i-1)],0.95,na.rm = T)),1]),]  
    export.plink(fl_2)
    cmd_plink=paste(plink,"--tfile C:\\Users\\MSI\\OneDrive\\Documents\\plink --make-bed --out",paste(path1,"out_all_",names(phdata(fl_1))[i],sep=""),sep=" ")
    system(cmd_plink)
  }else{
    export.plink(fl.data.merge)
    fl_2=fl.data.merge
    cmd_plink=paste(plink,"--tfile C:\\Users\\MSI\\OneDrive\\Documents\\plink --make-bed --out",paste(path1,"out_all_",names(phdata(fl_1))[i],sep=""),sep=" ")
    system(cmd_plink)
  }
  GCTA_estKin(gcta=gcta,output =paste(path1,"oak2",sep=""),input=paste(path1,"out_all_",names(phdata(fl_1))[i],sep=""))
  phe_file=phdata(fl_2)[,c(1,2,i)]
  phe_file[,2]=as.numeric(phe_file[,1])
  phe_file[,1]=1:nrow(phe_file)
  phe_file[,1]=as.numeric(phe_file[,1])
  file_path=paste(path1,names(phdata(fl_1))[i],"_phe_all.txt",sep="")
  write.table(phe_file,file = file_path,col.names = F,row.names = F)
  GCTA_REML(t1=1,output=paste(path1,"heritability\\",names(phdata(fl_1))[i],"all",sep=""),gcta=gcta,phenotype=paste(path1,names(phdata(fl_1))[i],"_phe_all.txt",sep=""),grm=paste(path1,"oak2",sep=""))
}

#summary heritability
herita=data.frame(matrix(nrow=27,ncol=6))
names(herita)=c("main_group_h2","se","p","all_h2","se2","p2")
row.names(herita)=names(phdata(fl_1))[3:29]
for(i in 3:29){
  path=paste(path1,"heritability\\",names(phdata(fl_1))[i],".hsq",sep="")
  her_trait=read.csv(file=path,stringsAsFactors = F,sep="\t")
  herita[(i-2),1:3]=as.numeric(c(her_trait[4,2:3],her_trait[9,2]))
  path=paste(path1,"heritability\\",names(phdata(fl_1))[i],"all.hsq",sep="")
  her_trait=read.csv(file=path,stringsAsFactors = F,sep="\t")
  herita[(i-2),4:6]=as.numeric(c(her_trait[4,2:3],her_trait[9,2]))
}

##visulization boxplot
library(ggplot2)
herita$name=row.names(herita)
##significant heritability
herita2=herita[which(herita$p<0.05),]
herita3=herita[which(herita$p2<0.05),]
#main
ggplot(herita2,aes(x=reorder(name,main_group_h2),y=main_group_h2,fill=name),ylab="",xlab="")+
  geom_bar(stat="identity",width=0.5,color="black")+
  geom_errorbar(aes(ymin=main_group_h2-se, ymax=main_group_h2+se), width=0,)+
  geom_text(aes(label = paste("p value =",p)),size=3,hjust = 1,vjust=3)+
  guides(fill=FALSE)+
  theme_bw()+
  coord_flip(ylim = c(0,1))+
  xlab("trait")+ylab("Heritability")+ggtitle("Heritability of main gorup")
ggsave(file=paste(path2,"heritability_main_significant.png",sep=""),dpi=300,width=6,height = 8)
#all
ggplot(herita,aes(x=reorder(name,all_h2),y=all_h2,fill=name),ylab="",xlab="")+
  geom_bar(stat="identity",width=0.5,color="black")+
  geom_errorbar(aes(ymin=all_h2-se2, ymax=all_h2+se2), width=0,)+
  geom_text(aes(label = paste("p value =",p2)),size=3,hjust = 1,vjust=2)+
  guides(fill=FALSE)+
  theme_bw()+
  coord_flip(ylim = c(0,1))+
  xlab("trait")+ylab("Heritability")+ggtitle("Heritability of all")
ggsave(file=paste(path2,"heritability_all.png",sep=""),dpi=300,width=6,height = 8)
###parallel mlm
phname<-names(phdata(fl_1))
funmlm <- function(i) {
  GCTA_mlm(gcta=gcta,bfile=paste(path1,"out_",phname[i],sep=""),grm=paste(path1,"mlm\\",phname[i],"sp_oak",sep=""),pheno=paste(path1,phname[i],"_phe.txt",sep=""),mlm=paste(path1,"mlm//",phname[i],"_mlm",sep=""),mark=1)
  GCTA_mlm(gcta=gcta,bfile=paste(path1,"out_all_",phname[i],sep=""),grm=paste(path1,"mlm\\",phname[i],"sp_oak2",sep=""),pheno=paste(path1,phname[i],"_phe_all.txt",sep=""),mlm=paste(path1,"mlm//",phname[i],"_mlm_all",sep=""),mark=1)  
}
clnum<-detectCores() #set parallel cores
cl <- makeCluster(getOption("cl.cores", clnum));
clusterExport(cl,deparse(substitute(funmlm)))
clusterEvalQ(cl,library(GenABEL))#add library into parallel enviroment
clusterExport(cl,c("phname","GCTA_mlm","gcta","path1"))#add variable into parrllel enviroment
parLapply(cl, 3:29,  funmlm)
stopCluster(cl)

###use gec to estimate effectional snp number

funforpar <- function(i){
  cmd <- paste(gec, "Xmxlg --effect-number --plink-binary",paste(path1,"out_",phname[i],sep=""), "--genome --out",paste(path1,"gec\\","gec_",phname[i],sep=""))
  system(cmd)
  print(i)
  cmd <- paste(gec, "Xmxlg --effect-number --plink-binary",paste(path1,"out_all_",phname[i],sep=""), "--genome --out",paste(path1,"gec\\","gec_all_",phname[i],sep=""))
  system(cmd)
  print(i)
}
#parallel computing
clnum<-detectCores() 
cl <- makeCluster(getOption("cl.cores", clnum));
clusterExport(cl,deparse(substitute(funforpar)))
clusterExport(cl,c("phname","path1","gec"))
parLapply(cl, 3:29,  funforpar)
stopCluster(cl)


#summary snp number
phname<-names(phdata(fl_1))
effsnp=data.frame(matrix(nrow=27,ncol=13))
i=3
data=read.csv(file=paste(path1,"gec\\","gec_all_",phname[i],".sum",sep=""),sep="\t")
names(effsnp)=c("trait",names(data),paste("all",names(data)))
for(i in 3:29){
  effsnp[(i-2),1]=phname[i]
  data=read.csv(file=paste(path1,"gec\\","gec_",phname[i],".sum",sep=""),sep="\t")
  effsnp[(i-2),2:7]=data[1,]
  data=read.csv(file=paste(path1,"gec\\","gec_all_",phname[i],".sum",sep=""),sep="\t")
  effsnp[(i-2),8:13]=data[1,]
  
}

###use cojo to filter significant snp
funcojo<-function(i){
  data=read.csv(file=paste(path1,"mlm\\",phname[i],"_mlm.fastGWA",sep=""),sep="\t",header=T)
  data=data[,c(2,4,5,7,8,9,10,6)]
  lam <- estlambda(data$P,plot = F) 
  chis <- qchisq(data$P,df=1,lower.tail = FALSE)
  ps <- pchisq(chis/lam$estimate,df=1, lower.tail = FALSE)
  data$P=ps
  names(data)=c("SNP","A1","A2","freq","b","se","p","N")
  write.table(data,file=paste(path1,"cojo\\","cojo_",phname[i],".ma",sep=""),sep="\t",row.names = F,quote = FALSE)
  cmd <- paste(gcta, "--bfile",paste(path1,"out_",phname[i],sep=""), "--cojo-file",paste(path1,"cojo\\","cojo_",phname[i],".ma",sep=""), "--cojo-p",as.character(effsnp[(i-2),6]*2), "--cojo-collinear 0.1", "--cojo-wind 200", "--cojo-slct","--maf 0.05", "--out", paste(path1,"cojo\\","cojo_result_",phname[i],sep=""))
  system(cmd)
  
  data=read.csv(file=paste(path1,"mlm\\",phname[i],"_mlm_all.fastGWA",sep=""),sep="\t",header=T)
  data=data[,c(2,4,5,7,8,9,10,6)]
  lam <- estlambda(data$P,plot = F) 
  chis <- qchisq(data$P,df=1,lower.tail = FALSE)
  ps <- pchisq(chis/lam$estimate,df=1, lower.tail = FALSE)
  data$P=ps
  names(data)=c("SNP","A1","A2","freq","b","se","p","N")
  write.table(data,file=paste(path1,"cojo\\","cojo_",phname[i],"_all.ma",sep=""),sep="\t",row.names = F,quote = FALSE)
  cmd <- paste(gcta, "--bfile",paste(path1,"out_all_",phname[i],sep=""), "--cojo-file",paste(path1,"cojo\\","cojo_",phname[i],"_all.ma",sep=""), "--cojo-p",as.character(effsnp[(i-2),12]*2), "--cojo-collinear 0.1", "--cojo-wind 200", "--cojo-slct","--maf 0.05", "--out", paste(path1,"cojo\\","cojo_result_all_",phname[i],sep=""))
  system(cmd)
}
clnum<-detectCores() 
cl <- makeCluster(getOption("cl.cores", clnum));
clusterExport(cl,deparse(substitute(funcojo)))
clusterEvalQ(cl,library(GenABEL))
clusterExport(cl,c("phname","effsnp","path1","gcta"))
parLapply(cl, 3:29,  funcojo)
stopCluster(cl)
##check cojo result number
#get file name in path
result_name=list.files(paste(path1,"cojo\\",sep=""))
result_name_match=data.frame(matrix(nrow=27,ncol=4))
result_name_match[,1]=paste(rep("cojo_result_",27),phname[3:29],rep(".cma.cojo"),sep="")
result_name_match[,3]=paste(rep("cojo_result_all_",27),phname[3:29],rep(".cma.cojo"),sep="")
for(i in 1:27){
  name=result_name_match[i,1]
  if(name%in%result_name){
    result_name_match[i,2]=1
  }
  name=result_name_match[i,3]
  if(name%in%result_name){
    result_name_match[i,4]=1
  }
}

###use qqman to visulization
##row data
draw_manha_raw<-function(i){
  data=read.csv(file=paste(path1,"mlm\\",phname[i],"_mlm.fastGWA",sep=""),sep="\t",header=T)
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
  
  data=read.csv(file=paste(path1,"mlm\\",phname[i],"_mlm_all.fastGWA",sep=""),sep="\t",header=T)
  data=data[!is.na(data$P),]
  data=data[,c(1,2,3,10)]
  names(data)=c("CHR","SNP","BP","P")
  lam <- estlambda(data$P,plot = F) 
  chis <- qchisq(data$P,df=1,lower.tail = FALSE)
  data$P <- pchisq(chis/lam$estimate,df=1, lower.tail = FALSE)
  png(paste(path2,"raw_data\\qq_all_",phname[i],".png",sep=""))
  lam <- estlambda(data$P,plot = T)
  text(5, 15, paste("lambda =",substr(as.character(lam$estimate),1,6)))
  dev.off()
  png(paste(path2,"raw_data\\all_",phname[i],".png",sep=""))
  qqman::manhattan(data,col = c("blue4", "orange3"))
  dev.off()
}
clnum<-detectCores() 
cl <- makeCluster(getOption("cl.cores", clnum));
clusterExport(cl,deparse(substitute(draw_manha_raw)))
clusterEvalQ(cl,library(GenABEL))
clusterExport(cl,c("phname","effsnp","path1","gcta","path2"))
parLapply(cl, 3:29,  draw_manha_raw)
stopCluster(cl)
##cojo data
draw_manha<-function(i){
  data=read.csv(file=paste(path1,"cojo\\","cojo_result_",phname[i],".cma.cojo",sep=""),sep="\t",header=T)
  data=data[,c(2,1,3,13)]
  data=data[which(!is.na(data$pC)),]
  names(data)<-c("SNP", "CHR", "BP", "P")
  png(paste(path2,"cojo\\",phname[i],"_maha_cojo.png",sep=""))
  qqman::manhattan(data,col = c("blue4", "orange3"),main=phname[i],annotateTop=T,suggestiveline = -log10(effsnp[(i-2),5]),genomewideline =-log10(effsnp[(i-2),6] ))
  dev.off()
  write(phname[i],file="C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.09.27\\Text.txt",append=TRUE)
  
  
  data=read.csv(file=paste(path1,"cojo\\","cojo_result_all_",phname[i],".cma.cojo",sep=""),sep="\t",header=T)
  data=data[,c(2,1,3,13)]
  data=data[which(!is.na(data$pC)),]
  names(data)<-c("SNP", "CHR", "BP", "P")
  png(paste(path2,"cojo\\",phname[i],"_maha_all_cojo.png",sep=""))
  qqman::manhattan(data,col = c("blue4", "orange3"),main=paste(phname[i],"all"),annotateTop=T,suggestiveline = -log10(effsnp[(i-2),11]),genomewideline =-log10(effsnp[(i-2),12] ))
  dev.off() 
  write(paste("all",phname[i]),file="C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.09.05\\Text.txt",append=TRUE)
}
clnum<-detectCores() 
cl <- makeCluster(getOption("cl.cores", clnum));
clusterExport(cl,deparse(substitute(draw_manha)))
clusterEvalQ(cl,library(GenABEL))
clusterEvalQ(cl,library(qqman))
clusterExport(cl,c("phname","effsnp","path1","path2"))
parLapply(cl, 3:29,  draw_manha)
stopCluster(cl)
##check result number
#get file name in path
result_name=list.files(paste(path2,"cojo\\",sep=""))
result_name_match=cbind.data.frame(result_name_match,data.frame(matrix(nrow=27,ncol=4)))
result_name_match[,5]=paste(phname[3:29],rep("_maha_cojo.png"),sep="")
result_name_match[,7]=paste(phname[3:29],rep("_maha_all_cojo.png"),sep="")
for(i in 1:27){
  name=result_name_match[i,5]
  if(name%in%result_name){
    result_name_match[i,6]=1
  }
  name=result_name_match[i,7]
  if(name%in%result_name){
    result_name_match[i,8]=1
  }
}
###locuszoom plot
#trans plink to VCF
names(phdata(fl_1))
for(i in c(3,15,16,17,18,19,22,23,26)){
  trait=names(phdata(fl_1))[i]
  cmd=paste(plink,"--bfile",paste(path1,"out_",trait,sep=""),"--export vcf","--out",paste(path1,"vcf\\","vcf_",trait,sep=""))
  system(cmd)  
}
for(i in c(8,11,12,14,16,17,18,20,22,25,27,29)){
  trait=names(phdata(fl_1))[i]
  cmd=paste(plink,"--bfile",paste(path1,"out_all_",trait,sep=""),"--export vcf","--out",paste(path1,"vcf\\","vcf_all_",trait,sep=""))
  system(cmd)  
}


rm(list=setdiff(ls(),c("fl_1","fl.data.merge","effsnp","her_trait","gtf","hmp","gwa")))
gc()
#install and library
#install.packages(c("devtools","remotes"))
#install.packages(c("ggplot2","ggrepel","reshape2"))
#BiocManager::install(c("SNPRelate","gdsfmt"))
#library(remotes)
#install_github("whweve/IntAssoPlot",build=TRUE,build_vignettes = TRUE)

library(IntAssoPlot)
#load gtf file
gtf=read.csv("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\工作文件\\Arabidopsis_thaliana.TAIR10.54.gtf",header = F,sep="\t")
head(gtf)
#load hapmap
hmp=read.csv(paste(path1,"hmp\\all_individual_and_SNP.hmp.txt",sep=""),header = T,sep="\t",stringsAsFactors = F)
head(hmp)
names(hmp)[1]="rs"
splitfun2=function(a){
  as.integer(strsplit(as.character(a[1]),":")[[1]][2])
}
memory.limit(1000000)
hmp$pos=apply(hmp,1,splitfun2)  

#load gwas file
names(phdata(fl_1))
trait=names(phdata(fl_1))[3]
trait=paste("all_",trait,sep="")
gwa=read.csv(paste(path1,"cojo\\cojo_result_",trait,".cma.cojo",sep=""),header=T,sep="\t")
head(gwa)
gwa2=gwa
gwa=gwa[,c(2,1,3,13)]
gwa[which(gwa$p==min(gwa$p,na.rm=T)),]
gwa=gwa[which(!is.na(gwa$pC)),]
names(gwa)=c("Marker","Locus","Site","p")
splitfun=function(a){
  as.integer(strsplit(as.character(a[1]),":")[[1]][2])
}
gwa$Site =apply(gwa,1,splitfun)  

#find significant snp
gwa1=gwa[which(gwa$Locus==1),]
#gwa1=gwa1[which((gwa1$Site>11000000)&(gwa1$Site<12000000)),]
gwa1[which(gwa1$p==min(gwa1$p,na.rm=T)),]
gwa1[order(gwa1$p),][1:10,]

marker=gwa1[which(gwa1$p<1e-6),c(1,2,3)]
names(marker)=c("rs","chrom","pos")

#draw
chr=1
site=10453719
range=20000
path=paste("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.10.06\\","",trait,"_",chr,"_",site,"_",range,".png",sep="")
threshold=-log10(1e-5)
png(path,res=300,width=3000, height=2000)
pal <- colorRampPalette(c("white", "red"))
IntRegionalPlot(chr=chr,left=(site-range),right=(site+range),gtf=gtf,association=gwa,hapmap=hmp,hapmap_ld=hmp,threshold=threshold,leadsnp_size=2,label_gene_name = TRUE,marker2label = marker,marker2label_angle = 0,marker2label_size = 3,colour02 = pal(5)[1],colour04 = pal(5)[2],colour06 = pal(5)[3],colour08 = pal(5)[4],colour10 = pal(5)[5])+
  ggtitle(trait)

dev.off()
gc()

#view box plot
name=paste(as.character(chr),as.character(site),sep=":")
table=as.data.frame(as.character(fl_1[,name]))
table=cbind.data.frame(table,phdata(fl_1)[,3])
names(table)=c("mut","value")
boxplot(value~mut,data=table,varwidth = T,xlab = paste(c(as.character(as.numeric(table(table$mut)))),names(table(table$mut)),collapse="   "),main=name)  

###draw SNP beta bar plot
SNP_list=c("1:10453719","2:5878663","3:10524813","5:816974","3:10524159","2:11543244","1:11016927","3:10824568","3:7916048","4:16360801","1:11009297","5:811307","1:18963726","5:22229967")
duplicated(SNP_list)
beta_data=data.frame(matrix(nrow=length(SNP_list),ncol=6*27))
row.names(beta_data)=SNP_list
colnames(beta_data)=c(paste(rep(phname[3:29],each=3),rep(c("","_se","_p"),27),sep=""),paste(rep("all_",3*27),rep(phname[3:29],each=3),rep(c("","_se","_p"),27),sep=""))
for(i in seq(1,(6*27),3)){
  trait=names(beta_data)[i]
  print(i)
  print(trait)
  gwa=read.csv(paste(path1,"cojo\\cojo_result_",trait,".cma.cojo",sep=""),header=T,sep="\t")
  gwa=gwa[,c(1,2,11,12,13)]
  gwa$SNP=as.character(gwa$SNP)
  for(j in 1:nrow(beta_data)){
    data=gwa[which(gwa$SNP==row.names(beta_data)[j]),]
    beta_data[j,i:(i+2)]=data[1,3:5]
  }
}

##summary each snp's data. i means the SNP
for(j in 1:14){
  data=beta_data[j,]
  snp_data=data.frame(matrix(nrow=9,ncol=18))
  row.names(snp_data)=substring(phname[3:11],4)
  names(snp_data)=c(paste(rep(rep(c("CD","HY","HY-CD"),each=3),2),rep(c("main","all"),each=9),rep(c("beta","se","p"),6),sep="_"))
  for(i in seq(1,(6*27),3)){
    if(i<28){
      snp_data[(i+2)/3,1:3]=data[1,i:(i+2)]
    }else if((28<=i)&(i<55)){
      snp_data[(i-25)/3,4:6]=data[1,i:(i+2)]
    }else if((55<=i)&(i<82)){
      snp_data[(i-52)/3,7:9]=data[1,i:(i+2)]
    }else if((82<=i)&(i<109)){
      snp_data[(i-79)/3,10:12]=data[1,i:(i+2)]
    }else if((109<=i)&(i<136)){
      snp_data[(i-106)/3,13:15]=data[1,i:(i+2)]
    }else if(136<=i){
      snp_data[(i-133)/3,16:18]=data[1,i:(i+2)]
    }
  }
  snp_data$name=row.names(snp_data)
  
  ##draw
  data=snp_data
  data2=data.frame(matrix(nrow=54,ncol=5))
  mark=1
  for(i in seq(1,18,3)){
    data2[mark:(mark+8),1:3]=data[1:9,i:(i+2)]
    data2[mark:(mark+8),4]=row.names(snp_data)
    mark=mark+9
  }
  data2[,5]=rep(c("CD_main","HY_main","HY-CD_main","CD_all","HY_all","HY-CD_all"),each=9)
  data2$X3=round(-log10(data2$X3),2)
  xx=factor(data2$X4,levels=c(unique(data2$X4)))
  ggplot(data2,aes(x=xx,y=X1,fill=xx),ylab="",xlab="")+
    geom_bar(stat="identity",width=0.5,color="black")+
    geom_errorbar(aes(ymin=X1-X2, ymax=X1+X2), width=0)+
    geom_text(aes(label = paste("p",X3)),size=3,hjust = -0.1,vjust=1.2)+
    xlab("trait")+ylab("Beta")+ggtitle(paste(row.names(beta_data)[j],"beta of all trait and all group"))+
    guides(fill=FALSE)+
    facet_grid(~ X5)+
    theme_bw()+coord_flip()
  ggsave(file=paste(path2,substr(row.names(beta_data)[j],1,1),"_",substring(row.names(beta_data)[j],3),"beta.png",sep=""),dpi=300,width=18,height = 6)
}


###snp allele box plot


par(mfrow=c(2,2))
name=as.character("5:811307")
trait=14
table=as.data.frame(as.character(fl_1[,name]))
table=cbind.data.frame(table,phdata(fl_1)[,trait])
names(table)=c("mut","value")
boxplot(value~mut,data=table,varwidth = T,xlab = paste(c(as.character(as.numeric(table(table$mut)))),names(table(table$mut)),collapse="   "),main=paste(name,"main"))

table=as.data.frame(as.character(fl.data.merge[,name]))
table=cbind.data.frame(table,phdata(fl.data.merge)[,trait])
names(table)=c("mut","value")
boxplot(value~mut,data=table,varwidth = T,xlab = paste(c(as.character(as.numeric(table(table$mut)))),names(table(table$mut)),collapse="   "),main=paste(name,"all"))


###location analysis
###load map info
library(ggplot2)
library(plyr)
library(maptools)
library(sp)


map_info=read.csv("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\工作文件\\1001_genome_accessions.txt",sep=",",header=T)
load("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.10.02\\worldmap.RData")
source("C:\\Users\\MSI\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.10.02\\map_plot.R")

names(map_info)[6:7]=c("lat2","long2")
name=as.character("5:811307")
table=as.data.frame(as.character(fl_1[,name]))
table=as.data.frame(as.character(fl.data.merge[,name]))
table=as.data.frame(as.character(fl.data[,name]))

#table=cbind.data.frame(table,phdata(fl_1)[,trait])
table$X88=row.names(table)
table$mark=1
table$mark[table[,1]==unique(table[,1])[1]]=2
table=merge.data.frame(table,map_info,by="X88")

png(paste(path2,substr(name,1,1),"_",substring(name,3),"_map.png",sep=""),res=300,width=2300, height=2000)
plot_maps(wrld_map= world_map, layers="", region="World",main=name)
points(map_info$long2,map_info$lat2,col="black",pch=19,cex=0.2)
points(table$long2[which(table$mark==2)],table$lat2[which(table$mark==2)],col="red",pch=19,cex=0.2)
points(table$long2[which(table$mark==1)],table$lat2[which(table$mark==1)],col="blue",pch=19,cex=0.2)
legend("topright",legend=c("1001","major","minor"),col=c("black","red","blue"), lty=1,lwd=2)
dev.off()

#caculate LD
library(genetics)
snp_names=c(as.character(gwa[which((gwa$Locus==chr)&(gwa$Site<(site+range))&(gwa$Site>(site-range))),1]))
LD_matrix=LD(as.genotype(fl_1[,snp_names]))$R^2
LD_matrix[is.na(LD_matrix)]=0
diag(LD_matrix)=1
#install.packages("pheatmap")
library(pheatmap)
pheatmap(LD_matrix, 
         # annotation_row=dfGene, # （可选）指定行分组文件
         # annotation_col=dfSample, # （可选）指定列分组文件
         show_colnames = F, # 是否显示列名
         show_rownames=F,  # 是否显示行名
         fontsize=2, # 字体大小
         color = colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(50), # 指定热图的颜色
         annotation_legend=TRUE, # 是否显示图例
         border_color=NA,  # 边框颜色 NA表示没有
         scale="row",  # 指定归一化的方式。"row"按行归一化，"column"按列归一化，"none"不处理
         cluster_rows = F, # 是否对行聚类
         cluster_cols = F # 是否对列聚类
)



###性状相关性分析
library(Hmisc)
a=rowSums(Add_phe[,2:ncol(Add_phe)])
dd=Add_phe[!is.na(a),2:ncol(Add_phe)]
head(dd)
res2 <- rcorr(as.matrix(dd))
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
chart.Correlation(dd, histogram=TRUE, pch=19)



###boxplot of significant peaks
data=read.csv(file=paste(path1,"cojo\\","cojo_result_all_",phname[i],".cma.cojo",sep=""),sep="\t",header=T)
data=data[which(data$Chr==1),]
data=data[,c(2:ncol(data))]
data$bp=apply(data,1,splitfun)  
data=data[which((data$bp>11000000)&(data$bp<11050000)),]
data=data[order(data$pC)[1:10],]


###down load expression data
BiocManager::install("GEOquery")
library(GEOquery)
BiocManager::install("edgeR")
install.packages("edgeR")
library(edgeR)
