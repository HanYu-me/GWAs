library(GenABEL)
library(readxl)
###load data
load("C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\工作数据\\2022.07.11\\fl.data.RData")
phe=as.data.frame(read_excel(file.choose(),sheet=1,col_names=T))
for(i in 12:ncol(phe)){
  phe[,i]=as.numeric(phe[,i])
}
table(phe$Commongarden)
CD_phe=phe[which(phe$Commongarden=="Chengdu"),]
HY_phe=phe[which(phe$Commongarden=="Hongyuan"),]
#load id
arab_name=as.data.frame(read_excel(file.choose(),sheet=1,col_names=T))
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
#remove the ph data of individuals which CV>0.05
#caculate the mean ph data
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


#去除重复
Add_phe_dup <- Add_phe[which(Add_phe$id %in% Add_phe[duplicated(Add_phe$id),1]),]
Add_phe=Add_phe[which(!Add_phe$id%in%Add_phe_dup$id),]
sort(table(Add_phe_dup$id))
sort(table(Add_phe$id))
#导入经纬度文件由于拟南芥生态型种源地地理信息文件和 生态型拟南芥种源地信息文件的经纬度一致，所以只对比1001文件
data2=read.csv(file.choose(),header=T,sep=";")
data2=data2[,c(1,3,6,7)]
data=arab_name[which(arab_name$id_1001 %in% Add_phe_dup$id),c(9,2,4,5)]
names(data)[1]="id"
data3=merge.data.frame(data,data2,by="id")
data3=data3[,c(1,2,3,6,4,7,5)]
k=data3[which((data3$latitude.x!=data3$latitude.y)|(data3$longitude.x!=data3$longitude.y)),]
unique(k$id)
#肉眼对比经纬度，0.1数量级差异的数据都剔除c(5023,5768,6898,6982,7123,7133)
Add_phe_dup=Add_phe_dup[which(!Add_phe_dup$id %in% unique(k$id)),]
name=unique(Add_phe_dup$id)
for( i in 1:length(name)){
  a=Add_phe_dup[which(Add_phe_dup$id==name[i]),]
  #人工筛查差异大的数据
  print(a)
  for(j in 1:ncol(a)){
    a[1,j]=median(a[,j],na.rm=T)
  }
  Add_phe=rbind.data.frame(Add_phe,a[1,])
}
#添加成都和红原相减的表型
for(i in 1:8){
  Add_phe=cbind.data.frame(Add_phe,(Add_phe[,(i+9)]-Add_phe[,(i+1)]))
  names(Add_phe)[17+i]=paste(names(Add_phe)[i+9],"_",names(Add_phe)[i+1],sep="")
}
#添加表型
fl.data.merge=add.phdata(fl.data,Add_phe)
phdata(fl.data.merge)[1:10,]


#删除表型
names(phdata(fl.data.merge)[1,])
fl.data.merge=del.phdata(fl.data.merge,c("FT10_mean","FT16_mean","FT10_sd","FT16_sd"))
names(phdata(fl.data.merge)[1,])
#刘伟QC去除了"108","9963","9998","10022"的个体
fl.data.merge=fl.data.merge[which(fl.data.merge@phdata$id%in%Add_phe$id),]
nids(fl.data.merge)
###Quality Control
#QC to all individual
qc1 <- check.marker(fl.data.merge, p.level=0)
fl.data.merge=fl.data.merge[qc1$idok,qc1$snpok]
#population structure and QC to main group
par(mfrow=c(1,1))
ibs_m=ibs(fl.data.merge,weight="freq")
dis_m=as.dist(0.5-ibs_m)
cmd=cmdscale(dis_m)
km <- kmeans(cmd,centers=4, nstart=1000,iter.max = 100)
plot(cmd[,1],cmd[,2])
points(cmd[which(km$cluster==1),1],cmd[which(km$cluster==1),2],col="red")
points(cmd[which(km$cluster==2),1],cmd[which(km$cluster==2),2],col="blue")
points(cmd[which(km$cluster==3),1],cmd[which(km$cluster==3),2],col="orange")
#choose the smallest group as outer group
sub_id2=names(km$cluster[which(km$cluster==1)])
sub_id1=names(km$cluster[which(km$cluster!=1)])
#将数据分成两个群体
fl_1=fl.data.merge[which(fl.data.merge@phdata$id%in% sub_id1),]
qc2 <- check.marker(fl_1, p.level=0)

fl_1=fl_1[qc2$idok,qc2$snpok]
fl_1=fl_1[which(!fl_1@phdata$id%in%c("108","9963","9998","10022")),]
save(fl_1,file="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\工作数据\\2022.08.23\\fl_1.RData")

fl.data.merge=fl.data.merge[which(!fl.data.merge@phdata$id%in%c("108","9963","9998","10022")),]
save(fl.data.merge,file="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\工作数据\\2022.08.23\\fl.data.merge.RData")

nids(fl_1)
nids(fl.data.merge)

###去除极端CV的个体进行后续分析
CD_CV$id=1:240
CD_CV$id=as.numeric(as.character(arab_name[which(arab_name$label==CD_CV$id),9]))
CD_CV=CD_CV[,c(9,1:8)]
All_CV=cbind.data.frame(CD_CV,HY_CV)

###计算遗传力和se
#导出数据
#plink 代码  plink --tfile plink --make-bed  --out out

#计算亲缘关系矩阵
gcta <- "C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\gcta_v1.94.0Beta_windows_x86_64"
GCTA_estKin <- function(gcta,output,input){
  cmd_k <- paste(gcta,"--bfile",input, "--make-grm" ,"--make-grm-alg 1","--out",output,sep=" ")
  system(cmd_k)
}
#计算遗传力
GCTA_REML <- function(t1,output,gcta,phenotype,grm){
  cmd_biGREML <- paste(gcta,"--pheno",phenotype,"--grm",grm,"--reml","--mpheno",t1, "--out",output,sep=" ")
  cat(cmd_biGREML,"\n")
  system(cmd_biGREML)
  #print(paste(output,"is finished"))
}
###进行mlm
#mlm函数
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
#导出表型文件 并且计算遗传力

for(i in 3:26){
  print(i)
  if(i<19){
    fl_2=fl_1[which(!fl_1@phdata$id%in%All_CV[which(All_CV[,(i-1)]>quantile(All_CV[,(i-1)],0.95,na.rm = T)),1]),]  
    export.plink(fl_2)
    plink <- "C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\plink_win64_20220402\\plink"
    cmd_plink=paste(plink,"--tfile C:\\Users\\98033\\OneDrive\\Documents\\plink --make-bed --out",paste("C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\out_",names(phdata(fl_1))[i],sep=""),sep=" ")
    system(cmd_plink)
  }else{
    export.plink(fl_1)
    fl_2=fl_1
    plink <- "C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\plink_win64_20220402\\plink"
    cmd_plink=paste(plink,"--tfile C:\\Users\\98033\\OneDrive\\Documents\\plink --make-bed --out",paste("C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\out_",names(phdata(fl_1))[i],sep=""),sep=" ")
    system(cmd_plink)
  }
  path1="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\"
  GCTA_estKin(gcta=gcta,output = paste(path1,"oak",sep=""),input=paste(path1,"out",sep=""))
  phe_file=phdata(fl_2)[,c(1,2,i)]
  phe_file[,2]=as.numeric(phe_file[,1])
  phe_file[,1]=1:nrow(phe_file)
  phe_file[,1]=as.numeric(phe_file[,1])
  file_path=paste(path1,names(phdata(fl_2))[i],"_phe.txt",sep="")
  write.table(phe_file,file = file_path,col.names = F,row.names = F)
  GCTA_REML(t1=1,output=paste(path1,names(phdata(fl_1))[i],sep=""),gcta=gcta,phenotype=paste(path1,names(phdata(fl_1))[i],"_phe.txt",sep=""),grm="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\oak")
  GCTA_mlm(gcta=gcta,bfile=paste(path1,paste("C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\out_",names(phdata(fl_1))[i],sep=""),sep=""),grm=paste(path1,"sp_oak",sep=""),pheno=paste(path1,names(phdata(fl_1))[i],"_phe.txt",sep=""),mlm=paste(path1,names(phdata(fl_1))[i],"_mlm",sep=""),mark=1)  
  
}

#计算全部个体所在群体的遗传力

for(i in 3:26){
  print(i)
  if(i<19){
    fl_2=fl.data.merge[which(!fl.data.merge@phdata$id%in%All_CV[which(All_CV[,(i-1)]>quantile(All_CV[,(i-1)],0.95,na.rm = T)),1]),]  
    export.plink(fl_2)
    plink <- "C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\plink_win64_20220402\\plink"
    cmd_plink=paste(plink,"--tfile C:\\Users\\98033\\OneDrive\\Documents\\plink --make-bed --out",paste("C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\out_all_",names(phdata(fl_1))[i],sep=""),sep=" ")
    system(cmd_plink)
  }else{
    export.plink(fl.data.merge)
    fl_2=fl.data.merge
    plink <- "C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\plink_win64_20220402\\plink"
    cmd_plink=paste(plink,"--tfile C:\\Users\\98033\\OneDrive\\Documents\\plink --make-bed --out",paste("C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\out_all_",names(phdata(fl_1))[i],sep=""),sep=" ")
    system(cmd_plink)
  }
  path1="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\"
  GCTA_estKin(gcta=gcta,output =paste(path1,"oak2",sep=""),input=paste(path1,"out2",sep=""))
  phe_file=phdata(fl_2)[,c(1,2,i)]
  phe_file[,2]=as.numeric(phe_file[,1])
  phe_file[,1]=1:nrow(phe_file)
  phe_file[,1]=as.numeric(phe_file[,1])
  file_path=paste(path1,names(phdata(fl_1))[i],"_phe_all.txt",sep="")
  write.table(phe_file,file = file_path,col.names = F,row.names = F)
  GCTA_REML(t1=1,output=paste(path1,names(phdata(fl_1))[i],"all",sep=""),gcta=gcta,phenotype=paste(path1,names(phdata(fl_1))[i],"_phe_all.txt",sep=""),grm="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\oak2")
  GCTA_mlm(gcta=gcta,bfile=paste(path1,paste("C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\out_all_",names(phdata(fl_1))[i],sep=""),sep=""),grm=paste(path1,"sp_oak2",sep=""),pheno=paste(path1,names(phdata(fl_1))[i],"_phe_all.txt",sep=""),mlm=paste(path1,names(phdata(fl_1))[i],"_mlm_all",sep=""),mark=1)  
}

#统计遗传力信息
herita=data.frame(matrix(nrow=24,ncol=6))
names(herita)=c("main group h2","se","p","all h2","se","p")
row.names(herita)=names(phdata(fl_1))[3:26]
for(i in 3:26){
  path=paste(path1,names(phdata(fl_1))[i],".hsq",sep="")
  her_trait=read.csv(file=path,stringsAsFactors = F,sep="\t")
  herita[(i-2),1:3]=as.numeric(c(her_trait[4,2:3],her_trait[9,2]))
  path=paste(path1,names(phdata(fl_1))[i],"all.hsq",sep="")
  her_trait=read.csv(file=path,stringsAsFactors = F,sep="\t")
  herita[(i-2),4:6]=as.numeric(c(her_trait[4,2:3],her_trait[9,2]))
}


###cojo
for(i in 3:26){
  data=read.csv(file=paste(path1,names(phdata(fl_1))[i],"_mlm.fastGWA",sep=""),sep="\t",header=T)
  data=data[,c(2,4,5,7,8,9,10,6)]
  lam <- estlambda(data$P,plot = T) 
  chis <- qchisq(data$P,df=1,lower.tail = FALSE)
  ps <- pchisq(chis/lam$estimate,df=1, lower.tail = FALSE)
  data$P=ps
  names(data)=c("SNP","A1","A2","freq","b","se","p","N")
  write.table(data,file=paste(path1,"cojo_",names(phdata(fl_1))[i],".ma",sep=""),sep="\t",row.names = F,quote = FALSE)
  cmd <- paste(gcta, "--bfile",paste(path1,"out_",names(phdata(fl_1))[i],sep=""), "--cojo-file",paste(path1,"cojo_",names(phdata(fl_1))[i],".ma",sep=""), "--cojo-p 5e-08", "--cojo-collinear 0.1", "--cojo-wind 200", "--cojo-slct","--maf 0.05", "--out", paste(path1,"cojo_result_",names(phdata(fl_1))[i],sep=""))
  system(cmd)
}
for(i in 3:26){
  data=read.csv(file=paste(path1,names(phdata(fl_1))[i],"_mlm_all.fastGWA",sep=""),sep="\t",header=T)
  data=data[,c(2,4,5,7,8,9,10,6)]
  lam <- estlambda(data$P,plot = T) 
  chis <- qchisq(data$P,df=1,lower.tail = FALSE)
  ps <- pchisq(chis/lam$estimate,df=1, lower.tail = FALSE)
  data$P=ps
  names(data)=c("SNP","A1","A2","freq","b","se","p","N")
  write.table(data,file=paste(path1,"cojo_",names(phdata(fl_1))[i],"_all_.ma",sep=""),sep="\t",row.names = F,quote = FALSE)
  cmd <- paste(gcta, "--bfile",paste(path1,"out_all_",names(phdata(fl_1))[i],sep=""), "--cojo-file",paste(path1,"cojo_",names(phdata(fl_1))[i],"_all_.ma",sep=""), "--cojo-p 5e-08", "--cojo-collinear 0.1", "--cojo-wind 200", "--cojo-slct","--maf 0.05", "--out", paste(path1,"cojo_result_all_",names(phdata(fl_1))[i],sep=""))
  system(cmd)
}

###
for(i in c(11:25)){
  data=read.csv(file=paste(path1,"cojo_result_",names(phdata(fl_1))[i],".cma.cojo",sep=""),sep="\t",header=T)
  data=data[,c(2,1,3,13)]
  data=data[which(!is.na(data$pC)),]
  names(data)<-c("SNP", "CHR", "BP", "P")
  path="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.09.01\\"
  png(paste(path,names(phdata(fl_1))[i],"_maha.png",sep=""))
  qqman::manhattan(data,col = c("blue4", "orange3"),main=names(phdata(fl_1))[i],annotateTop=T)
  dev.off()
}
for(i in c(20:21,23:24)){
  data=read.csv(file=paste(path1,"cojo_result_all_",names(phdata(fl_1))[i],".cma.cojo",sep=""),sep="\t",header=T)
  data=data[,c(2,1,3,13)]
  data=data[which(!is.na(data$pC)),]
  names(data)<-c("SNP", "CHR", "BP", "P")
  path="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.09.01\\"
  png(paste(path,names(phdata(fl_1))[i],"_maha_all.png",sep=""))
  qqman::manhattan(data,col = c("blue4", "orange3"),main=paste(names(phdata(fl_1))[i],"all"),annotateTop=T)
  dev.off()
}

data=read.csv(file.choose(),sep="\t",header=T)
data=data[,c(2,4,5,7,8,9,10,6)]
write.table(data,file=paste(path1,"cojo_",names(phdata(fl_1))[i],sep=""),sep="\t",row.names = F,quote = FALSE)

cmd <- paste(gcta, "--bfile",paste(path1,"out_",names(phdata(fl_1))[i],sep=""), "--cojo-file",paste(path1,"cojo_",names(phdata(fl_1))[i],sep=""), "--cojo-p 2.512472e-08", "--cojo-collinear 0.1", "--cojo-wind 200", "--cojo-slct", "--out", "cojo")
system(cmd)

data=read.csv(file.choose(),sep="\t")
data2=data[,c(2,1,3,13)]
names(data2)<-c("SNP", "CHR", "BP", "P") 
png("C:\\Users\\98033\\OneDrive\\桌面\\maha.png")
qqman::manhattan(data2,col = c("blue4", "orange3"))
dev.off()



###可视化
#所有性状的两种群体
lambda_all=data.frame(matrix(ncol=2,nrow=24))
names(lambda_all)=c("main_raw","all_raw")
row.names(lambda_all)=names(phdata(fl_1))[3:26]
for(i in 3:26){
  path=paste(path1,names(phdata(fl_1))[i],"_mlm.fastGWA",sep="")
  res3=read.csv(file=path,sep="\t",header=T)
  res3=data.frame(res3$SNP, res3$CHR, res3$POS, res3$P) 
  names(res3)<-c("SNP", "CHR", "BP", "P") 
  path="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.08.23\\"
  #qq
  png(paste(path,names(phdata(fl_1))[i],"_QQplot.png",sep=""))
  lam <- estlambda(res3$P,plot = T) 
  text(5, 15, paste("lambda =",substr(as.character(lam$estimate),1,6)))
  dev.off()
  chis <- qchisq(res3$P,df=1,lower.tail = FALSE)
  ps <- pchisq(chis/lam$estimate,df=1, lower.tail = FALSE)
  kk=ps
  
  lambda_all[(i-2),1]=lam$estimate
  
  #manhattan
  png(paste(path,names(phdata(fl_1))[i],"_maha.png",sep=""))
  qqman::manhattan(res3,col = c("blue4", "orange3"))
  dev.off()
  
  res3$P=kk
  png(paste(path,names(phdata(fl_1))[i],"_QQplot2.png",sep=""))
  lam <- estlambda(res3$P,plot = T) 
  text(5, 15, paste("lambda =",substr(as.character(lam$estimate),1,6)))
  dev.off()
  png(paste(path,names(phdata(fl_1))[i],"_maha2.png",sep=""))
  qqman::manhattan(res3,col = c("blue4", "orange3"))
  dev.off()
  print(paste(names(phdata(fl_1))[i],"finished"))
  
  ##all
  path=paste(path1,names(phdata(fl_1))[i],"_mlm_all.fastGWA",sep="")
  res3=read.csv(file=path,sep="\t",header=T)
  res3=data.frame(res3$SNP, res3$CHR, res3$POS, res3$P) 
  names(res3)<-c("SNP", "CHR", "BP", "P") 
  path="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.08.23\\"
  #qq
  png(paste(path,names(phdata(fl_1))[i],"_QQplot_all.png",sep=""))
  lam <- estlambda(res3$P,plot = T) 
  text(5, 15, paste("lambda =",substr(as.character(lam$estimate),1,6)))
  dev.off()
  chis <- qchisq(res3$P,df=1,lower.tail = FALSE)
  ps <- pchisq(chis/lam$estimate,df=1, lower.tail = FALSE)
  kk=ps
  
  lambda_all[(i-2),2]=lam$estimate
  
  #manhattan
  png(paste(path,names(phdata(fl_1))[i],"_maha_all.png",sep=""))
  qqman::manhattan(res3,col = c("blue4", "orange3"))
  dev.off()
  
  
  res3$P=kk
  png(paste(path,names(phdata(fl_1))[i],"_QQplot_all2.png",sep=""))
  lam <- estlambda(res3$P,plot = T) 
  text(5, 15, paste("lambda =",substr(as.character(lam$estimate),1,6)))
  dev.off()
  png(paste(path,names(phdata(fl_1))[i],"_maha_all2.png",sep=""))
  qqman::manhattan(res3,col = c("blue4", "orange3"))
  dev.off()
  print(paste(names(phdata(fl_1))[i],"all finished"))
}

###保存遗传力和lambda

write.table(herita,file = "C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.08.11\\heritability.txt",col.names = T,row.names = T)
write.table(lambda_all,file = "C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.08.11\\lambda.txt",col.names = T,row.names = T)

###绘制CV变化图
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

for(i in 1:8){
  path="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.08.11\\"
  png(paste(path,"CD_",names(CD_CV)[i],"_CV.png",sep=""))
  par(mfrow=c(2,2))
  hist(CD_CV2_main[,i],breaks=90,main=paste("Main group raw CV of CD",names(CD_CV)[i]),xlab = "value")
  hist(CD_CV2_all[,i],breaks=90,main=paste("All group raw CV of CD",names(CD_CV)[i]),xlab = "value")
  hist(CD_CV_main[,i],breaks=90,main=paste("Main group adjusted CV of CD",names(CD_CV2)[i]),xlab = "value")  
  hist(CD_CV_all[,i],breaks=90,main=paste("All group adjusted CV of CD",names(CD_CV2)[i]),xlab = "value")  
  dev.off()
}
for(i in 1:8){
  path="C:\\Users\\98033\\OneDrive\\桌面\\拟南芥or柳树\\结果保存\\2022.08.11\\"
  png(paste(path,"HY_",names(HY_CV)[i],"_CV.png",sep=""))
  par(mfrow=c(2,2))
  hist(HY_CV2_main[,i],breaks=90,main=paste("Main group raw CV of HY",names(HY_CV)[i]),xlab = "value")
  hist(HY_CV2_all[,i],breaks=90,main=paste("All group raw CV of HY",names(HY_CV)[i]),xlab = "value")
  hist(HY_CV_main[,i],breaks=90,main=paste("Main group adjusted CV of HY",names(HY_CV2)[i]),xlab = "value")  
  hist(HY_CV_all[,i],breaks=90,main=paste("All group adjusted CV of HY",names(HY_CV2)[i]),xlab = "value")  
  dev.off()
}

###性状相关性分析
library(Hmisc)
a=rowSums(Add_phe[,2:ncol(Add_phe)])
dd=Add_phe[!is.na(a),2:ncol(Add_phe)]
head(dd)
res2 <- rcorr(as.matrix(dd))
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
chart.Correlation(dd, histogram=TRUE, pch=19)






