library(hglm)
library(ggplot2)
library(GenABEL)
library(data.table)
library(plyr)
library(tidyr)
###function
c.z.hglm <- function(kin){
  relmat <- kin
  relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]#上三角变成下三角
  svd <- svd(relmat)
  Z <- svd$u %*% diag(sqrt(svd$d)) #左边的矩阵乘奇异值
  return(Z)
}
var.ex <- function(snp,gwaadata,ph,kin,record=1,geno){
  if(!require(hglm))
    require(hglm)
  # creat Z
  na <- !is.na(ph)
  Z <- c.z.hglm(kin = kin[na,na])
  # fit only IBS
  
  onlyIBS.hglm <- hglm(X = matrix(rep(1, length(gwaadata@gtdata@idnames[na]))), y = ph[na], Z = Z)
  #cat("fitting Null model","\n")
  # fit snp model
  
  geno=geno[,snp]
  na <- complete.cases(cbind(ph,geno))
  withsnp.hglm <- hglm(X = cbind(1, geno[na,]), y = ph[na], Z = Z)
  #cat("fitting SNP model","\n")
  r1 <- onlyIBS.hglm$varRanef/(onlyIBS.hglm$varRanef+ onlyIBS.hglm$varFix)
  r2 <- withsnp.hglm$varRanef/(withsnp.hglm$varRanef+ withsnp.hglm$varFix)
  #variation by snp
  R2.all <- 1- var(ph[na] - cbind(1, geno[na,])%*%withsnp.hglm$fixef)/var(ph[na])
  kin_R2=1-var(ph[na]-cbind(1, geno[na,])%*%withsnp.hglm$fixef-Z%*%withsnp.hglm$ranef)/var(ph[na])
  R2 <- numeric(length(snp))
  for( i in 1:length(snp)){
    R2[i] <- 1- var(ph[na] - cbind(1, geno[na,snp[i]])%*%withsnp.hglm$fixef[c(1,i+1)])/var(ph[na])
  }
  names(R2) <- snp
  # phe by kinship
  r2.tot <- r2*(1-R2.all)
  r.h2 <- (r1-r2.tot)/r1
  r.tot <- (r1-r2.tot)/R2.all
  eff <- withsnp.hglm$fixef
  s.g <- summary(withsnp.hglm)
  p <- s.g$FixCoefMat[,4]
  sd <- s.g$FixCoefMat[,2]
  eff[2:length(snp)] <- eff[2:length(snp)]
  names(eff)[1] <- "u" 
  if(record==1){
    return(list("Var.kin"=r.h2,"h2"=r1,"r2"=r2,"Var.add"=r.tot,"p"=p,"eff"=eff,"R2"=R2,"R2.all"=R2.all,"model.snp"=withsnp.hglm))
  }else if(record==2){
    #return all pve without correlation
    return(c(r2.tot,R2.all))
  }
  else if(record==3){
    #return all pve with correlation
    return(c((kin_R2-R2.all),R2.all))
  }
}
var.ex_cubic <- function(snp,bed_m2=bed_m2,ph,kin,gwaadata=data_cubic,recode=F){
  recode=recode
  if(!require(hglm))
    require(hglm)
  # creat Z
  na <- complete.cases(cbind(ph,bed_m2))
  Z <- c.z.hglm(kin = kin[na,na])
  # fit only IBS
  geno <- bed_m2[,snp]
  onlyIBS.hglm <- hglm(X = matrix(rep(1, length(gwaadata@gtdata@idnames[na]))), y = ph[na], Z = Z)
  cat("fitting Null model","\n")
  # fit snp model
  
  withsnp.hglm <- hglm(X = cbind(1, geno[na,]), y = ph[na], Z = Z)#geno就是marker，不接受NA
  cat("fitting SNP model","\n")
  r1 <- onlyIBS.hglm$varRanef/(onlyIBS.hglm$varRanef+ onlyIBS.hglm$varFix)#个体间差异对于表型的影响
  r2 <- withsnp.hglm$varRanef/(withsnp.hglm$varRanef+ withsnp.hglm$varFix)#排除了这些SNP后个体间差异对表型的影响
  #variation by snp ***
  R2.all <- 1- var(ph[na] - cbind(1, geno[na,])%*%withsnp.hglm$fixef)/var(ph[na])#SNP对于表型的影响
  R2 <- numeric(length(snp))
  for( i in 1:length(snp)){
    R2[i] <- 1- var(ph[na] - cbind(1, geno[na,snp[i]])%*%withsnp.hglm$fixef[c(1,i+1)])/var(ph[na])
  }
  names(R2) <- snp
  # phe by kinship
  r2.tot <- r2*(1-R2.all)
  r.h2 <- (r1-r2.tot)/r1
  r.tot <- (r1-r2.tot)/R2.all
  eff <- withsnp.hglm$fixef
  s.g <- summary(withsnp.hglm)
  p <- s.g$FixCoefMat[,4]
  sd <- s.g$FixCoefMat[,2]
  eff[2:length(snp)] <- eff[2:length(snp)]
  names(eff)[1] <- "u" 
  if (recode==1){
    return(list("Var.kin"=r.h2,"h2"=r1,"Var.add"=r.tot,"p"=p,"eff"=eff,"sd"=sd,"coding"= as.character(geno.re$coding),"geno"=geno,"R2"=R2,"R2.all"=R2.all,"model.snp"=withsnp.hglm))
  }else if (recode==2){
    return(list("Var.kin"=r.h2,"h2"=r1,"Var.add"=r.tot,"p"=p,"sd"=sd,"eff"=eff,"R2"=R2,"R2.all"=R2.all,"model.snp"=withsnp.hglm,"geno"=geno,"Z"=Z,"na"=na))
  }else if (recode==3){
    return(c(R2.all,R2))
  }
}
###load data
load("C:\\Users\\MSI\\OneDrive\\桌面\\花期可塑性\\工作数据\\genabel文件和QTLs\\FT&QTL_data.RData")
###variance of set
##kinship
ibs=ibs(data_ft,weight="freq")
##data initialise
phe_all=phdata(data_ft)
geno=as.double.gwaa.data(data_ft)
na1=complete.cases(geno)
memory.limit(16000000)

pve_all=data.frame(matrix(nrow=29,ncol=3))
pb=txtProgressBar(style=3)
for(i in 3:31){
  setTxtProgressBar(pb,(i-2)/length(3:31))
  result=var.ex(snp=qtl_all,geno = geno,ph=phe_all[,i],kin=ibs,gwaadata=data_ft,record=3)
  pve_all[(i-2),]=c(names(phe_all)[i],result)
}
close(pb)

###visulization

names(pve_all)=c("trait","kinship","QTL")
pve_all$trait[11:29]=paste("bio",1:19)
R2_all=reshape2::melt(pve_all,id.vars="trait")
R2_all$trait=factor(R2_all$trait,levels=pve_all$trait)
R2_all$value=as.numeric(R2_all$value)
ggplot(R2_all,aes(x=trait,y=value))+
  geom_bar(aes(fill=factor(variable,levels=c("kinship","QTL"))),stat="identity", width = 0.6,position = "stack")+
  ylim(0,1)+xlab("")+ylab("h2")+guides(fill=guide_legend(title="Type"))+
  theme_bw()+theme(panel.grid=element_blank())+theme_classic()


