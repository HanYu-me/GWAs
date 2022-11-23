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
  geno=geno[,snp]
  na <- complete.cases(cbind(ph,geno))
  withsnp.hglm <- hglm(X = cbind(1, geno[na,]), y = ph[na], Z = Z)
  r2 <- withsnp.hglm$varRanef/(withsnp.hglm$varRanef+ withsnp.hglm$varFix)
  #variation by snp
  R2.all <- 1- var(ph[na] - cbind(1, geno[na,])%*%withsnp.hglm$fixef)/var(ph[na])
  kin_R2=1-var(ph[na]-cbind(1, geno[na,])%*%withsnp.hglm$fixef-Z%*%withsnp.hglm$ranef)/var(ph[na])
  #return all pve with correlation
  return(c((kin_R2-R2.all),R2.all))
}
##use ibs function to caculate kinship of genabel data data_ft
ibs=ibs(data_ft,weight="freq")
##data initialise
phe_all=phdata(data_ft)
geno=as.double.gwaa.data(data_ft)
na1=complete.cases(geno)
setTxtProgressBar(pb,(i-2)/length(3:31))
#qtl_all is fixed snp list
result=var.ex(snp=qtl_all,geno = geno,ph=phe_all[,i],kin=ibs,gwaadata=data_ft,record=3)




