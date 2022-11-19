# GWAs
This repo is a simple fast GWAs pipeline.
## Table of contents
*[Dependencies](README.md#dependencies)
*[Input data description](README.md#Input-data-description)
*[Dependencies](README.md#dependencies)
*[Dependencies](README.md#dependencies)
*[Dependencies](README.md#dependencies)
*[Dependencies](README.md#dependencies)
## Dependencies
Code can be run on win, and also can be run on linux by changing the file path and software to linux form.
### R(3.6<version<4.0)
    packages
    - **GenABEL** for quality control and finding population structure.
    - qqman for drawing the manhattan plot.
    - IntAssoPlot for draw locuszoom plot.
    - parallel for multiple thread calculation.
### Softwares
<font color=gray>These softwares all can do GWAs independently.</font>
#### [PLINK](https://www.cog-genomics.org/plink/)
Convert GenABEL form data to binary form data.
#### [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#fastGWA)
Use a fast MLM-based Genome-Wide Association method to do GWAs.
#### [TASSEL](https://tassel.bitbucket.io/)
Convert data vcf file to diplo-hmm file.

## Input data description
### gwwa.data-class
In [**GenABEL-package**](https://github.com/HanYu-me/GWAs/blob/main/GenABEL-tutorial.pdf), a special data class, gwaa.data-class is used to store GWA data. You may consider an object of gwaa.data-class as a ’black box’ from which
you can get specific data using specific functions.
### Phenotype data
 You must have a file containing these information separated by tabs organized as such:
 ```txt
 id     phe1    phe2
 1001   100     122
 1002   101     125
 ...    ...     ...

 ```
 Where the id must be the same as gwaa.data-class individual.

## Add phenotype to Genabel data
Use GenABEL function `add.phdata()` to add phenotype value into data.
## Population structure
Population structure was described as cmd distance.
```R
ibs_m=ibs(fl.data.merge,weight="freq")
dis_m=as.dist(0.5-ibs_m)
cmd=cmdscale(dis_m)
```
And use k-means cluster method can drop outlier
```R
#choose best k
wss <- (nrow(cmd)-1)*sum(apply(cmd,2,var))
for(i in 2:10){
  wss[i]=kmeans(cmd,centers=i)$withinss  
}
plot(1:10, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
```
![population structure](imgs/k-means_population_structure.png)
## Heritability
Use GCTA to calculate narrow sense heritability：
```R
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
```
And And the result file contains the information as below
```txt
Source	Variance	SE
V(G)	0.807125	0.083061
V(e)	0.230305	0.095955
Vp	1.037431	0.071641
V(G)/Vp	0.778004	0.085178
logL	-362.320
logL0	-593.985
LRT	463.329
df	1
Pval	0.0000e+00
n	504
```
## GWAs and visualization
Also, use GCTA to get GWAs result by fast MLM-based Genome-Wide Association method.
```R
GCTA_mlm <- function(gcta,bfile,grm,pheno,mlm,mark=0){
    cmd_grm=paste(gcta,"--bfile",bfile,"--make-grm --sparse-cutoff 0.05 --thread-num 10 --out",grm,sep=" ")
    cat(cmd_grm,"\n")
    system(cmd_grm)
}
```
The result file contains the information as below:
```
CHR	SNP	POS	A1	A2	N	AF1	BETA	SE	P
1	1:73	73	A	C	504	0.0615079	-0.0236617	0.406874	0.953625
1	1:92	92	C	A	504	0.40873	0.0986783	0.198852	0.619725
```
P is the raw independent p-value, which need to be corrected beofore  visualization. Bonferroni-correct is the most used method.

***
After get the result file, we can draw QQplot and manhattan plot:
```R
lam <- estlambda(data$P,plot = F) 
```
![QQplot](imgs/qqplot.png)
```R
qqman::manhattan(data,col = c("blue4", "orange3"))
```
![manhattan plot](imgs/manhattan.png)

## Locuszoom
I use R package `IntAssoPlot` to draw locuszoom plot. Diplo-hmp,gtf,and GWAs result file must be prepared in advance.
```R
gwa=read.csv(file.choose(),header=T,sep="\t",stringsAsFactors = F)
gtf=read.csv(file.choose(),header = F,sep="\t",stringsAsFactors = F)
hmp=read.csv(file.choose(),header = T,sep="\t",stringsAsFactors = F) 
IntRegionalPlot(chr=5,left=(26223729-10000),right=(26223729+10000),gtf=gtf,association=gwa,hapmap=hmp,hapmap_ld=hmp,threshold=5,leadsnp_size=2,label_gene_name = TRUE,marker2label = marker,marker2label_angle = 0,marker2label_size = 3)
```
![locuszoom](imgs/locuazoom.png)

## Parallel in R
Users can use R package `parallel` to accelerate any step cost long time.
```R
fun<-function(i){

}
clnum<-detectCores()
cl <- makeCluster(getOption("cl.cores", clnum));
clusterExport(cl,deparse(substitute(fun)))
clusterEvalQ(cl,library(anypackage))
clusterExport(cl,c("variable1","variable2"))#  must load the variables, that the function will used.
parLapply(cl, 3:31, fun)
stopCluster(cl)
```