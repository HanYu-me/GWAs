# GWAs
This repo is a GWAs simple pipeline.

## Dependencies
Code can be run on win, and also can be run on linux by changing the file path and software to linux form.
#### R(3.6<=version<4.0)
    packages
    - **GenABEL** for quality control and finding population structure.
    - qqman for drawing the manhattan plot.
    - IntAssoPlot for draw locuszoom plot.
    - parallel for multiple thread calculation.
#### [Plink](https://www.cog-genomics.org/plink/)
Convert GenABEL form data to binary form data.
#### [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#fastGWA)
Use a fast MLM-based Genome-Wide Association method to do GWAs.

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
 Where the id must be same as gwaa.data-class individual.

## Add phenotype to Genabel data
Use GenABEL function `add.phdata()` to add phenotype value into data.
## Population structure
Population structure was describe as cmd distance.
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



## Heritability

## GWAs and visualization
