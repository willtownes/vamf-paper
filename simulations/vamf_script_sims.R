# Run VAMF against the simulated data
# usage: Rscript vamf_script_sims.R 3
# will run the script on the data in simulation index 3

i <- as.integer(commandArgs(TRUE)[1])
#i<-1
NRESTARTS<-3
options(mc.cores=min(parallel::detectCores(),NRESTARTS))

library(modules)
vamf<-import("../algs/vamf_stan")

d1<-paste0("data/G500/vamf_out/",i)
d2<-paste0("data/G1000/vamf_out/",i)

dir.create(d1,showWarnings=FALSE,recursive=TRUE)
dir.create(d2,showWarnings=FALSE,recursive=TRUE)

vamf_idx<-read.csv("data/vamf_param_index.csv")

for(j in seq.int(nrow(vamf_idx))){
  vfunc<-function(Y){
    vamf$vamf(Y,vamf_idx[j,"L"],nrestarts=NRESTARTS,log2trans=FALSE,svmult=vamf_idx[j,"svmult"])
  }
  Y<-do.call(Matrix::sparseMatrix,read.table(paste0("data/G500/input/Y",i,".tsv"),header=TRUE))
  ref<-read.table(paste0("data/G500/input/ref",i,".tsv"),header=TRUE)
  res<-vfunc(Y)$factors
  colnames(res)<-paste("vamf",colnames(res),sep="_")
  write.table(cbind(ref,res),file=paste0(d1,"/",j,".tsv"),row.names=FALSE,quot=FALSE)
  
  Y<-do.call(Matrix::sparseMatrix,read.table(paste0("data/G1000/input/Y",i,".tsv"),header=TRUE))
  ref<-read.table(paste0("data/G1000/input/ref",i,".tsv"),header=TRUE)
  res<-vfunc(Y)$factors
  colnames(res)<-paste("vamf",colnames(res),sep="_")
  write.table(cbind(ref,res),file=paste0(d2,"/",j,".tsv"),row.names=FALSE,quot=FALSE)
}
