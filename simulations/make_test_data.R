#script to produce simulated data for testing VAMF ARD dimension learning
#scenario I: N=160, G=500, Gsignal=50
#scenario II: N=160, G=1000, Gsignal=50

library(modules) #devtools::install_github('klmr/modules')
sims<-import("../simulations/sims") #simulate missing data mechanisms
#import_package("Matrix",attach=TRUE)
set.seed(101)
Nsims<-10

dir.create("data/G500/input",recursive=TRUE,showWarnings=FALSE)
dir.create("data/G1000/input",recursive=TRUE,showWarnings=FALSE)

for(i in seq.int(Nsims)){
  s<-sims$latent_clusters(160,500,50)
  write.table(Matrix::summary(s$Y_obs),file=paste0("data/G500/input/Y",i,".tsv"),row.names=FALSE,quot=FALSE)
  write.table(s$ref,file=paste0("data/G500/input/ref",i,".tsv"),row.names=FALSE,quot=FALSE)
  
  s<-sims$latent_clusters(160,1000,50)
  write.table(Matrix::summary(s$Y_obs),file=paste0("data/G1000/input/Y",i,".tsv"),row.names=FALSE,quot=FALSE)
  write.table(s$ref,file=paste0("data/G1000/input/ref",i,".tsv"),row.names=FALSE,quot=FALSE)
}
