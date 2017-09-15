# Functions for running existing methods such as tSNE, PCA, ZIFA, etc on censored data
library(modules)
#irlba<-import_package("irlba")
#MASS<-import_package("MASS")
import_package("Matrix",attach=TRUE)
#tsne_pack<-import_package("tsne")
#Rtsne<-import_package("Rtsne")
#simlr_pack<-import_package("SIMLR")
bioc_pcaMethods<-import_package("pcaMethods")

rm_zero_rowcol<-function(Y){
  #remove all rows and columns containing all zeros
  Y<-Y[rowSums(Y>0)>0,] #remove rows with zeros all the way across
  Y<-Y[,colSums(Y>0)>0]
  Y
}

norm<-function(v){sqrt(sum(v^2))}

colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  apply(x,2,norm)
}

pca_irlba<-function(Y,L=2,center=TRUE,scale=TRUE){
  #Y is data matrix, L is desired number of latent dimensions
  #scale is flag of whether or not to center and scale Y before applying irlba
  #Y<-as.matrix(Y)
  Y<-rm_zero_rowcol(Y)
  #if(scale) Y<-t(scale(t(Y)))
  #svd1<-irlba::irlba(Y,L)
  #factors<-t(svd1$d * t(svd1$v)) #using recycling
  factors<-irlba::prcomp_irlba(t(Y),n=L,center=center,scale=scale)
  #colnames(factors)<-paste0("pca_irlba",1:L)
  colnames(factors)<-paste0("dim",1:L)
  as.data.frame(factors$x)
}

pca<-function(Y,L=2,center=TRUE,scale=TRUE){
  Y<-rm_zero_rowcol(Y)
  #if(scale) Y<-scale(Y)
  factors<-prcomp(as.matrix(t(Y)),center=center,scale=scale)$x
  factors<-factors[,1:L]
  colnames(factors)<-paste0("dim",1:L)
  as.data.frame(factors)
}

mds<-function(Y,L=2,metric=TRUE,distance="euclidean",scale=TRUE){
  #Multidimensional Scaling
  #Y is data
  #L is desired latent dimension
  #metric=TRUE means use cmdscale(). metric=FALSE means use isoMDS()
  #see http://www.statmethods.net/advstats/mds.html
  #distance is choice of distance function passed to dist()
  Y<-rm_zero_rowcol(Y)
  Yt<-scale(t(Y),scale=scale)
  d <- dist(as.matrix(Yt),distance) # euclidean distances between the cols
  if(metric){
    fit<-cmdscale(d,k=L)
  } else {
    fit<-MASS::isoMDS(d,k=L)$points
  }
  colnames(fit)<-paste0("dim",1:L)
  as.data.frame(fit)
}

tsne<-function(Y,L=2,center=TRUE,scale=TRUE,method="Rtsne",...){
  Yt<-t(rm_zero_rowcol(Y))
  if(center || scale){
    Yt<-scale(Yt,center=center,scale=scale)
  }
  Yt<-as.matrix(Yt)
  if(method=="Rtsne"){
    fit<-Rtsne::Rtsne(Yt,dims=L,...)$Y
  } else if(method=="tsne"){
    fit<-tsne::tsne(Yt,k=L,...)
  }
  colnames(fit)<-paste0("dim",1:L)
  as.data.frame(fit)
}

simlr<-function(Y,nclust,L=2,...){
  res<-SIMLR::SIMLR(Y,nclust,no.dim=L,...)
  #the clustering results are in res$y
  factors<-res$ydata
  colnames(factors)<-paste0("dim",1:L)
  factors
}

ppca<-function(Y,L=2,...){
  stop("implementation not yet finished")
  pcaMethods::pp(Y,L,...)
}