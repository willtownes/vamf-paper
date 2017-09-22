#fast approximate version of VAMF based on coordinate ascent optimization

#input Y a matrix with genes in rows cells in cols
#L number of dimensions

logit<-function(p){log(p)-log(1-p)}

norm<-function(v){sqrt(sum(v^2))}

colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  apply(x,2,norm)
}

Y<-as.matrix(Y_obs)
L<-5
nIter<-100
eps<-.1 #gradient ascent step size
s2v<-1; s2u<-.1; s2w<-10;
N<-ncol(Y); G<-nrow(Y)
s2y<-sd(Y)^2
y0<-mean(Y)
V<-matrix(rnorm(L*G,0,sv),nrow=L)
V<-t(svd(V)$v)
U<-matrix(rnorm(L*N,0,su),nrow=L)
w<-rowMeans(Y)-y0
Z<-Y>0
Zrs<-rowSums(Z)
Yrs<-rowSums(Y)
#rough estimate for dropout parameters a,b
a<-.75
pdet<-colMeans(Z) #detection rates for each cell
aec<-colMeans(Y) #average expression for each cell
b<-aec + logit(pdet)/a #logit shift parameters for each cell's curve, a vector
#boxplot(b~ref$batch)
R<-crossprod(V,U)+y0+w #recycling of w vector
Lambda<-t(t(R)>b)
for(i in 1:nIter){
  #update U
  ZL<-Z-Lambda
  M<-V%*%(a*ZL+(1/s2y)*(Y-Z*(y0+w))) #LxN temporary
  penalty<-(1/s2u)*diag(L)
  for(n in 1:N){
    #Vtilde_n<-t(t(V)*Z[,n])
    U[,n]<-U[,n]-eps*(crossprod(t(V)*Z[,n])/s2y+penalty)%*%U[,n]+eps*M[,n]
  }
  zUv<-vapply(1:G,function(g){Z[g,]%*%crossprod(U,V[,g])},FUN.VALUE=1.0)
  y0<-y0-eps*(sum(Zrs)/s2y)*y0+eps*(a*sum(ZL)+(1/s2y)*(sum(Yrs)-sum(w%*%Z)-sum(zUv)))
  w<-w-eps*(Zrs/s2y+1/s2w)*w+eps*(a*rowSums(ZL)+(Yrs-y0*Zrs-zUv)/s2y)
  M2<-U%*%(a*t(ZL)+(1/s2y)*t(Y-Z*(y0+w))) #LxG temporary
  penalty2<-(1/s2v)*diag(L)
  for(g in 1:G){
    #Utilde_g<-t(t(U)*Z[g,])
    V[,g]<-V[,g]-eps*(crossprod(t(U)*Z[g,])/s2y+penalty2)%*%V[,g]+eps*M2[,g]
  }
  #orthogonalize V
  Vsvd<-svd(V)
  AD<-Vsvd$u %*% diag(Vsvd$d) #LxL
  V<-t(Vsvd$v) #LxG
  U<-crossprod(AD,U)
  #sort factors in decreasing order
  fo<-order(apply(U,1,norm),decreasing = TRUE)
  U<-U[fo,]
  V<-V[fo,]
  R<-y0+w+crossprod(V,U)
  s2y<-sd(Y-R*Z)^2
  Lambda<-t(t(R)>b)
}
rownames(U)<-paste0("dim",1:L)
plt_clst(t(U[1:2,]),ref)

