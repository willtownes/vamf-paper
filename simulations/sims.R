# Create simulated single cell RNA-seq data to test out the various algorithms
# A key criterion is ability to recover batch effects even after censoring.
# Also want to see that the latent variables are correctly estimated.
# Data is generated from the true model.

library(modules)
import_package("Matrix",attach=TRUE) #for sparse matrix handling
#handy creation of dummy variables for multiple batch factors
import_package("caret",attach=TRUE) 
MASS<-import_package("MASS")

rm_zero_rowcol<-function(Y){
  #remove all rows and columns containing all zeros
  Y<-Y[rowSums(Y>0)>0,] #remove rows with zeros all the way across
  Y<-Y[,colSums(Y>0)>0]
  Y
}

genCluster<-function(n,id,mu=c(0,0),sigma=diag(2)){
  d<-data.frame(MASS$mvrnorm(n,mu,sigma))
  d$id<-id
  return(d)
}

simulate_data<-function(L,G,N,Gsignal=G,U=NULL,y0=10,std=list(a=2,w=2,u=2,v=1/sqrt(L))){
  #L is dimension of latent subspace
  #G is number of rows, N is number of columns
  #Gsignal is number of rows that have signal, rest are noise
  #U is optionally prespecified matrix of latent factors (MxN)
  #default if NULL: all columns will be in same batch and no batch effect
  #y0: overall mean (bias/intercept term)
  #std: list of standard deviation parameters for random effects
  #RETURN: a list of parameters:
  # U is MxN matrix of latent factors (ie column/ cell factors)
  # V is MxG "factor loading" matrix (ie row/ gene factors)
  # a is Nx1 column biases
  # w is Gx1 row biases
  # y0 is overall bias/intercept term
  # Y is "true" latent expression matrix (uncensored)
  # E[Y_ng] = y0+a_n+w_g+u_n'v_g
  stopifnot(Gsignal>0 && Gsignal<=G)
  a<-if(std$a==0) 0 else rnorm(N,0,std$a)
  w<-rnorm(G,0,std$w)
  if(is.null(U)) U<-matrix(rnorm(N*L,0,std$u),nrow=L,ncol=N)
  V<-matrix(rnorm(Gsignal*L,0,std$v),nrow=L,ncol=Gsignal)
  Y<-t(V)%*%U
  if(Gsignal<G){
    Ynoise<-matrix(0,nrow=G-Gsignal,ncol=N)
    Y<-rbind(Y,Ynoise)
  }
  #rows=genes,columns=cells
  Y<-Y+w #take advantage of recycling to add same value to each row
  Y<-t(t(Y)+a) #again taking advantage of recycling
  Y<-y0+Y
  mget(c("U","V","a","w","y0","Y")) #combine into list of parameters
}

censor_dat_mar<-function(dat,censor_rate,s_y=1){
  noise<-rnorm(prod(dim(dat)),0,s_y)
  Z<-rsparsematrix(nrow=nrow(dat),ncol=ncol(dat),1-censor_rate,rand.x=NULL)
  #binary "pattern" matrix
  (dat+noise)*Z
}

expit<-function(x){1/(1+exp(-x))}
inv_probit<-function(x){pnorm(x)}
inv_cauchit<-function(x){.5+atan(x)/pi}
INV_LINK_FNS<-list(logit=expit,probit=inv_probit,cauchit=inv_cauchit)
QUANTILE_FNS<-list(logit=qlogis,probit=qnorm,cauchit=qcauchy)

rbeta2<-function(N,mu,phi=max(1/mu,1/(1-mu))+0.1){
  #return beta distributed variates with a given mean and precision
  #If phi not provided, shape chosen such that there is maximal dispersion for the given mean while still perserving a unimodal distribution
  # large values of phi mean smaller variance (higher precision)
  # more details on this parameterization
  # https://cran.r-project.org/web/packages/betareg
  #stopifnot(mu>0 && mu<1)
  #guarantee both shape pars >1 ==> unimodal
  rbeta(N, phi*mu, phi*(1-mu))
}

# det2thresh<-function(n,ptrue,Q,k,b1){
#   #n is index of cell
#   #for other params, see censor_dat_mnar1()
#   #Q is detection rate (1-censor_rate)
#   #(1-k) is fraction of total data MCAR
#   #b1 is slope param for MNAR model (ie capture efficiency)
#   un<-ptrue$U[,n]
#   s2<-var(ptrue$w)+un%*%cov(t(ptrue$V))%*%un
#   term1<- b1*(ptrue$y0+ptrue$a[n]) 
#   term2<- sqrt(s2+b1^2)*qnorm(Q/k)
#   dn<-(term1+term2)/sqrt(s2)
#   dn
# }

censor_dat_mnar1<-function(ptrue,cens_rates=0.5,mcar=0.3,capture_eff=0.5,sigma_y=2){
  # apply censoring using following formula
  # ptrue is a list of true parameters from simulate_data()
  # cens_rates is overall rate of missing data, including MCAR and MNAR parts
  # mcar is the fraction of missing data that is MCAR
  # ie, data mcar as a fraction of total data = mcar*censor_rate = "(1-k)"
  # capture_eff is slope parameter for probit MNAR model. Typically >0
  # sigma_y is stdev of observed Y values (measurement error noise)
  # MNAR model determined as follows
  # P(missing|eta_ng) = k*inv_probit(dn + capture_eff * eta_ng)
  # k=1-mcar*censor_rate
  # eta_ng = y0+an+wg+un'vg (from simulate_data function)
  # analytically integrating over all gene specific parameters (wg and vg), we get a formula for marginal missingness for cell n. This is set equal to the censor_rate
  # invert this formula to get implied dn values for each cell
  # use average of these as the intercept in the probit model
  # this leads to heterogeneity of the censoring rates, roughly centered around the desired censor_rate
  # see also http://arxiv.org/abs/1603.06045v1
  N<-ncol(ptrue$Y)
  G<-nrow(ptrue$Y)
  Yt<-t(ptrue$Y)
  b1<-capture_eff
  Qn<-1-cens_rates #overall detection rates for each cell
  kn<-1-mcar*cens_rates #1-k = overall MCAR rates
  stopifnot(all(kn>Qn))
  #dn<-sapply(1:N,function(n){det2thresh(n,ptrue,Q,k,b1)})
  dng<-t(qnorm(Qn/kn)-b1*Yt)
  dn<-colMeans(dng)
  probs<-t(kn*pnorm(dn + b1*Yt)) #implicit vectorization
  U<-matrix(runif(N*G),nrow=G)
  Z<-Matrix(U<probs)
  #hist(colMeans(Z))
  #plot(Qn,colMeans(Z))
  noise<-matrix(rnorm(N*G,0,sigma_y),nrow=G)
  res<-(ptrue$Y+noise)*Z
}

noise_only<-function(N,G,Gsignal,y0=5,std=list(a=0,w=2,v=1/sqrt(2)),cens=c(.56,.95),mcar=0.1,capture_eff=1,sigma_y=1){
  #creates a dataset based on the noise-only model (two batches with variable detection rates)
  stopifnot(N%%2==0) #even number of cells required
  dat<-matrix(rnorm(N*2,0,1/sqrt(2)),N)
  colnames(dat)<-paste0("dim",1:2)
  ptrue<-simulate_data(2,G,N,Gsignal,U=t(dat),y0=y0,std=std)
  batch1<-1:floor(N/2)
  batch2<-(floor(N/2)+1):N
  cens_rates<-rep(NA,N)
  cens_rates[batch1]<-rbeta2(length(batch1),cens[1],phi=200)
  cens_rates[batch2]<-rbeta2(length(batch2),cens[2],phi=200)
  ref<-data.frame(dat,cens_rates=cens_rates)
  ref$batch<-factor(rep(c("high_detection","low_detection"),each=N/2))
  Y_obs<-censor_dat_mnar1(ptrue,cens_rates,mcar=mcar,capture_eff=capture_eff,sigma_y=sigma_y)
  list(Y_obs=rm_zero_rowcol(Y_obs),ref=ref)
}

latent_clusters<-function(N,G,Gsignal,y0=5,std=list(a=0,w=2,v=1/sqrt(2)),cens=c(.56,.95),mcar=0.1,capture_eff=1,sigma_y=1){
  #generate samples from the four clusters, two batches scenario
  stopifnot(N%%4==0)
  #sigma<-matrix(c(1.5,1.3,1.3,1.5),nrow=2)
  sigma<-matrix(c(1,0,0,1),nrow=2)
  mu<-list(c(-2.5,2.5),c(2.5,2.5),c(2.5,-2.5),c(-2.5,-2.5))
  dat<-do.call("rbind",lapply(1:4,function(x){genCluster(N/4,x,mu[[x]],sigma)}))
  colnames(dat)[1:2]<-paste0("dim",1:2)
  ref<-data.frame(dat[,1:2],id=as.factor(as.character(dat$id)))
  ptrue<-sims$simulate_data(2,G,N,Gsignal,U=t(as.matrix(ref[,1:2])),y0=y0,std=std)
  #randomly assign half of cells to a low-censoring regime and other half to high-censoring
  batch1<-seq(from=1,to=N,by=2)
  batch2<-seq(from=2,to=N,by=2)
  cens_rates<-rep(NA,N)
  cens_rates[batch1]<-rbeta2(length(batch1),cens[1],phi=200)
  cens_rates[batch2]<-rbeta2(length(batch2),cens[2],phi=200)
  ref$batch<-factor(rep(c("high_detection","low_detection"),N/2))
  ref$cens_rates<-cens_rates
  Y_obs<-censor_dat_mnar1(ptrue,cens_rates,mcar=mcar,capture_eff=capture_eff,sigma_y=sigma_y)
  Y_obs<-rm_zero_rowcol(Y_obs)
  list(Y_obs=Y_obs,ref=ref)
}

### Test/Visualize the Simulations
# batch_labs<-data.frame(tumor=factor(rep(c("1","2"),each=N/2)),
#                        plate=factor(rep(rep(c("A","B"),each=N/4),2)))
# p<-simulate_data(15,100,48,U=NULL,batch_labs=batch_labs)
# #apply censoring process- missing at random. Choose random 30% to not be censored
# p$Y_obs<-censor_dat_mar(p$Y,.7)
#Z_obs<-Y_obs!=0
#image.matrix(Y_obs) #still see batch effect across the columns
#compute column means treating zeros as NA
#plot(1:N,colSums(Y_obs)/colSums(Z_obs)) #obvious batch effect
#plot(1:G,rowSums(Y_obs)/rowSums(Z_obs)) #no pattern
#if(CACHE) save(p,batch_labs,file="sim_data.RData")