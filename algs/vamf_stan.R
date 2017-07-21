# Varying-censoring Aware Matrix Factorization (VAMF)
# inference via Stan

library(modules)
import_package("Matrix",attach=TRUE)
irlba<-import_package("irlba")
gdata<-import_package("gdata")
parallel<-import_package("parallel")
library(rstan)
rstan_options(auto_write=TRUE)

STMOD<-stan_model("../algs/vamf.stan")

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

effective_dimension<-function(x,thresh=.05){
  #x is a matrix representing latent factors in the columns.
  #finds the effective dimension by first computing L2 norms of each column
  #columns with L2 norm greater than the maximum norm times 'thresh' are nonzero
  #columns with L2 norm less than the cutoff (above) are 'zero'
  #returns a count of the number of nonzero columns
  l2norms<-colNorms(x)
  sum( l2norms > max(l2norms)*thresh )
}

init_ss<-function(Y,L,log2trans=TRUE,pseudocount=0.0,b1_range=c(0.3,.7)){
  #Y is data matrix, typically non-negative values such as counts or TPM
  #L is desired latent dimension (upper bound)
  #log2_trans= should nonzero values of Y be converted to log-scale?
  #transformation is g(y)=log2(pseudocount+y) for all nonzero y values
  #b1_range is a rough guess of typical slopes for the dropout mechanism
  Y<-Matrix(Y,sparse=TRUE)
  Y<-rm_zero_rowcol(Y) #remove rows with all zeros
  Z<-Y>0
  #convert sparse matrix to stan-friendly long format
  Ys<-summary(Y) #converts to triplet format
  if(log2trans) Ys[,3]<-log2(pseudocount+Ys[,3])
  #stan data variables & hyperparameters
  ss<-list(gg=Ys[,1],nn=Ys[,2],y=Ys[,3],nvals=nrow(Ys),L=L,N=ncol(Y),G=nrow(Y))
  ss$Z<-t(matrix(as.integer(as.matrix(Z)),nrow=ss$G)) #dense integer matrix for stan
  #note ss$Z is transpose of Z
  ss$ymn<-median(ss$y) #automatically calculated in stan program
  Yctr<-Y-ss$ymn*Z
  #a<-colSums(Yctr)/colSums(Z)
  #Yctr<-t(t(Yctr)-a*t(Z)) #take advantage of recycling
  w<-rowSums(Yctr)/rowSums(Z)
  #ss$sa<-mad(a)
  ss$sw<-mad(w)
  ### this block inefficient, improve later!
  Yctr<-as.matrix(Yctr)-w
  #sd_cols<-apply(Yctr,2,function(x){mad(x[x>0])})
  #ss$su<-mean(sd_cols[!is.na(sd_cols)])
  #rough estimate for variation in the row factors
  sd_rows<-apply(Yctr,1,function(x){mad(x[x>0])})
  ss$sv<-mean(sd_rows[!is.na(sd_rows)])
  ### end inefficient block
  ss$Q<-colMeans(Z) #detection rates
  ss$b1_mn<-mean(b1_range)
  ss$b1_sd<-diff(b1_range)/4 #ie, range is +/- 2 standard devs from mean
  ss
}

vb_wrap<-function(svmult,stmod,ss,resnames,output_samples){
  #convenience function for parallelizing calls to VB
  #Also pipes stdout to string for later parsing (to get ELBO, etc)
  #here svmult is a scalar not a vector
  ss$sv_rate<- 1/(svmult*ss$sv) #empirical bayes hyperparam set
  #since Gamma shape is 2, the mode is svmult*ss$sv
  logtxt<-capture.output(stan_vb_fit<-vb(stmod,data=ss,pars=resnames,output_samples=output_samples),type="output",split="true")
  mget(c("stan_vb_fit","logtxt"))
}

vamf_stan<-function(ss, svmult=rep.int(1.0,4), output_samples=100){
  # ss is a list of data and hyperparameters
  # length of svmult determines number of parallel restarts to run
  # svmult values are multiplier for estimate of sigma_v hyperparameter. Large value= less shrinkage
  # output_samples controls how many samples from variational distr are used to compute approximate posterior mean.
  resnames<-c("U","w","sy","y0","V","sv","b0","b1")
  #variational bayes
  res<-parallel$mclapply(svmult,vb_wrap,STMOD,ss,resnames,output_samples)
  #get rid of model runs that resulted in errors
  return(Filter(function(x){class(x)!="try-error"},res))
}

extract_elbo<-function(logtxt){
  #extract the ELBO value at the final iteration based on log file output
  #assumes logtxt is a list of strings (one element per line of output)
  elbo<- -Inf
  x<-grep("ELBO CONVERGED",logtxt,value=TRUE,fixed=TRUE)
  if(length(x)!=1){
    if(length(x)==0){
      warning("Failed to converge!")
    } else {
      warning("Multiple convergence points found!")
    }
    return(elbo)
  }
  x<-unlist(strsplit(x,"\\s+",perl=TRUE))
  elbo<-as.numeric(x[3])
  #parse out scientific notation eg -3e+04
  #x<-grep("\\de\\+\\d+",x,perl=TRUE,value=TRUE) 
  #elbo<-as.numeric(x)
  return(elbo)
}

extract_factors<-function(stan_vb_fit,ss){
  #stan_vb_fit<-vamf_stan_fit$stan_vb_fit
  resnames<-stan_vb_fit@sim$pars_oi
  resnames<-resnames[resnames != "lp__"]
  vars_tmp<-rstan::summary(stan_vb_fit)$summary[,"mean"]
  varnames<-names(vars_tmp)
  vars<-lapply(resnames,function(x){vars_tmp[gdata$startsWith(varnames,x)]})
  #vars<-lapply(resnames,function(x){get_posterior_mean(stan_vb_fit,x)})
  names(vars)<-resnames
  #stan output indexes first by rows, then cols (inconvenient for recycling)
  vars$U<-t(matrix(vars$U,ncol=ss$L)) #output dim: LxN
  vars$V<-matrix(vars$V,nrow=ss$L) #output dim: LxG
  vars
}

ortho_vamf<-function(vars){
  #convert factors to orthonormal basis
  v<-vars$V #LxG
  u<-vars$U #LxN, recycles the sv vector
  svd_v<-svd(v)
  A<-svd_v$u #LxL?
  D<-if(length(svd_v$d)>1) diag(svd_v$d) else svd_v$d #LxL?
  Q<-svd_v$v #GxL?
  loadings<-t(Q) #LxG?
  factors<-data.frame(crossprod(u,A%*%D))
  colnames(factors)<-paste0("dim",1:ncol(factors))
  c(vars,mget(c("factors","loadings")))
}

ortho_extract<-function(stfit,ss){
  ortho_vamf(extract_factors(stfit$stan_vb_fit,ss))
}

vamf<-function(Y, L, nrestarts=4, log2trans=TRUE, pseudocount=0.0, output_samples=100, save_restarts=FALSE,svmult=1){
  #convenience wrapper for running instances of selection factorization
  ss<-init_ss(Y,L,log2trans=log2trans,pseudocount=pseudocount)
  svmult_vals<-rep_len(svmult,nrestarts)
  stan_fits<-vamf_stan(ss, svmult_vals, output_samples=output_samples)
  elbos<-vapply(stan_fits, function(x){extract_elbo(x$logtxt)}, 1.0)
  if(save_restarts){
    factor_list<-lapply(stan_fits, ortho_extract, ss)
  } else { #keep only max elbo values
    good<- elbos==max(elbos)
    factor_list<-lapply(stan_fits[good], ortho_extract, ss)
    svmult_vals<-svmult_vals[good]
    elbos<-elbos[good]
  }
  eff_dim_table<-rep.int(0,L)
  eff_dim_idx<-rep.int(NA,length(elbos))
  for(i in seq_along(elbos)){
    factor_list[[i]]$elbo<-elbos[i]
    factor_list[[i]]$svmult<-svmult_vals[i]
    L_eff<-effective_dimension(factor_list[[i]]$factors)
    factor_list[[i]]$effdim<-L_eff
    eff_dim_table[L_eff]<-eff_dim_table[L_eff]+1
    eff_dim_idx[i]<-L_eff
  }
  if(save_restarts){
    return(factor_list)
  } else if(length(factor_list)==1){
    return(factor_list[[1]])
  } else { #if ELBOs are tied, choose result with most common effective dimension
    L_eff_popular<-which.max(eff_dim_table) #most popular effective dimension
    return(factor_list[eff_dim_idx==L_eff_popular][[1]])
  }
}
