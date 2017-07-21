#quick tests for vamf.stan program
#assumes data file already created by alg_test.Rmd
library(ggplot2)
source("../algs/vamf_stan.R")

Y_obs<-do.call(sparseMatrix,read.table("data/noise_only.tsv",header=TRUE))
batch<-factor(rep(c("High Detection","Low Detection"),each=ncol(Y_obs)/2))

ss<-init_ss(Y_obs,4,log2trans=FALSE)
stan_fits<-vamf_stan(ss)
elbos<-vapply(stan_fits, function(x){extract_elbo(x$logtxt)}, 1.0)
stan_vb_fit<-stan_fits[[1]]$stan_vb_fit

factor_list<-lapply(stan_fits[elbos==max(elbos)], ortho_extract, ss)
#if(length(factor_list)<1) stop("all cores failed")
factor_list[[1]]#$factors

plot(factor_list[[1]]$factors[,1:2],col=as.integer(batch))
barplot(existing$colNorms(factor_list[[1]]$factors))
barplot(existing$colNorms(t(factor_list[[1]]$loadings)))
