#!/usr/bin/env Rscript
#example usage
#Rscript ../algs/vamf_script.R -i data/MGH26/batch_separability/subsets_scran/test.tsv -o data/MGH26/batch_separability/vamf_out/test.tsv

library(modules)
#library(optparse)
import_package("optparse",attach=TRUE)

opt_list<-list(
  make_option(c("-i","--infile"), type="character", default=NULL, help="Required. Filename containing the sparse matrix to analyze in three column COO format (row_id, col_id, data_value). First row assumed to be header. Assumes observations are columns are rows are features."),
  make_option(c("-o","--outfile"), type="character", default=NULL, help="Required. Filename for output of the inferred latent factors. A header is included by default"),
  make_option(c("-d","--dims"), type="integer", default=2, help="Desired number of latent dimensions. Default is 2. VAMF can learn the dimensionality using automatic relevance determination. For this, specify a large number and check the norms of each dimension."),
  make_option(c("-p","--parallel"), type="integer", default=3, help="Number of parallel cores used to run replicate variational bayes procedures. Sometimes helps if one VB run converges to a local optimum or crashes. VAMF will only save the instance with the highest ELBO value"),
  make_option(c("-s","--svmult"),type="character", default="1", help="Hyperparameter multiplier(s) for the ARD prior, if more than one, provide as comma separated list eg '.5, 1, 1.5'. Larger values mean less regularization of the latent factors. Smaller values cause ARD to aggressively prune off extra dimensions")
)
opt_parser<-OptionParser(option_list = opt_list)
opt <- parse_args(opt_parser)
if(is.null(opt$infile)){
  print_help(opt_parser)
  stop("Input file required.")
} else if(is.null(opt$outfile)){
  print_help(opt_parser)
  stop("Output file required.")
} else {
  import_package("Matrix",attach=TRUE)
  vamf<-import("../algs/vamf_stan") #our algorithm
  #parallel<-import_package("parallel")
  #print(paste("number of cores=",parallel$detectCores()))
}

options(mc.cores=min(parallel::detectCores(),opt$parallel)) #use parallel cores

Y<-do.call(sparseMatrix,read.table(opt$infile,header=TRUE))
svmult<-as.numeric(strsplit(opt$svmult,",")[[1]])
res<-vamf$vamf(Y,opt$dims,nrestarts=opt$parallel,svmult=svmult)
write.table(res$factors,file=opt$outfile,row.names=FALSE,quot=FALSE)