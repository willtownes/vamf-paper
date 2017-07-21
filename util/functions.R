library(modules)
import_package("ggplot2",attach=TRUE)
import_package("mgcv",attach=TRUE)
import_package("scam",attach=TRUE)
survival<-import_package("survival")
plyr<-import_package("plyr")
reshape2<-import_package("reshape2")
AnnotationDbi<-import_package("AnnotationDbi")
misc<-import_package("miscTools")
foreach<-import_package("foreach")
Matrix<-import_package("Matrix")
colorRamps<-import_package("colorRamps")
#parallel<-import_package("parallel")

ssNorm <- function(cnt) {
  cntPlus <- cnt
  cntPlus[cnt==0] <- 1
  alpha <- scaleFactor(cnt)
  cntNorm <- sweep(cntPlus, 2, alpha, FUN="/")
  censored <- cnt==0
  list(cntNorm=cntNorm, censored=censored)
}

scaleFactor <- function(cnt, plot=FALSE) {
  # Compute target
  libSize <- colSums(cnt)
  vl <- t(log2(t(cnt + 0.5)/(libSize + 1) * 1e+06)) # Same as voom E slot
  target <- Matrix$rowMeans(vl)
  binCutoffs <- c(-Inf, seq(0,10), Inf)
  targetBin <- cut(target, binCutoffs, right=FALSE)
  table(targetBin)
  l <- log2(cnt+0.5)
  med <- apply(l, 2, tapply, targetBin, median)
  detected <- med > (-1)
  if (plot==TRUE) {
    df <- reshape2$melt(med, varnames = c("Log2Count", "sampleid"))
    df$class <- pd$class[match(df$sampleid, pd$sample_id)]
    p <- ggplot(df, aes(Log2Count, value, group=sampleid)) + geom_line() + facet_wrap(~class) + theme_bw() + ggtitle("Median expression")
    plot(p)
  }
  detectionLimitCat <- apply(detected, 2, function(x) min(which(sapply(1:length(x), function(i) all(x[i:length(x)])))))
  detectionLimit <- binCutoffs[detectionLimitCat]
  subThreshold <- sapply(detectionLimit, function(x) target<x)
  subThreshold[is.na(subThreshold)] <- TRUE
  lr <- sweep(l, 1, target, FUN="-") # Sample to target log ratio
  lr[subThreshold] <- NA
  2^misc$colMedians(lr, na.rm=TRUE)
}

meanOrderNumber <- function(x, cens) {
  x <- -x
  o <- order(x)
  xO <- x[o]
  censO <- cens[o]
  n <- length(x)
  pmon <- 0 # Previous mean order number
  meanOrderNumber <- foreach$foreach(i=1:length(x), .combine=c) %do% {
    if (censO[i]==1) {
      mon <- pmon + (n-i+2)/2
    } else {
      mon <- pmon + (n+1-pmon)/(1+n-i+1)
      pmon <- mon
    }
    mon
  }
  ties <- unique(xO[duplicated(xO)])
  if (length(ties)>0) message("Warning: Ties")
  ret <- rep(NA, length(x))
  ret[o] <- meanOrderNumber
  ret
}

meanOrderNumberSS <- function(x, cens) {
  x <- -x
  o <- order(x)
  xO <- x[o]
  censO <- cens[o]
  n <- length(x)
  pmon <- 0 # Previous mean order number
  pmonc <- 0 # Previous censored mean order number
  cntNoncens <- 0 # Number of noncensored numbers since last censored value
  meanOrderNumber <- foreach$foreach(i=1:length(x), .combine=c) %do% {
    if (censO[i]==1) {
      if (sum(censO[1:i-1])==0) { #get the first pmonc (Martin's original algorithm)
      mon <- pmon + (n-i+2)/2
      pmonc <- mon
      }
      else { #algorithm for subsequent censored values
        mon <- pmonc + (cntNoncens)*(n+1-pmonc)/((1+n-i+1)+cntNoncens)
        pmonc <- mon
      }
      cntNoncens <- 0
    } else { #algorithm for non-censored values
      mon <- pmon + (n+1-pmon)/(1+n-i+1)
      pmon <- mon
      cntNoncens <- cntNoncens + 1
    }
    mon
  }
  ties <- unique(xO[duplicated(xO)])
  if (length(ties)>0) message("Warning: Ties")
  ret <- rep(NA, length(x))
  ret[o] <- meanOrderNumber
  ret
}

censRank <- function(x, cls, cens, rho=0) {
  df <- length(unique(cls))-1
  foreach$foreach(i=1:nrow(x), .combine="c") %dopar% {
    fit <- survival$survdiff(Surv(-x[i,], !cens[i,]) ~ cls, rho=rho)
    pchisq(fit$chisq, df,lower.tail=FALSE)
  }
}

################# Functions created by Will ####################
rm_bad_genes<-function(dat,threshold=5,log2_transform=TRUE){
  # remove bad genes from the dataset and return a cleaned-up version.
  # assumes the bad cells were already discarded
  detected_genes<-Matrix$rowSums(dat>0)
  dat<-dat[detected_genes >threshold,] #genes must show up in at least 5 cells
  if(log2_transform){
    return(log2(1+dat))
  } else {
    return(dat)
  }
}

tscale<-function(x){
  # apply scaling (subtract mean and divide by standard deviation) to rows instead of columns
  t(scale(t(x))) #scale by rows
}

get_cell_detection<-function(y){
  #calculate the number of detected genes in each column
  #apply(y,2,function(x){sum(x>0)})
  Matrix$colSums(y>0)
}

plot_qc<-function(meta,total_reads_cutoff=500000,align_cutoff=.2,gene_cutoff=2000,batch_col="sample"){
  #generate a list of ggplot objects to visualize the QC thresholds for cells and genes
  ggbase<-ggplot(data=meta)+theme_bw()+facet_wrap(batch_col)
  plt<-list()
  plt[["qc1_total_reads"]]<-ggbase+geom_histogram(aes(x=total_reads))+geom_vline(xintercept=total_reads_cutoff,colour="red")+ggtitle("Total Reads per Cell")+ylab("Number of Cells")

  plt[["qc2_pct_aligned_reads"]]<-ggbase+geom_histogram(aes(x=pct_aligned))+ ggtitle("Percent of Reads Aligned per cell")+ ylab("Number of Cells")+ geom_vline(xintercept=align_cutoff,colour="red")

  plt[["qc3_detection_hist"]]<-ggbase+geom_histogram(aes(x=detection))+geom_vline(xintercept=gene_cutoff,colour="red")+ggtitle("Detection Rates: #genes nonzero per cell")+ylab("Number of Cells")

  plt[["qc4_detection_sampleid"]]<-ggplot(data=meta)+theme_bw()+ geom_point(aes_string(x=1:nrow(meta),y="detection",colour=batch_col))+ ylab("Detection Rates (number of genes)")+ xlab("cell id")+ ggtitle("Detection Rates across Samples")+ geom_hline(yintercept=gene_cutoff,colour="red")

  plt[["qc5_combined_scatter"]]<-ggbase+geom_point(aes(x=total_reads,y=pct_aligned,colour=detection),size=5)+geom_hline(yintercept=align_cutoff,colour="red")+geom_vline(xintercept=total_reads_cutoff,colour="red")+ggtitle("Total Reads vs Percent Aligned")+scale_colour_gradient(low="green",high="black")

  plt[["qc6_detection_vs_total_reads"]]<-ggplot(data=meta)+theme_bw()+geom_point(aes_string(x="total_reads",y="detection",colour=batch_col))+xlab("Total reads per cell")+ylab("Number of Genes Detected")+ggtitle("Detection Rates Vary with Sequencing Depth")+geom_hline(yintercept=gene_cutoff,colour="red")+geom_vline(xintercept=total_reads_cutoff,colour="red")

  return(plt)
}

apply_qc<-function(dat,qc_vals,total_reads_cutoff=500000,align_cutoff=.2,gene_cutoff=2000,verbose=TRUE){
  # apply QC filters to the matrix "dat" (could be tpm or counts) with rows=genes, cols=cells.
  # qc_vals is a data frame with metadata. Should have ID column "cell_id" with ids matching exactly to the columns of dat. Also should have columns for total_reads and aligned_reads
  # filter criteria:
  # 1. more than 500,000 reads total
  # 2. more than 20% reads aligned successfully
  # 3. more than 2000 genes/transcripts detected
  # all genes/transcripts with zeros across all samples are removed for data compression
  stopifnot(all(colnames(dat) %in% qc_vals$cell_id))
  initial_cells<-ncol(dat)
  rownames(qc_vals)<-qc_vals$cell_id
  #re-sort so rows of qc_vals align with column order of dat
  qc_vals<-qc_vals[colnames(dat),]
  qc1<-qc_vals[,"total_reads"]>total_reads_cutoff
  #print(sum(qc1))
  #if(!("pct_aligned" %in% colnames(qc_vals))){
  #  qc_vals$pct_aligned<-with(qc_vals,aligned_reads/total_reads)
  #}
  qc2<-with(qc_vals,aligned_reads>align_cutoff*total_reads)
  #print(sum(qc2))
  qc3<-get_cell_detection(dat)>gene_cutoff
  #print(sum(qc3))
  qc_vals$passed_qc<-qc1&qc2&qc3
  dat<-dat[,qc_vals$passed_qc] #remove bad cells
  #remove zero genes/transcripts
  det_genes<-apply(dat,1,function(x){any(x>0)})
  dat<-dat[det_genes,]
  if(verbose) print(paste(ncol(dat),"cells passed QC out of",initial_cells))
  mget(c("dat","qc_vals"))
}

cleanup_transform<-function(x,batch_cols="batch",gene_thresh=3,log_trans=TRUE){
  #x is output from apply_qc function
  #x has elements "dat" and "qc_vals"
  #remove all "bad" cells/genes
  #transform data matrix to log2(1+Y) scale
  meta<-with(x,droplevels(qc_vals[qc_vals$passed_qc==TRUE,]))
  if(is.null(batch_cols) || length(batch_cols)==0){
    batch<-NULL
  } else {
    batch<-as.data.frame(meta[,batch_cols])
    colnames(batch)<-batch_cols
  }
  Y<-x$dat[,meta$cell_id]
  Z<-Y>0
  Zg<-Matrix$rowSums(Z)
  #hist(Zg)
  #sum(Zg>3) #number of genes appearing in >3 good cells
  Y<-Y[Zg>gene_thresh,] #dim(Y)=26487 x 69
  Z<-Z[Zg>gene_thresh,]
  #sum(Y>0)/prod(dim(Y))) #about 55% of values are observed.
  #log-transform
  Y[Z]<-log2(1+Y[Z])
  mget(c("Y","meta","batch"))
}

getDetectionLimit <- function(cnt) {
  #same function as scaleFactor except returns detection limits instead of scale factors
  # Compute target
  libSize <- Matrix$colSums(cnt)
  vl <- t(log2(t(cnt + 0.5)/(libSize + 1) * 1e+06)) # Same as voom E slot
  target <- Matrix$rowMeans(vl)
  binCutoffs <- c(-Inf, seq(0,10), Inf)
  targetBin <- cut(target, binCutoffs, right=FALSE)
  table(targetBin)
  l <- log2(cnt+0.5)
  med <- apply(l, 2, tapply, targetBin, median)
  detected <- med > (-1)
  detectionLimitCat <- apply(detected, 2, function(x) min(which(sapply(1:length(x), function(i) all(x[i:length(x)])))))
  detectionLimit <- binCutoffs[detectionLimitCat]
  return(detectionLimit)
}

names2batch<-function(data_names,batch_ids){
  # convert the vector of data column names to the corresponding batch ID.
  # batch_id is a named vector whose names are a superset of the data column names
  # often batch_id can be obtained from the getBatch() function in confounders.R
  # returns a vector of length ncol(dat) containing the batch ID of each column
  # this is useful for labeling each cell with its batch after QC filtering, which tends to reduce the number of cells, preventing a simple rbind or cbind operation
  df_dat<-data.frame(cnames=data_names)
  df_batch<-data.frame(cnames=names(batch_ids),batch=batch_ids)
  df<-plyr$join(df_dat,df_batch,by="cnames",type="left")
  res<-df$batch
  names(res)<-df$cnames
  return(res)
}

# getCC<-function(){
#   ### Cell Cycle Annotated Genes
#   #cycle base
#   if(!file.exists("data/human_periodic.tsv")){
#     download.file("http://download.jensenlab.org/Cyclebase3/human_periodic.tar","data/human_periodic.tar")
#     untar("data/human_periodic.tar",exdir="data")
#     file.remove("data/human_periodic.tar")
#   }
#   dataCB<-read.table("data/human_periodic.tsv",header=TRUE,stringsAsFactors = FALSE)
#   dataCB<-dataCB[1:min(600,nrow(dataCB)),"gene"] #restrict to top 600 genes
#   #dataCB is in ENSEMBL protein ID format
#   p2entrez<-unlist(AnnotationDbi$intraIDMapper(dataCB,"HOMSA","ENSEMBLPROT","EG"))
#   cb<-unlist(AnnotationDbi$intraIDMapper(p2entrez,"HOMSA","EG","ENSEMBL"))
#   #cb in ENSEMBL gene format

#   #G.O.
#   ens_ids_cc <- scLVM$getEnsembl('GO:0007049',species="Hs")
#   return(list(cb=cb,go=ens_ids_cc))
# }

#try color-blind friendly pallete
#http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/#a-colorblind-friendly-palette

image.matrix<-function(m){
  #make a heat map of a matrix using ggplot
  #overides default plotting method
  #note the transpose operation to present matrix in its correct orientation (unlike image.default())
  ggplot(reshape2$melt(t(m)), aes(Var1,Var2, fill=value))+geom_raster()+theme_bw()+scale_fill_gradientn(colours=colorRamps$matlab.like2(10))
  #+scale_fill_continuous(low=low_color,high=high_color)
}

#test code for image.matrix
#m <- matrix(rnorm(20),5)
#image(m)

split_by_col<-function(x,labels){
  #x is a matrix
  #labels are groupings of the columns of x
  #returns a list of submatrices, could regenerate x by cbind() if labels are ordered
  rnames<-rownames(x)
  cnames<-colnames(x)
  if(is.null(cnames)) cnames<-seq.int(ncol(x))
  if(is.matrix(x)){
    matfn<-matrix #regular matrix
  } else {
    matfn<-Matrix::Matrix #sparse matrix
  }
  label_lengths<-as.list(table(labels))
  d<-split(Matrix::t(x),labels) #list of vectors, uses recycling
  cnames<-split(cnames,labels)
  #reshape vectors into matrices
  fn<-function(q,s,cn){
    res<-Matrix::t(matfn(q,nrow=s))
    rownames(res)<-rnames
    colnames(res)<-cn
    res
  }
  mapply(fn,d,label_lengths,cnames)
}

sparse2disk<-function(Y,filename){
  #given a sparse Matrix Y, write it to disk as COO format
  if(!is(Y,"sparseMatrix")) Y<-Matrix(Y,sparse=TRUE)
  write.table(Matrix$summary(Y), file=filename, row.names=FALSE, quot=FALSE)
}

# Smoothing Functions used in Shalek and Trapnell datasets

midpt<-function(x){
  # given a sequence of points, return sequence of midpoints between adjacent points, it will have length one less than x
  (x[1:(length(x)-1)]+x[-1])/2
}

binsmooth<-function(y,x,bins=25,xname="x",cname_base="y",value_name="avg_detection"){
  xqt<-unique(quantile(x,probs=seq(0,1,length.out=bins)))
  xmid<-midpt(xqt)
  cutpts<-cut(x,breaks=xqt,include.lowest=TRUE,right=FALSE)
  fn<-function(y){tapply(y,cutpts,mean)}
  if(is.null(dim(y))){ #case of vector
    z<-fn(y)
    res<-data.frame(x=xmid,y=z)
    colnames(res)<-c(xname,cname_base)
    return(res)
  } else if(length(dim(y))==2){ #case of matrix, apply to each column
    z<-apply(y,2,fn)
    res<-cbind(xmid,z)
    colnames(res)<-c(xname,paste0(cname_base,1:ncol(z)))
    res<-data.frame(res)
    res<-reshape2$melt(res,id.vars=xname,variable.name=cname_base,value.name=value_name)
    return(res)
  } else {
    stop("y must have dimension 1 or 2")
  }
}

mspline1<-function(x,y,k=10,lower=NA,upper=NA){
  #fits a monotonic spline to data
  #small values of k= more smoothing (flatter curves)
  #large values of k= more flexible (wiggly curves)
  #k is related to effective degrees of freedom and number of knots
  #use unconstrained gam to get rough parameter estimates
  #lower, upper optional bounds on the function
  #basically a slight modification of an example in the mgcv::pcls documentation
  dat<-data.frame(x=x,y=y)
  init_gam <- gam(y~s(x,k=k,bs="cr"))
  # Create Design matrix, constraints etc. for monotonic spline....
  sm <- smoothCon(s(x,k=k,bs="cr"),dat,knots=NULL)[[1]]
  mc <- mono.con(sm$xp,lower=lower,upper=upper) # monotonicity constraints
  M <- list(X=sm$X, y=y, #design matrix, outcome
            C=matrix(0,0,0), #equality constraints (none)
            Ain=mc$A, bin=mc$b, #inequality constraints
            sp=init_gam$sp, p=sm$xp, #initial guesses for param estimates
            S=sm$S, #smoothness penalty matrix
            w=y*0+1, off=0 #weights, offset
  )
  #fit spine using penalized constrained least squares
  p<-pcls(M)
  return(list(sm=sm,p=p))
}

predict.mspline1<-function(msp,x){
  #using the monotone spline msp, predict values for the vector x
  Predict.matrix(msp$sm,data.frame(x=x))%*%msp$p
}

mspline2<-function(x,Y,k=10,family=binomial()){
  #fits monotone increasing spline to pairwise comparisons between 
  #vector x and each column of Y separately
  #low value of k means more smoothing, high value more flexible curve
  #returns a data frame with columns:
  #"x" are span of 50 points with same range as original x
  #"f" is predicted values for each function
  #"df" is predicted derivative for each function
  #"id" is id_base and a number indicating the column of Y
  xo<-order(x)
  x<-x[xo]
  Y<-Y[xo,]
  if(is.null(colnames(Y))){
    cix<-seq.int(ncol(Y))
  } else {
    cix<-colnames(Y)
  }
  if(length(x)>50){
    xnz<-x[x>min(x)] #nonzero
    xct<-cut(xnz,50) #fifty break points
    xr<-tapply(xnz,xct,max) #find segment values actually present in data
    xr<-c(min(x),xr) #include "zero" point (special category)
  } else {
    xr<-x
  }
  #dat<-data.frame(x=x,Y)
  dat_list<-lapply(cix,function(j){data.frame(x=x,y=Y[,j],id=j)})
  fn<-function(dt,k,family,xr){
    fmbase<-paste0("~s(x,k=",k,",bs='mpi')")
    fit<-scam(as.formula(paste0("y",fmbase)),family=family,data=dt)
    res<-data.frame(dt,f=fitted(fit))#,df=derivative.scam(fit)$d) #only works with gaussian
    res$y<-NULL
    res<-res[res$x %in% xr,]
    res<-res[!duplicated(res),]
    res$df<-c(diff(res$f)/diff(res$x),0) #compute derivatives manually by finite difference
    res
  }
  res<-do.call("rbind",lapply(dat_list,fn,k,family,xr))#,mc.cores=ncores))
  res$id<-as.factor(res$id)
  res
}

inflection_func<-function(x,f,fstar=0.5){
  #given a sequence of points x and their (nondecreasing) function values f(x)
  #find the point xstar at which f crosses the fstar mark
  #useful for finding inflection point of monotone splines
  i<-which.max(x[f<=.5])
  x1<-x[i]; y1<-f[i]; x2<-x[i+1]; y2<-f[i+1]
  (.5-y1)*(x2-x1)/(y2-y1) + x1
}

summarize_mspline2<-function(dat,id="id"){
  #input: a data frame with columns x,y,id,f,df (such as output from mspline2 function
  #output: for each id value, minimum and maximum of f, maximum slope df,...
  #... and the point at which the function crosses the 0.5 mark (roughly the inflection point)
  #f must be monotonically increasing
  fn<-function(d){
    slope<-max(d$df[!is.nan(d$df)])
    lo<-min(d$f)
    hi<-max(d$f)
    inflection<-inflection_func(d$x,d$f)
    if(length(inflection)==0) inflection<-NA
    data.frame(slope,lo,hi,inflection)
  }
  plyr$ddply(dat,id,fn)
  #plyr$ddply(dat,id,plyr$summarise, slope=max(df[!is.nan(df)]), lo=min(f), hi=max(f), inflection=inflection_func(x,f))
}