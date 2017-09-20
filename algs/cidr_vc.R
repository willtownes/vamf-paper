library(minpack.lm)
#cidr_vc_dissim<-function(){}
logit<-function(p){log(p)-log(1-p)}
expit<-function(x){1/(1+exp(-x))}
wthresh<-function(object,cutoff=.5){
  delete <- which(rowSums(object@dropoutCandidates)==object@sampleSize)
  if(length(delete)>0){
    nData <- object@nData[-delete,]
    dropoutCandidates <- object@dropoutCandidates[-delete,]
  } else {
    nData <- object@nData
    dropoutCandidates <- object@dropoutCandidates
  }
  
  N <- object@sampleSize
  dropoutRates <- rowSums(dropoutCandidates)/N
  nzmean <- function(x){mean(x[x!=0])}
  averLcpm <- apply(nData*!dropoutCandidates, 1, nzmean)
  qu <- nlsLM(dropoutRates ~ 1/(1+exp(a*(averLcpm-b))),
              start=list(a=1, b=round(median(averLcpm))), trace=FALSE)
  a <- coef(qu)[1]
  b <- coef(qu)[2]
  threshold <- 1/a*log(1/cutoff-1)+b
  #include cell-specific thresholds
  pdet<-1-colMeans(dropoutCandidates) #detection rates for each cell
  aec<-colMeans(nData) #average expression for each cell
  #aec<-apply(nData*!dropoutCandidates,2,nzmean)
  bn<-aec + logit(pdet)/a #logit shift parameters for each cell's curve
  threshc<-bn - logit(cutoff)/a #cell-specific thresholds for 50% detection
  #object@cThreshold <- threshc
  #object@dropoutCoefB_cells <- bn
  #end new part
  object@wThreshold <- threshold
  object@pDropoutCoefA <- a
  object@pDropoutCoefB <- b
  #return(object)
  return(list(cidr_object=object,cell_thresholds=threshc))
}

dist_inner<-function(x,y,xdropout,ydropout,xthresh,ythresh){
  #x,y are vectors representing two different cells
  #xdropout,ydropout are vectors of indicators of "dropout candidates" in x,y
  #xthresh,ythresh are scalars indicating the expression level at which 50% detection probability.
  x_questionable<- !xdropout & x<xthresh
  y_questionable<- !ydropout & y<ythresh
  dvec<-x-y
  #dvec[x_dropout & y_dropout]<-0 #already zero if zeros are only dropouts
  dvec[xdropout & y_questionable]<- 0
  dvec[x_questionable & ydropout]<- 0
  #sum(abs(dvec)) #L1 metric
  sum(dvec^2) #L2 metric
}

dist_inner_wt<-function(x,y,xdropout,ydropout,px,py){
  #x,y are vectors representing two different cells
  #xdropout,ydropout are vectors of indicators of "dropout candidates" in x,y
  #px,py are vectors of probabilities of detection
  dvec<-x-y
  #dvec[x_dropout & y_dropout]<-0 #already zero if zeros are only dropouts
  dvec[xdropout]<- dvec[xdropout]*py[xdropout]
  dvec[ydropout]<- dvec[ydropout]*px[ydropout]
  dvec[xdropout & ydropout] <- 0
  #sum(abs(dvec)) #L1 metric
  sum(dvec^2) #L2 metric
}

dist_expr_only<-function(x,y,xdropout,ydropout){
  good_genes<-!xdropout & !ydropout
  dvec<-x[good_genes] - y[good_genes]
  sum(dvec^2)
}

scdist<-function(object,threshc,distance=c("cutoff","weighted","expressed_only")){
  #threshc is vector of length Ncells with expresison levels at which 50% detection prob.
  distance<-match.arg(distance)
  N <- ncol(object@nData) #number of cells
  D <- array(0, dim=c(N, N))
  pdetect<-expit(object@pDropoutCoefA*(t(t(object@nData)-threshc)))
  for(i in 1:(N-1)){
    x<-object@nData[,i]
    xdropout<-object@dropoutCandidates[,i]
    px<-pdetect[,i]
    for(j in (i+1):N){
      y<-object@nData[,j]
      ydropout<-object@dropoutCandidates[,j]
      if(distance=="cutoff"){
        D[j,i]<-dist_inner(x,y,xdropout,ydropout,threshc[i],threshc[j])
      } else if(distance=="weighted"){
        D[j,i]<-dist_inner_wt(x,y,xdropout,ydropout,px,pdetect[,j])
      } else if(distance=="expressed_only"){
        D[j,i]<-dist_expr_only(x,y,xdropout,ydropout)
      }
    }
  }
  D <- sqrt(D)
  D <- D+t(D)
  object@dissim <- D
  return(object)
}

Y<-2^as.matrix(Y_obs)-1
sData <- cidr::scDataConstructor(Y)
sData <- cidr::determineDropoutCandidates(sData,zerosOnly=TRUE)
hist(sData@dThreshold) #these are the lower bounds
#sData <- cidr::wThreshold(sData) #replace with custom function
#sData <- cidr::scDissim(sData) #replace with custom function
res<-wthresh(sData)
sData<-res$cidr_object
threshc<-res$cell_thresholds
sData<-scdist(sData,threshc,distance="expressed_only")
sData <- cidr::scPCA(sData)
sData <- cidr::nPC(sData)
factors<-as.data.frame(sData@PC[,1:2])
colnames(factors)<-paste0("dim",1:2)
#plt_noise(factors,batch)
plt_clst(factors,ref)

#hist(threshc)
#abline(v=sData@wThreshold)
#plot(threshc,sData@dThreshold)