library(minpack.lm)
#cidr_vc_dissim<-function(){}
logit<-function(p){log(p)-log(1-p)}
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

Y<-2^as.matrix(Y_obs)-1
sData <- cidr::scDataConstructor(Y)
sData <- cidr::determineDropoutCandidates(sData,zerosOnly=TRUE)
hist(sData@dThreshold) #these are the lower bounds
res<-wthresh(sData)
sData<-res$cidr_object
threshc<-res$cell_thresholds
hist(threshc)
abline(v=sData@wThreshold)
plot(threshc,sData@dThreshold)

#sData <- cidr::wThreshold(sData) #replace with custom function
sData <- cidr::scDissim(sData) #replace with custom function

sData <- cidr::scPCA(sData)
sData <- cidr::nPC(sData)
factors<-as.data.frame(sData@PC[,1:L])