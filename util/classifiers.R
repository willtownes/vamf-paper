# Functions for classifying data and visualizing decision boundaries
# some funcs based on http://stackoverflow.com/questions/24260576/plot-decision-boundaries-with-ggplot2

library(modules)
MASS<-import_package("MASS")
gp<-import_package("ggplot2")
svm<-import_package("e1071")

rm_const_cols<-function(X,eps=1e-15){
  #X is a matrix. 
  #This function finds the columns of X that are constant across rows
  #returns a matrix with only the non-constant columns
  #useful for fixing errors in LDA and QDA since they don't work with constant cols
  #eps is numerical threshold for something being essentially zero
  ranges<-apply(X,2,range)
  zranges<-abs(ranges[2,]-ranges[1,]) < eps
  as.matrix(X[,!zranges])
}

fit_classifier<-function(X,labels,cmethod=c("lda","qda","svm")){
  #using the two columns of X, fit classifier to labels
  #methods possible include 
  #return the classifier object fitting
  labels<-factor(labels)
  cmethod<-match.arg(cmethod)
  if(cmethod=="svm"){
    cmethod<-svm$svm
    #X<-as.matrix(X)
  } else {
    cmethod<-MASS[[cmethod]]
  }
  tryCatch(cmethod(X,labels),
           error=function(e){
             lda_nonconst<-grepl("appears to be constant within groups",e$message,fixed=TRUE)
             qda_nonconst<-grepl("rank deficiency in group",e$message,fixed=TRUE)
             if(!(lda_nonconst || qda_nonconst)){
               stop(e) #unknown error
             }
             #otherwise, case where X has constant values in some column
             #remedy: run LDA/QDA etc on all non-constant columns
             X2<-rm_const_cols(X)
             if(0 %in% dim(X2)) return(NULL) #all columns were constant, can't fit
             cmethod(X2,labels)
           })
}
error_rate<-function(fit,X,labels){
  #evaluate the misclassification error (0/1) of the classifier "fit"
  #X is a 2-column matrix of predictors
  #labels is the outcome data
  #returns NA if the original fit didn't work (fit==NULL)
  if(is.null(fit)) return(NA)
  preds<-tryCatch(predict(fit,X),
                  error=function(e){
                    if(!grepl("wrong number of variables",e$message,fixed=TRUE)){
                      stop(e) #unknown error
                    }
                    #otherwise case where X has constant values in some column
                    predict(fit,rm_const_cols(X))
                  })
  if(class(fit) != "svm") preds<-preds$class
  conf<-table(preds,factor(labels)) #confusion matrix
  1-sum(diag(prop.table(conf))) #misclass. rate
}

cls_factory<-function(cmethod=c("lda","qda")){
  #produces wrapper functions for running LDA/QDA and getting error rates
  cmethod<-match.arg(cmethod)
  fn<-function(X,labels){
    fit<-fit_classifier(X,labels,cmethod)
    error_rate(fit,X,labels)
  }
  return(fn)
}
lda_wrap<-cls_factory("lda")
qda_wrap<-cls_factory("qda")

vec2seq<-function(x,len=30){seq(from=min(x),to=max(x),length=len)}

plot_classifier<-function(fit,X,labels,pointsize=1){
  #only works for LDA, QDA
  if(class(fit)=="svm") stop("Plotting not supported for SVM")
  cnames<-colnames(X)
  X<-as.data.frame(X)
  group<-labels
  dcon<-expand.grid(data.frame(lapply(X,vec2seq)))
  dcon$probs<-predict(fit,dcon)$posterior[,1] #assume only two class labels
  gp$ggplot(X,gp$aes_string(x=cnames[1],y=cnames[2])) + gp$geom_point(size=pointsize,gp$aes(colour=group)) + gp$geom_contour(data=dcon,gp$aes_string(x=cnames[1],y=cnames[2],z="probs"),breaks=0.5)
}
