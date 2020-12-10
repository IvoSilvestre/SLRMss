fitted.SLRMss <-
  function(object,H0=FALSE){
  if(H0==FALSE){
  fitted=object$y.fitted
  }else{
  fitted=object$X[,colnames(object$X)%in%rownames(object$beta.coefficients.h0)]%*%object$beta.coefficients.h0[,1]
  }
  return(fitted)
 }
