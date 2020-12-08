fitted.SLRMss <-
  function(fit,H0=FALSE){
  if(H0==FALSE){
  fitted=fit$y.fitted
  }else{
  fitted=fit$X[,colnames(fit$X)%in%rownames(fit$beta.coefficients.h0)]%*%fit$beta.coefficients.h0[,1]
  }
  return(fitted)
 }
