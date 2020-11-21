residuals.SLRMss <-
function(fit,H0=FALSE,std=FALSE){
  if(H0==FALSE){
    if(std==F){
    return(fit$residuals)
    }else{
    return(fit$std.residuals)
    }
  ]else{
    if(std==F){
    res=(fit$y-fited(fit,H0=T))
    return(res)
    }else{
    return(res/fit$phi.h0)
    }
  }

}
