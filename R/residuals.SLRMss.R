residuals.SLRMss <-
function(fit,H0=FALSE,std=FALSE){
   res=(fit$y-fitted(fit,H0=H0))
   if(std==T){
   return(res/coef(fit,H0=H0)$phi)  
   }else{
   return(res)  
   }
}
