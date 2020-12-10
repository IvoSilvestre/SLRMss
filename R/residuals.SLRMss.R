residuals.SLRMss <-
function(object,H0=FALSE,std=FALSE){
   res=(object$y-fitted(object,H0=H0))
   if(std==T){
   return(res/coef(object,H0=H0)$phi)  
   }else{
   return(res)  
   }
}
