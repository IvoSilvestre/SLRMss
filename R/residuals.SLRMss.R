residuals.SLRMss <-
function(fit,H0=FALSE,std=FALSE){
   res=(fit$y-fited(fit,H0=H0))
   if(std==T){
   return(res/coef(fit,H0=H0))  
   }else{
   return(res)  
   }
}
