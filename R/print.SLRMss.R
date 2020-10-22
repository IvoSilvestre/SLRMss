print.SLRMss<-
  function(fit){
  cat("Call:\n")
  print(fit$call)
  cat("\nBeta Coefficients:\n")
  print(round(fit$beta.coefficients[,1],4))
  cat("\nPhi:\n") 
  print(round(fit$phi[,1],4))
  cat("\nAICc:\n")
  print(fit$AIC)
  cat("\nBIC:\n")
  print(fit$BIC)
}