print.SLRMss<-
  function(fit){
  cat("Call:\n")
  print(fit$call)
  cat("\nBeta Coefficients:\n")
  print(fit$beta.coefficients[,1])
  cat("\nPhi:\n") 
  print(fit$phi[,1])
  cat("\nAICc:\n")
  print(fit$AIC)
  cat("\nBIC:\n")
  print(fit$BIC)
}