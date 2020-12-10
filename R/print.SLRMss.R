print.SLRMss<-
  function(x){
  cat("Call:\n")
  print(x$call)
  cat("\nBeta Coefficients:\n")
  print(round(x$beta.coefficients[,1],4))
  cat("\nPhi:\n") 
  print(round(x$phi[,1],4))
  cat("\nAICc:\n")
  print(x$AIC)
  cat("\nBIC:\n")
  print(x$BIC)
  invisible(x)
}
