plot.SLRMss<-
  function(x){
    par(mfrow=c(2,2))
    resid=residuals(x,std=T)
    plot(fitted(x),resid,
         xlab="Fitted Values",ylab="Standardized Residuals",
         main="Against Fitted Values",pch=16)
    plot(resid,xlab="Index",
         ylab="Standardized Residuals",main="Against Index",pch=16)
    qqnorm(resid,main="Normal Q-Q Plot",pch=16)
    abline(0,1,col="red")
    plot(density(resid),xlab="Standardized Residuals",
         ylab="Density",main="Density Estimate")
    points(resid,rep(0,length(resid)),pch=16)
    par(mfrow=c(1,1))
  }
