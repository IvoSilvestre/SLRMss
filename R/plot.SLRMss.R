plot.SLRMss<-
  function(x,H0=FALSE,
             xlab=c("Fitted Values","Index","Theoretical Quantiles","Standardized Residuals"),
             ylab=c("Standardized Residuals", "Standardized Residuals","Sample Quantiles","Density"),
             main=c("Against Fitted Values","Against Index","Normal Q-Q Plot","Density Estimate")){
    par(mfrow=c(2,2))
    resid=residuals(x,std=T,H0=H0)
    plot(fitted(x),resid,
         xlab=xlab[1],ylab=ylab[1],
         main=main[1],pch=16)
    plot(resid,xlab=xlab[2],
         ylab=ylab[2],main=main[2],pch=16)
    qqnorm(resid,xlab=xlab[3],ylab=ylab[3],main=main[3],pch=16)
    abline(0,1,col="red")
    plot(density(resid),xlab=xlab[4],
         ylab=ylab[4],main=main[4])
    points(resid,rep(0,length(resid)),pch=16)
    par(mfrow=c(1,1))
  }
