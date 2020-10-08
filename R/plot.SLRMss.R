plot.SLRMss <-
  function(fit){
    family=fit$family
    if(family=="Normal"){
      qfam = function(x,mu,sigma){qnorm(x,mu,sigma)}
      pfam = function(x,mu,sigma){pnorm(x,mu,sigma)}
      rfam = function(x,mu,sigma){rnorm(x,mu,sigma)}
    }else{
      if(family=="Student"){
        qfam = function(x,mu,sigma){sigma*qt(x,df=fit$xi)+mu}
        pfam = function(x,mu,sigma){pt((x-mu)/sigma,df=fit$xi)}
        rfam = function(x,mu,sigma){sigma*rt(x,df=fit$xi)+mu}
      }else{
        if(family=="Logistic"){
          qfam = function(x,mu,sigma){qlogis(x,mu,sigma)}
          pfam = function(x,mu,sigma){plogis(x,mu,sigma)}
          rfam = function(x,mu,sigma){rlogis(x,mu,sigma)}
        }else{
          qfam = function(x,mu,sigma){qnormp(x,mu,sigma,p=2/(fit$xi + 1))}
          pfam = function(x,mu,sigma){pnormp(x,mu,sigma,p=2/(fit$xi + 1))}
          rfam = function(x,mu,sigma){rnormp(x,mu,sigma,p=2/(fit$xi + 1))}
        }
      }
    }
      set.seed(2612)
      y=fit$y
      n=length(y)
      J=100
      fv=fit$beta.fitted
      sige=fit$phi[1,1]
      rqobs <- qfam(pfam(fit$std.residuals, 0, 1),mu=0,sigma=1)
      mrq <- matrix(NA, J, n)
      for (j in 1:J) {
        Yj <- rfam(n, fv, sige)
        form=as.formula(paste("Yj ~ ",paste(colnames(fit$X)[-1],
                                            collapse="+")))
        mj <- slrmss(form,
                     data = data.frame(Yj,fit$X),statistic="Wald",
                     testingbeta = colnames(fit$X)[2],family=family,
                     xi=fit$xi)
        mrq[j,] <- qfam(pfam(mj$std.residuals, 0, 1),mu=0,sigma=1)
        mrq[j,] <- sort(mrq[j,])
      }
      conf <- 0.95
      infsup <- apply(mrq, 2, quantile, probs = c((1 - conf) / 2,
                                                  (1 + conf) / 2), type = 6)
      media <- colMeans(mrq)
      faixay <- range(mrq, rqobs)
      qq0 <- qqnorm(rqobs, main = "Envelope plot",xlab="Quantile N(0,1)",
                    pch = 20,col = "blue",ylim = faixay)
      eixox <- sort(qq0$x)
      lines(eixox, media)
      lines(eixox, infsup[1,])
      lines(eixox, infsup[2,]) 
  }
