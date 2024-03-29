envplot <-
function (object, J=100, conf = 0.95, seed = NULL, H0=FALSE,colors=c("red","green"),pch=16,lty=2,xlab,ylab,main) 
{   
    if(class(object)!="SLRMss"){
    stop("Object class must be 'SLRMss'.")
    }
    fit=object
    family = fit$family
    if (family == "Normal") {
        rfam = function(x, mu, sigma) {
            rnorm(x, mu, sigma)
        }
    }
    else {
        if (family == "Student") {
            rfam = function(x, mu, sigma) {
                sigma * rt(x, df = fit$xi) + mu
            }
        }
        else {
            if (family == "Logistic") {
                rfam = function(x, mu, sigma) {
                  rlogis(x, mu, sigma)
                }
            }
            else {
                rfam = function(x, mu, sigma) {
                  rnormp(x, mu, sigma, p = 2/(fit$xi + 1))
                }
            }
        }
    }
    set.seed(seed)
    y = fit$y
    n = length(y)
    fv = fitted(fit,H0=H0)
    sige = coef(fit,H0=H0)$phi
    rqobs <- residuals(fit,H0=H0,std=TRUE)
    mrq <- matrix(NA, J, n)
    for (j in 1:J) {
        Yj <- rfam(n, fv, sige)
        form = as.formula(paste("Yj ~ ", paste(colnames(fit$X)[-1], 
            collapse = "+")))
        mj <- SLRMss(form, data = data.frame(Yj, fit$X), statistic = "Wald", 
            testingbeta = fit$testingbeta, family = family, 
            xi = fit$xi)
        mrq[j, ] <- residuals(mj,H0=H0,std=TRUE)
        mrq[j, ] <- sort(mrq[j, ])
    }
    infsup <- apply(mrq, 2, quantile, probs = c((1 - conf)/2, 
        (1 + conf)/2), type = 6)
    media <- colMeans(mrq)
    faixay <- range(infsup)
    if(missingArg(xlab)) xlab="Quantile N(0,1)"
    if(missingArg(ylab)) ylab="Standardized Residuals"
    if(missingArg(main)) main=paste0("Envelope plot - ",100*conf,"% confidence")
    ylim=c(min(faixay,min(rqobs)),max(faixay,max(rqobs)))
    if(length(colors)>=2){
        if(length(colors)>2) warning("Only the first two colors entries were used.")
        colors=ifelse(sort(rqobs)<infsup[1,] | sort(rqobs) > infsup[2,],colors[1],colors[2])
    }
    if(length(pch)>=2){
        if(length(pch)>2) warning("Only the first two pch entries were used.")
        pch=ifelse(sort(rqobs)<infsup[1,] | sort(rqobs) > infsup[2,],pch[1],pch[2])
    }
    qq0 <- qqnorm(sort(rqobs), main = main, xlab = xlab,ylab=ylab, 
        col=colors,pch=pch, ylim = faixay)
    eixox <- sort(qq0$x)
    if(length(lty)>2)warning("Only the first two lty entries were used")
    if(length(lty)==1) lty[2]=lty
    lines(eixox, media,lty=lty[1])
    lines(eixox, infsup[1, ],lty=lty[2])
    lines(eixox, infsup[2, ],lty=lty[2])
    out=NULL
    out$x=qq0$x
    out$y=qq0$y
    invisible(out)
}
