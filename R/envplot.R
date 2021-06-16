envplot <-
function (x, J=100, conf = 0.95, seed = NULL, H0=FALSE,xlab=NULL,ylab=NULL,main=NULL) 
{   
    fit = x
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
    if(is.null(xlab)) xlab="Quantile N(0,1)"
    if(is.null(ylab)) ylab="Sample Quantiles"
    if(is.null(main)) main=paste0("Envelope plot - ",100*conf,"% confidence")
    ylim=c(min(faixay,min(rqobs)),max(faixay,max(rqobs)))
    qq0 <- qqnorm(rqobs, main = main, xlab = xlab,ylab=ylab, 
        pch = 1, col = "white", ylim = faixay)
    eixox <- sort(qq0$x)
    for (i in 1:length(qq0$x)) {
        if (sort(qq0$y)[i] < infsup[1, i] | sort(qq0$y)[i] > 
            infsup[2, i]) {
            points(eixox[i], sort(qq0$y)[i], col = "red", 
                pch = 20)
        }
        else {
            points(eixox[i], sort(qq0$y)[i], col = "green", 
                pch = 20)
        }
    }
    lines(eixox, media)
    lines(eixox, infsup[1, ])
    lines(eixox, infsup[2, ])
    out=NULL
    out$x=qq0$x
    out$y=qq0$y
    invisible(out)
}
