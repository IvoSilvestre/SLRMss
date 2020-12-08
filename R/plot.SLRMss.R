plot.SLRMss <-
function (fit, conf = 0.95, seed = 2612, H0=F) 
{
    family = fit$family
    if (family == "Normal") {
        qfam = function(x, mu, sigma) {
            qnorm(x, mu, sigma)
        }
        pfam = function(x, mu, sigma) {
            pnorm(x, mu, sigma)
        }
        rfam = function(x, mu, sigma) {
            rnorm(x, mu, sigma)
        }
    }
    else {
        if (family == "Student") {
            qfam = function(x, mu, sigma) {
                sigma * qt(x, df = fit$xi) + mu
            }
            pfam = function(x, mu, sigma) {
                pt((x - mu)/sigma, df = fit$xi)
            }
            rfam = function(x, mu, sigma) {
                sigma * rt(x, df = fit$xi) + mu
            }
        }
        else {
            if (family == "Logistic") {
                qfam = function(x, mu, sigma) {
                  qlogis(x, mu, sigma)
                }
                pfam = function(x, mu, sigma) {
                  plogis(x, mu, sigma)
                }
                rfam = function(x, mu, sigma) {
                  rlogis(x, mu, sigma)
                }
            }
            else {
                qfam = function(x, mu, sigma) {
                  qnormp(x, mu, sigma, p = 2/(fit$xi + 1))
                }
                pfam = function(x, mu, sigma) {
                  pnormp(x, mu, sigma, p = 2/(fit$xi + 1))
                }
                rfam = function(x, mu, sigma) {
                  rnormp(x, mu, sigma, p = 2/(fit$xi + 1))
                }
            }
        }
    }
    set.seed(seed)
    y = fit$y
    n = length(y)
    J = 100
    fv = fitted(fit,H0=H0)
    sige = coef(fit,H0=H0)$phi
    rqobs <- residuals(fit,H0=H0,std=T)
    mrq <- matrix(NA, J, n)
    for (j in 1:J) {
        Yj <- rfam(n, fv, sige)
        form = as.formula(paste("Yj ~ ", paste(colnames(fit$X)[-1], 
            collapse = "+")))
        mj <- SLRMss(form, data = data.frame(Yj, fit$X), statistic = "Wald", 
            testingbeta = fit$testingbeta, family = family, 
            xi = fit$xi)
        mrq[j, ] <- residuals(mj,H0=H0,std=T)
        mrq[j, ] <- sort(mrq[j, ])
    }
    conf <- 0.95
    infsup <- apply(mrq, 2, quantile, probs = c((1 - conf)/2, 
        (1 + conf)/2), type = 6)
    media <- colMeans(mrq)
    faixay <- range(mrq, rqobs)
    qq0 <- qqnorm(rqobs, main = paste("Envelope plot -",100*conf,"% confidence"), xlab = "Quantile N(0,1)", 
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
}
