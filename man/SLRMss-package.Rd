\name{SLRMss-package}
\alias{SLRMss-package}
\alias{SLRMss}
\docType{package}
\title{ 
\packageTitle{SLRMss}
}
\description{
\packageDescription{SLRMss}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{SLRMss}
\packageIndices{SLRMss}
}
\author{
\packageAuthor{SLRMss}

Maintainer: \packageMaintainer{SLRMss}
}
\references{
Medeiros, F. M. C and Ferrari, S. L. P. (2017 Small-sample testing inference in symmetric and log-symmetric linear regression models,
Statistica Neerlandica 2017.
}

\keyword{ package }
\examples{
library(SLRMss)
data(orange)
fit=SLRMss(emulsion~arabicgum+xanthangum+orangeoil, family="Student",
xi=3,testingbeta="xanthangum",statistic="Wald",data=orange)
summary(fit)
}
