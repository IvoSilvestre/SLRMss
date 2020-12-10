\name{orange}
\alias{orange}
\docType{data}
\title{Orange Data Set
%%   ~~ data name/kind ... ~~
}
\description{ The Mirhosseini's e Tan's data set gives the measurements of the physicochemical properties of emulsion in orange drinks.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("orange")}
\format{
  A data frame with 20 observations on the following 6 variables.
  \describe{
    \item{\code{emulsion}}{the emulsion density in \eqn{g/cm^3}}
    \item{\code{arabicgum}}{the arabic gum quantity in \eqn{g}}
    \item{\code{xanthangum}}{the xanthan gum quantity in \eqn{g}}
    \item{\code{orangeoil}}{the orange oil quantity in \eqn{g}}
  }
}

\source{
Mirhosseini, H., Tan, C.P. (2010). Discrimination of orange beverage emulsions with different formulations
using multivariate analysis. Journal of the Science of Food and Agriculture.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(orange)
plot(orange)
fit=SLRMss(emulsion~arabicgum+xanthangum+orangeoil,family="Student",
xi=3,testingbeta="xanthangum",statistic="LR",data=orange)
summary(fit)
plot(fit)
}
\keyword{datasets}
