\name{cheese}
\alias{cheese}
\docType{data}
\title{Cheese Data set
%%   ~~ data name/kind ... ~~
}
\description{ The Nateghi's cheese data set gives measurements of the cheddar cheese textures.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("cheese")}
\format{
  A data frame with 16 observations on the following 4 variables.
  \describe{
    \item{\code{cohe}}{the cohesiveness}
    \item{\code{fat}}{the fat percentage}
    \item{\code{xangum}}{the xanthan gum percentage}
    \item{\code{sodcase}}{the sodium caseinate percentage}
  }
}

\source{
Nateghi, L., Roohinejad, S., Totosaus, A., Mirhosseini, H., Shuhaimi, M., Meimandipour, A., Omidizadeh, A., Manap, M.Y.A (2012). Optimization of textural properties and formulation of reduced fat Cheddar cheeses containing fat replacers. Journal of Food, Agriculture & Environment.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(cheese)
plot(cheese)
plot(exp(cheese$cohe),cheese$fat)
fit=SLRMss(exp(cohe)~fat+xangum,family="Normal",testingbeta="xangum",
statistic="Gradient",data=cheese)
summary(fit)
}
\keyword{datasets}
