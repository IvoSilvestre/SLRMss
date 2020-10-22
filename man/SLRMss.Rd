\name{SLRMss}
\alias{SLRMss}
\alias{SLRMss.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Symmetric Linear Regression Models for small samples
%%  ~~function to do ... ~~
}
\description{Computes Wald, Likelihood-ratio, Score or Gradient of symmetrical linear regression models. Also computes modified Likelihood-ratio, Score and Gradient for small samples.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SLRMss(formula, family, xi, statistic, testingbeta, data=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
  \item{family}{a description of the error distribution to be used in the model. There are four supported families, Normal, t-Student, Power Exponential and Logistic ("Normal", "Student", "Powerexp" and "Logistic", 
respectively)
%%     ~~Describe \code{family} here~~
}
  \item{xi}{an extra parameter of some specified error distribution. For t-Student is a positive value and for Power Exponential is a real number between -1 and 1/3.
%%     ~~Describe \code{xi} here~~
}
  \item{statistic}{the statistic which will be used. It includes "Wald", "LR", "Score" or "Gradient".
%%     ~~Describe \code{statistic} here~~
}
  \item{testingbeta}{an vector containing the names of the variables to be testing.
%%     ~~Describe \code{testingbeta} here~~
}
  \item{data}{an optional data frame containing the variables in the model.
%%     ~~Describe \code{data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{beta.coefficients}{a matrix with the estimated position parameters under alternative hypothesis.}
\item{phi}{a numeric value with the estimated precision paramater under alternative hypothesis.}
\item{beta.coefficients.h0}{a matrix with the estimated position parameters under null hypothesis.}
\item{phi.h0}{a numeric value with the estimated precision paramater under null hypothesis.}
\item{beta.fitted}{a vector with the fitted values of the model.}
\item{null.hypothesis}{the description of the null hypothesis.}
\item{statistics}{a matrix with the selected statistics and theirs p-values.}
\item{statistic.distribution}{the name of the statistics' distribution used to test null hypothesis. It always return "Chi-Squared".}
\item{df}{the degrees of freedom of the statistics' distribution. It's the length of the testingbeta vector.}
\item{std.residuals}{the standardized residuals.}
\item{AICc}{the corrected Akaike Information Criterion for small samples.}
\item{BIC}{the Bayesian Information Criterion.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Medeiros, F. M. C and Ferrari, S. L. P. (2017 Small-sample testing inference in symmetric and log-symmetric linear regression models,
Statistica Neerlandica 2017.
%% ~put references to the literature/web site here ~
}
\author{
Ivonaldo Silvestre da Silva-Junior <ivonaldojunior@hotmail.com.br> and Francisco Moises C. de Medeiros
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
#Using the orange data frame to fit the emulsion density and testing the xanthan gum quantity.

data(orange)

SLRMss(emulsion~arabicgum+xanthangum,family="Normal",
testingbeta="xanthangum",statistic="LR",data=orange)

SLRMss(emulsion~arabicgum+xanthangum+orangeoil,family="Normal",
testingbeta="xanthangum",statistic="Score",data=orange)

fit<-SLRMss(emulsion~arabicgum+xanthangum+orangeoil,family="Student",xi=3,
testingbeta="xanthangum",statistic="LR",data=orange)
summary(fit)
plot(fit)
###############################################################
###############################################################
#Using the cheese data frame to fit the cohesiveness and test xanthan gum percentage and sodium caseinate percentage.

data(cheese)

fit<-SLRMss(exp(cohe)~fat+xangum+sodcase,family="Normal",
testingbeta=c("xangum","sodcase"),statistic="Gradient",data=cheese)
summary(fit)
plot(fit)
}