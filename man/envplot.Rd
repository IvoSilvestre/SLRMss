\name{envplot}
\alias{envplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Envelope Plot
%%  ~~function to do ... ~~
}       
\usage{
envplot(object, J=100, conf=0.95, seed=NULL, H0=FALSE, colors=c("red","green"), pch=16, xlab, ylab, main)
}  
\arguments{
\item{object}{an object of class 'SLRMss'.}
\item{J}{the number of Monte Carlo replications. 100 by default.}
\item{conf}{the confidence level. 0.95 by default.}
\item{seed}{a optional seed.}
\item{H0}{if TRUE, choose the fitted values under null hypothesis, if FALSE, choose the fitted values under alternative hypothesis (default).}
\item{colors}{a vector with two characters. The first one represents the color of the points above the upper limit or below the lower limit and the second one represents the color of the points under the limits.}
\item{xlab}{a title for the x axis.}
\item{ylab}{a title for the y axis.} 
\item{main}{a title for the plot.}
} 
\description{
Shows an envelope plot of a SLRMss fitted model.
}
