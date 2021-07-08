\name{envplot}
\alias{envplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Envelope Plot
%%  ~~function to do ... ~~
}       
\usage{
envplot(object, J=100, conf=0.95, seed=NULL, H0=FALSE, colors=c("red","green"), pch=16, lty=2, xlab, ylab, main)
}  
\arguments{
\item{object}{an object of class 'SLRMss'.}
\item{J}{the number of Monte Carlo replications. 100 by default.}
\item{conf}{the confidence level. 0.95 by default.}
\item{seed}{a optional seed.}
\item{H0}{if TRUE, choose the fitted values under null hypothesis, if FALSE, choose the fitted values under alternative hypothesis (default).}
\item{colors}{a vector with one or two characters. If it has one character, that represents the color of the plotted points. If it has two characters, the first one represents the color of the points out of the limits and the second one represents the color of the points under the limits. Red and green by default.}
\item{pch}{a vector with one or two numerics. If it has one numeric, that represents the plot pch. If it has two numerics, the first one represents the pch of the points out of the limits and the second one represents the pch of the points under the limits. 16 by default.}
\item{lty}{a vector with one or two numerics. If it has one numeric, that represents the lty of all lines. If it has two numerics, the first one represents the lty of the middle line and the second one represents the lty of the limits line. 2 by default.}
\item{xlab}{a title for the x axis.}
\item{ylab}{a title for the y axis.} 
\item{main}{a title for the plot.}
} 
\description{
Shows an envelope plot of a SLRMss fitted model.
}
