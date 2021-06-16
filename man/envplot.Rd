\name{envplot}
\alias{envplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{envplot
%%  ~~function to do ... ~~
\usage{
envplot(object, J=100, conf=0.95, seed=NULL, H0=FALSE, xlab=NULL, ylab=NULL, main=NULL)
}  
\arguments{
\item{object}{an object of class 'SLRMss'.}
\item{J}{the number of Monte Carlo replications. 100 by default.}
\item{conf}{the confidence level. 0.95 by default.}
\item{seed}{a optional seed.}
\item{H0}{if TRUE, choose the fitted values under null hypothesis, if FALSE, choose the fitted values under alternative hypothesis (default).}
\item{xlab}{a title for the x axis.}
\item{ylab}{a title for the y axis.} 
\item{main}{a title for the plot.}
} 
\description{
Shows an envelope plot of a SLRMss fitted model. The green points are the observations among the limits, the red points are the observations above or below these limits.
}
