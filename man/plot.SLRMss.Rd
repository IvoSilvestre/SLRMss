\name{plot.SLRMss}
\alias{plot.SLRMss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.SLRMss
%%  ~~function to do ... ~~
}
\usage{
\S3method{plot}{SLRMss}(x, conf=0.95, seed=2612, H0=FALSE, xlab=NULL, ylab=NULL, main=NULL)
}  
\arguments{
\item{x}{an object of class 'SLRMss'.}
\item{conf}{the confidence level. 0.95 by default.}
\item{seed}{the seed to be used. 2612 by default.}
\item{H0}{if TRUE, choose the fitted values under null hypothesis, if FALSE, choose the fitted values under alternative hypothesis (default).}
\item{xlab}{a title for the x axis.}
\item{ylab}{a title for the y axis.} 
\item{main}{a title for the plot.}
} 
\description{
Shows an envelope plot of a SLRMss fitted model. The green points are the observations among the limits, the red points are the observations above or below these limits.
}
