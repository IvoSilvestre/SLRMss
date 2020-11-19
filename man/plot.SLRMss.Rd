\name{plot.SLRMss}
\alias{plot.SLRMss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.SLRMss
%%  ~~function to do ... ~~
}
\usage{
\S3method{plot}{SLRMss}(fit, conf=0.95, seed=2612)}
\arguments{
\item{fit}{an object of class 'SLRMss'.}
\item{conf}{the confidence level. 0.95 by default.}
\item{seed}{the seed to be used. 2612 by default.}
} 
\description{
Shows an envelope plot of a slrmss fitted model. The green points are the observations among the limits, the red points are the observations above or below these limits.
}
