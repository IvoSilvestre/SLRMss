\name{plot.SLRMss}
\alias{plot.SLRMss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.SLRMss
%%  ~~function to do ... ~~
}
\usage{
\S3method{plot}{SLRMss}(x,H0=FALSE,
                        xlab=c("Fitted Values","Index",
                                 "Theoretical Quantiles","Standardized Residuals"),
                        ylab=c("Standardized Residuals", "Standardized Residuals",
                               "Standardized Residuals","Density"),
                        main=c("Against Fitted Values","Against Index",
                               "Normal Q-Q Plot","Density Estimate"))
}
\arguments{
\item{x}{an object of class 'SLRMss'.}
\item{H0}{if TRUE, plot the graphics under null hypothesis, if FALSE, plot the graphics under alternative hypothesis (default).}
\item{xlab}{a vector containing the four x-axis titles.}
\item{ylab}{a vector containing the four y-axis titles.}
\item{main}{a vector containing the four main plot titles.}
}
\description{
Shows four plot for residual analysis. The first plot shows the standardized residuals against the fitted values. The second one shows the standardized residuals against by their index. The third one presents QQ-normal plot of them and the last one shows their density estimate.
}
