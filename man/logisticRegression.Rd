% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{logisticRegression}
\alias{logisticRegression}
\title{Estimate logistic model via optimization}
\usage{
logisticRegression(y, X)
}
\arguments{
\item{y}{A \code{double} value of the 1 column matrix containing the response of interest.}

\item{X}{An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.}
}
\description{
This function generates the outputs listed in the final project description.  First, it calculates
the initial values for optimization obtained from the least-sqaures formula.  Next, it optimizes the values using
numerical optimization.  After this, it plots the fitted logistic curve to the actual values. Also, it outputs
Prevalence, Accuracy, Sensitivity, Specificity, False Discovery Rate, Diagnostic Odds Ratio.
}
