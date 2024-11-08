% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CHAD.R
\name{new_mean}
\alias{new_mean}
\title{Constructor of subclass 'meanDetector' of class 'CHAD'}
\usage{
new_mean(p, leading_constant, baseline_sd, constant_penalty, estimate_mean)
}
\arguments{
\item{p}{Dimension of the data}

\item{leading_constant}{A numeric vector of length 2. The first entry is
the leading constant of the penalty function \eqn{r(k,p,t)} for
\eqn{k\geq\sqrt{p\log n}}. The second entry is the leading constant for the
penalty function when \eqn{k\leq\sqrt{p\log n}}}

\item{baseline_sd}{Baseline standard deviation for each component of the
data}
}
\value{
An object of subclass 'meanDetector' of class 'CHAD'
#' @details Uses the test statistic of Liu et al. (2021) to test for a change
in mean online
}
\description{
Constructor of subclass 'meanDetector' of class 'CHAD'
}
\examples{
detector <- new_mean(p=100, leading_constant = c(2,2),
                     constant_penalty = FALSE)
}
\references{
\itemize{
\item Liu et al. (2021) Minimax rates in sparse, high-dimensional change
 point detection. Ann. Statist. 49(2): 1081-1112.}
}
