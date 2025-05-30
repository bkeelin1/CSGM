% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{trim_scores}
\alias{trim_scores}
\title{Trim VIP scores}
\usage{
trim_scores(scores, lv_method)
}
\arguments{
\item{scores}{data table consisting of two variables: Response, Predictor, and Score.
This is an object generated by the \code{Bio.VIP} function}

\item{lv_method}{Character string indicating the method to trim the predictor
variables by its VIP score. Options are "trim" (all predictor variables with
a score greater than 1), "median" (all predictor variables which are greater
than the median VIP score), and "mean" (all predictor variables which are
greater than the mean VIP score).}
}
\value{
a vector containing the names of the VIP selected predictor variables.
}
\description{
Trim Variable Importance in Projection (VIP) scores of predictor
variables by the mean, median, or scores greater than one.
}
\details{
This function estimates important variable predictors using several
standard summary statistics such as the mean and median score as well as the
rule of thumb as being greater than one (the latent variable block variable
contribution exceeds the mean).
}
\author{
Brian Anthony Keeling
}
\keyword{internal}
