% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_coef}
\alias{get_coef}
\title{Extract Coefficients from Multivariate Adaptive Regression Splines}
\usage{
get_coef(coefficients, Predictor_names)
}
\arguments{
\item{coefficients}{A named vector of coefficients from a MARS model}

\item{Predictor_names}{A character vector containing the names of predictor
variables used in the model}
}
\value{
Returns a matrix with:
\itemize{
\item Rows corresponding to predictor variables
\item A single column containing the aggregated coefficients
\item Row names matching the original predictor names
}
}
\description{
This function extracts and processes coefficients from Multivariate Adaptive
Regression Splines (MARS) models, specifically handling hinge functions and
aggregating coefficients for each predictor. It is particularly useful for
interpreting MARS models and understanding variable importance.
}
\details{
Extract and Process Coefficients from MARS Models

The function processes MARS model coefficients in several steps:
\itemize{
\item Identifies and extracts base variable names from hinge functions
\item Aggregates coefficients for each predictor across multiple terms
\item Handles both linear and non-linear (hinge) terms
\item Creates a standardized output format for easy interpretation
}

The function is particularly useful for models created using the earth package
or similar MARS implementations, where coefficients may be spread across
multiple hinge functions for the same predictor.
}
\examples{
\dontrun{
# Create sample coefficients
coef_vector <- c(
  "(Intercept)" = 1.5,
  "h(x1-2)" = 0.5,
  "h(3-x1)" = -0.3,
  "x2" = 0.7
)

# Define predictor names
pred_names <- c("x1", "x2")

# Extract coefficients
result <- get_coef(coef_vector, pred_names)
print(result)
}

}
\seealso{
\code{\link[earth]{earth}}
}
\author{
Brian Anthony Keeling
}
\keyword{internal}
