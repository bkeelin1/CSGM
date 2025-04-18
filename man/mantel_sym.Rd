% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mantel_sym.R
\name{mantel_sym}
\alias{mantel_sym}
\title{Mantel test for shared factors influencing asymmetric and symmetric landmark variation}
\usage{
mantel_sym(x, pairs, method = "pearson", perm = 1000, ...)
}
\arguments{
\item{x}{Either an p x k x n array (landmarks, dimensions, observations) or
as an object output by the \emph{bilat.symmetry} function from the \strong{geomorph} package.
Alternatively, a list with the following components will also be accepted:
\itemize{
\item a p x k x n array called symm.shape containing the symmetric shapes
for each individual.
\item a p x k x n array called asymm.shape containing the asymmetric shapes
for each individual.
}}

\item{pairs}{a vector of the same length as the number of landmarks in
argument x and only can contain the integers 1 or 2 where the integer 1
signifies landmarks on the left side and the integer 2 for the right side.}

\item{method}{a character string indicating the method used for the
mantel correlation test. For more details, please see the \emph{mantel} function
documentation from the \strong{stats} package for more details.}

\item{perm}{an integer indicating the number of permutations desired for the
mantel test.}

\item{...}{Any of the arguments from the \emph{bilat.sym} function from the
\strong{geomorph} package can be applied here. If no arguments are passed,
the defaults for the \emph{bilat.symmetry} function in the \strong{geomorph} will be applied.}
}
\value{
a list object containing the following:

\itemize{
\item A 2 x 7 data table containing the statistics for both the left and
right sided mantel tests.
\item a list object from the \emph{mantel} function containing statistics of the
mantel test for the left side landmarks.
\item a list object from the \emph{mantel} function containing statistics of the
mantel test for the right side landmarks.
}
}
\description{
This functions conducts a permuted mantel test between the
variance-covariance matrices of paired symmetric and asymmetric landmark
configurations as described in Lazic et al., 2018. This mantel test can be
used to test if the factors influencing symmetric shape variation are also
influencing the assymmetric shape variation. Additionally, this test can be
used to test for left and right side independence.
}
\details{
This function accepts either a p x n x k array or a \emph{bilat.symmetry} output
from the \strong{geomorph} package. If an array is input, the symmetric and asymmetric shape
components will be estimated using the \emph{bilat.symmemry} functionn in the \strong{geomorph}
package. Extra arguments from the \emph{bilat.symmetry} can be passed through mantel_sym.
The variance-covariance matrices of the symmetric and asymmetric shape components
are calculated and are separated into left and right sides. A mantel test
for each side between the symmetric and asymmetric shape variance-covariance
matrices is then conducted using the \emph{mantel} function from the \strong{vegan} package.
}
\examples{

# Import geomorph and vegan packages
library(geomorph)
library(vegan)

# Import sample dataset of mandibular corpus shape.
data(Corpus_Land)

# Create a vector *pairs* to satisfy the argument for the mantel_sym function
# We will split the left and right sides
pairs <- cbind(left = c(1:120, 121:240), right = c(361:480, 481:600))

# Conduct the test

output <- mantel_sym(Corpus_Land, pairs)

# View the results

output$combined

}
\references{
Lazić MM, Carretero MA, Crnobrnja-Isailović J, Kaliontzopoulou A.
Effects of environmental disturbance on phenotypic variation: an
integrated assessment of canalization, developmental stability,
modularity, and allometry in lizard head shape. The American
Naturalist 185, 44-58 (2015).
}
\seealso{
\code{\link[=mantel(), bilat.symmetry]{mantel(), bilat.symmetry()}}
}
\author{
(s)
Brian Keeling
}
