#' @title Generate Extreme Principal Component Shapes
#'
#' @description
#' This function generates a list containing the mean, positive, and negative
#' shapes for each Principal Component of a Principal Component Analysis.
#'
#' @param Data_PCA A PCA object which contains the following elements:
#' x (principal component scores), sdev (eigenvalues), rotation (eigenvectors),
#' center (2D matrix of mean shape of individuals within the PCA).
#' @param k An integer value indicating the number of dimensions in the
#' landmark configuration. Inputs of 2 or 3 are accepted.
#' @param PCs An integer value indcating the number of meaningful
#' Principal Components within the PCA by which the function should
#' calculate mean shapes for.
#'
#' @returns a list object containing the following:
#'
#' \itemize{
#'    \item {negative PC shape for PC 1}
#'    \item {positive PC shape for PC 1}
#'    \item {negative PC shape for PC i}
#'    \item {positive PC shape for PC i}
#'    \item {... PCs alternate continues until all interested PCs are estimated}
#' }
#'
#' @details
#' This function takes a list object containing several important variables
#' commonly calculated during principal component analysis
#' (i.e., principal component scores, eigenvalues and eigenvectors,
#' number of interested principal components, number of dimensions)
#' and estimates the landmark configurations for the minimum and maximum observations.
#' Data_PCA requires four components supplied by the *gm.prcomp* function from the
#' **Geomorph** package:
#' \itemize{
#'    \item{x           principal component scores}
#'    \item{sdev        eigenvalues}
#'    \item{rotation    eigenvectors}
#'    \item{center      consensus (center of PCA) shape as a 2D matrix}
#' }
#'
#' @author Keeling et al., 2025
#'
#' @seealso \code{\link{gm.prcomp}}
#'
#' @importFrom geomorph arrayspecs
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' library(geomorph)
#' # Import sample dataset human_corpus
#'
#' data(human_corpus)
#'
#' Data_PCA = gm.prcomp(human_corpus)
#'
#' output = get_extreme_shapes(Data_PCA,
#'                             k = 3,
#'                             PCs = 2)
#'
#' # view the shape extremes
#'
#' lolshape(output[[[1]]])
#'
#' }
get_extreme_shapes <- function(Data_PCA, # PCA output
                               k, # dimensions
                               PCs # number of principal compnents to estimate
                               ) {

  loadings <- Data_PCA$rotation # contribution of each x,y,z point to the principal components
  mean_shape <- Data_PCA$center # mean landmark configuration
  sd_pc <- sqrt(Data_PCA$sdev) # square root of the eigenvalues
  scores <- Data_PCA$x # principal component scores for each observation

  extreme_shapes <- list()

  for (i in 1:PCs) { # Calculate standard deviations for landmark configurations for each principal component

    neg_shape <- mean_shape - abs(loadings[,i] * min(scores[,i]))
    neg_shape <- arrayspecs(t(data.frame(neg_shape)), length(neg_shape)/k, k)

    pos_shape <- mean_shape + abs(loadings[,i] * max(scores[,i]))
    pos_shape <- arrayspecs(t(data.frame(pos_shape)), length(pos_shape)/k, k)

    extreme_shapes[[paste("PC", i, "Negative")]] <- neg_shape
    extreme_shapes[[paste("PC", i, "Positive")]] <- pos_shape
  }

  return(extreme_shapes)
}
