#' @title A Principal Component Analysis Statistics Compiler
#'
#' @description
#' This function takes the results of a principal component analysis object
#' and outputs standard statistics of the principal component analysis in a
#' table format. This function generates a table "PCA_var" which can be used
#' with the function 'Lolligen' which provides information on the variance
#' explained by each principal component as well as an analysis of the
#' meaningful number of Principal Components based on the data.
#'
#' @param Data_PCA A PCA object which contains the following elements:
#' x (principal component scores), sdev (eigenvalues), rotation (eigenvectors),
#' center (mean shape of individuals within the PCA).
#'
#' @param PCs An optional logical value indicating whether to calculate the number
#' of meaningful PCs. As this involves factor analysis, this may take some time to
#' calculate. The default is to perform Principal Component factor analysis.
#'
#' @details
#' This function takes a list object containing several important variables
#' commonly calculated during principal component analysis
#' (i.e., principal component scores, eigenvalues and eigenvectors,
#' number of interested principal components, number of dimensions)
#' and generates several objects: a table containing the percentages of
#' explained variance for each principal components, a scree plot of the
#' explained variance by each principal component, and the number of meaningful
#' principal components. The meaningful number of principal components uses
#' factor analysis via the function **fa.parallel** from the *psych* package.
#' This function examines the percentages of variance across each principal
#' component in comparison to a random data matrix of the same size to assess
#' the most important number of principal components that can explain the
#' most significant percentages of variance. If the factor analysis is unable
#' to detect the most important number of principal components, by default the
#' number of principal components that can explain at least 80% of the total
#' cumulative variance will be selected.
#'
#' @returns a list object containing the following:
#'
#' \itemize{
#'    \item *pca_var*:   data table of dimensions pc x 3 where pc is the number of principal
#'    components obtained from the prinicipal component analysis.
#'    \item *scree*:   principal component scree plot generated from the data table.
#'    \item *PCs*:   number of meaningful principal components.
#' }
#'
#' @author Keeling et al., 2025
#'
#' @export
#'
#' @importFrom geomorph two.d.array
#' @importFrom FactoMineR PCA
#' @importFrom factoextra get_pca
#' @importFrom psych fa.parallel
#'
#' @examples
#'
#' \dontrun{
#'
#' # Import sample dataset human_corpus
#' data(Corpus_Land)
#'
#' # Conduct a Principal Component Analysis on the 3D array.
#' Data_PCA = gm.prcomp(Corpus_Land)
#'
#' # execute the function
#' output = pca_variances(Data_PCA) # object from the PCA
#'
#' # View the scree plot
#' output$Scree
#'
#' # View the data table
#' output$PCA_var
#'
#' # View the number of meaningful PCs
#' output$PCs
#'
#' }

PCA_variances <- function(Data_PCA, PCs = TRUE) {

  if(is.null(Data_PCA$A)) {
    scores <- Data_PCA$x
    loadings <- Data_PCA$rotation
    center <- Data_PCA$center
    scale <- Data_PCA$scale

    # Reconstruct the original data
    reconstructed_data <- scores %*% t(loadings)

    # Reverse scaling if it was applied
    if (!is.null(scale)) {
      reconstructed_data <- sweep(reconstructed_data, 2, scale, '*')
    }

    # Add the mean to reverse centering
    if (!is.null(center)) {
      reconstructed_data <- sweep(reconstructed_data, 2, center, '+')
    }

    Data = reconstructed_data

  } else {

    if (length(dim(Data_PCA$A)) == 3) {
      Data = data.frame(two.d.array(Data_PCA$A))
    } else {
      Data = Data_PCA$A
    }
  }

  PCA_variances <- list()

  PCA_variances$PCA_Data <- Data_PCA

  variances <- Data_PCA$sdev^2
  total_variance <- sum(variances)
  perc_var <- (variances / total_variance) * 100
  cum_var <- cumsum(perc_var)

  # Create a data frame
  PCA_var <- data.frame(PC = seq_along(perc_var),
                        var = perc_var, cum_var = cum_var)

  PCA_variances$PCA_var <- PCA_var

  PCA_variances$Scree <- plot(PCA_var$PC, PCA_var$var, type = "b", xlab = "Principal Component",
                              ylab = "Percentage of Variance Explained", main = "Scree Plot")
  if (length(dev.list()) > 0) {
    graphics.off()
  }

  if(isTRUE(PCs)) {
    cumvar = min(which(PCA_var$cum_var >= 80))
    parallel_which_PC = suppressMessages(suppressWarnings(psych::fa.parallel(Data, fa = "pc", n.iter = 1, plot = FALSE)))

    closeAllConnections()

    if (parallel_which_PC$ncomp == 0 | is.na(parallel_which_PC$ncomp)) {
      PCs = cumvar
    } else {
      PCs = min(as.numeric(parallel_which_PC$ncomp), as.numeric(cumvar))
    }
    PCA_variances$PCs <- PCs
  } else {
    PCA_variances$PCs <- NULL
  }

  return(PCA_variances)
}
