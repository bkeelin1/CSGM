#' @title Mantel test for shared factors influencing asymmetric and symmetric landmark variation
#'
#' @description
#' This functions conducts a permuted mantel test between the
#' variance-covariance matrices of paired symmetric and asymmetric landmark
#' configurations as described in Lazic et al., 2018. This mantel test can be
#' used to test if the factors influencing symmetric shape variation are also
#' influencing the assymmetric shape variation. Additionally, this test can be
#' used to test for left and right side independence.
#'
#' @param x Either an p x k x n array (landmarks, dimensions, observations) or
#' as an object output by the *bilat.symmetry* function from the **geomorph** package.
#' Alternatively, a list with the following components will also be accepted:
#' \itemize{
#'    \item a p x k x n array called symm.shape containing the symmetric shapes
#'    for each individual.
#'    \item a p x k x n array called asymm.shape containing the asymmetric shapes
#'    for each individual.
#' }
#'
#' @param pairs a vector of the same length as the number of landmarks in
#' argument x and only can contain the integers 1 or 2 where the integer 1
#' signifies landmarks on the left side and the integer 2 for the right side.
#'
#' @param method a character string indicating the method used for the
#' mantel correlation test. For more details, please see the *mantel* function
#' documentation from the **stats** package for more details.
#'
#' @param perm an integer indicating the number of permutations desired for the
#' mantel test.
#'
#' @param ... Any of the arguments from the *bilat.sym* function from the
#' **geomorph** package can be applied here. If no arguments are passed,
#' the defaults for the *bilat.symmetry* function in the **geomorph** will be applied.
#'
#' @details
#' This function accepts either a p x n x k array or a *bilat.symmetry* output
#' from the **geomorph** package. If an array is input, the symmetric and asymmetric shape
#' components will be estimated using the *bilat.symmemry* functionn in the **geomorph**
#' package. Extra arguments from the *bilat.symmetry* can be passed through mantel_sym.
#' The variance-covariance matrices of the symmetric and asymmetric shape components
#' are calculated and are separated into left and right sides. A mantel test
#' for each side between the symmetric and asymmetric shape variance-covariance
#' matrices is then conducted using the *mantel* function from the **vegan** package.
#'
#' @returns a list object containing the following:
#'
#' \itemize{
#'    \item A 2 x 7 data table containing the statistics for both the left and
#'    right sided mantel tests.
#'    \item a list object from the *mantel* function containing statistics of the
#'    mantel test for the left side landmarks.
#'    \item a list object from the *mantel* function containing statistics of the
#'    mantel test for the right side landmarks.
#' }
#'
#' @references
#' Lazić MM, Carretero MA, Crnobrnja-Isailović J, Kaliontzopoulou A.
#' Effects of environmental disturbance on phenotypic variation: an
#' integrated assessment of canalization, developmental stability,
#' modularity, and allometry in lizard head shape. The American
#' Naturalist 185, 44-58 (2015).
#'
#' @seealso \code{mantel} \code{bilat.symmetry}
#'
#' @author Keeling et al., 2025
#'
#' @importFrom geomorph bilat.symmetry
#' @importFrom vegan mantel
#' @importFrom stats cov cov2cor
#'
#' @examples
#'
#' # Import geomorph and vegan packages
#' library(geomorph)
#' library(vegan)
#'
#' # Import sample dataset of mandibular corpus shape.
#' data(Corpus_Land)
#'
#' # Create a vector *pairs* to satisfy the argument for the mantel_sym function
#' # We will split the left and right sides
#' pairs <- cbind(left = c(1:120, 121:240), right = c(361:480, 481:600))
#'
#' # Conduct the test
#'
#' output <- mantel_sym(Corpus_Land, pairs)
#'
#' # View the results
#'
#' output$combined
#'


mantel_sym <- function(x,
                       pairs,
                       method = "pearson",
                       perm = 1000,
                       ...
                       ){

  Mantel_results = list()

  if(is.array(x) || is.matrix(x)) {
    dots <- list(...)
    if(!"ind" %in% names(dots)) {
      dots$ind <- seq_len(dim(x)[3])
    }
    x = do.call(geomorph::bilat.symmetry,
                c(list(A = x,
                       land.pairs = pairs,
                       print.progress = FALSE),
                  dots))
  }

  left_side_sym = x$symm.shape[c(pairs[,1]), , ]
  left_cov_sym = cov(two.d.array(left_side_sym))
  left_cor_sym = cov2cor(left_cov_sym)
  left_var_sym = apply(left_cor_sym, 2, var)
  left_vcv_sym = left_var_sym * left_cor_sym

  left_side_asym = x$asymm.shape[c(pairs[,1]), , ]
  left_cov_asym = cov(two.d.array(left_side_asym))
  left_cor_asym = cov2cor(left_cov_asym)
  left_var_asym = apply(left_cor_asym, 2, var)
  left_vcv_asym = left_var_asym * left_cor_asym

  # Compute the Mantel Test of the left VCV matrices of the symmetric and asymmetric shape components
  mantel_test_left <- vegan::mantel(left_cov_sym, left_cov_asym, method = method, permutations = perm, na.rm = TRUE)
  mantel_test_left

  # Second, calculate the VCV matrices for the right side symmetric and asymmetric components
  right_side_sym = x$symm.shape[c(pairs[,2]), , ]
  right_cov_sym = cov(two.d.array(right_side_sym))
  right_cor_sym = cov2cor(right_cov_sym)
  right_var_sym = apply(right_cor_sym, 2, var)
  right_vcv_sym = right_var_sym * right_cor_sym

  right_side_asym = x$asymm.shape[c(pairs[,2]), , ]
  right_cov_asym = cov(two.d.array(right_side_asym))
  right_cor_asym = cov2cor(right_cov_asym)
  right_var_asym = apply(right_cor_asym, 2, var)
  right_vcv_asym = right_var_asym * right_cor_asym

  # Compute the Mantel Test of the right VCV matrices of the symmetric and asymmetric shape components
  mantel_test_right <- vegan::mantel(right_cov_sym, right_cov_asym, method = method, permutations = perm, na.rm = TRUE)
  mantel_test_right

  Mantel_results$combined <- rbind(mantel_test_left, mantel_test_right)
  Mantel_results$left <- mantel_test_left
  Mantel_results$right <- mantel_test_right

  return(Mantel_results)
} # end of function
