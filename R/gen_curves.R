#' @title Generate a curve matrix for landmark data
#'
#' @description
#' This function generates a 'n x 3' curve matrix containing the positions of
#' semilandmarks relative to fixed landmarks. This curve matrix can be
#' subsequently used to calculate equidistant semilandmark or landmark points.
#' This function also satisfies the "gpagen" function's argument 'curves' from
#' the *geomorph* package.
#'
#' @param ncurves A numeric value indicating the number of semilandmarks in your
#' landmark configuration.
#' @param nfixland A numeric value indicating the number of fixed landmarks that
#' should not be slid
#' @param nsemiland A numeric value indicating the total number of semilandmarks,
#' including semilandmark anchor points that form part of a curve. This value does not
#' include landmarks in your configuration that are unrelated to a
#' semilandmark curve.
#' @param def_points An optional vector containing numerical string representing the
#' landmark numbers associated with each fixed landmark point in your landmark
#' configuration.
#' @param def_semi An optional vector containing a numerical string representing the
#' landmark numbers associated with each semilandmark point in your landmark
#' configuration (not anchor points in a semilandmark configuration).
#' This argument is paired with def_points. Only applicable
#' when the landmark configurations contain fixed points independent of
#' the semilandmark curves.
#' @param closed A logical TRUE/FALSE statement indicating whether the curve
#' is closed in a loop.
#'
#' @details This function is primarily for landmark and semilandmark
#' configurations. A n x 3 landmark curve matrix is generated as a result which
#' can be used as an argument for functions to manage landmark spacing. This
#' function is only recommended for situations when the landmark configuration
#' contains (1) only semi-landmarks whose spacing needs further management,
#' (2) only landmarks that need to be equidistantly spaced, or (3) the landmark
#' configuration contains both semilandmarks and landmarks.
#'
#' @returns This function returns a n x 3 landmark curve matrix where n is
#' dependent on the number of landmarks and semilandmarks in the
#' landmark configuration.
#'
#' @author Keeling et al., 2025
#'
#' @seealso \pkg{\link{geomorph}}
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # You should know the exact positions of your semilandmarks in your landmark configuration before using this function
#'
#' # **Simple Example 1**
#'
#' # In this example, we only have 1 semilandmark curve totaling 11 semilandmarks including three fixed landmark anchor points.
#'
#' Landmark_matrix <- gen_curves(ncurves = 1, # There is only 1 semilandmark curve
#'                               nfixland = 3, # we have 3 anchor points
#'                               nsemiland = 11, # we have 11 semilandmarks including the anchor points
#'                               def_points = NULL,
#'                               def_semi = NULL,
#'                               closed = 0 # This is an open semilandmark configuration.
#'                               )
#'
#'
#'  #**Simple Example 2**
#'
#' # In this example, we only have 10 type 1 landmarks followed by 20 semilandmarks with two anchor points (first and last)
#'
#' Landmark_matrix <- gen_curves(ncurves = 1, # There is only 1 semilandmark curve
#'                               nfixland = 12, # we have 2 anchor points plus 10 fixed landmarks
#'                               nsemiland = 20, # we have 20 semilandmarks including the anchor points
#'                               def_points = c(1:10, # type 1 landmarks
#'                                              11, 30, # Keep in mind that these are the anchor points for the open curve
#'                                              ),
#'                               def_semi = c(11:30), # You must include anchor points as well
#'                               closed = 0 # There are no closed semilandmarks
#'                               )
#'
#'
#' # **Complex Example**
#'
#' # In this example, let's assume we have a landmark configuration mapping a human cranium:
#'
#' \itemize{
#'   \item There are 32 type 1 landmarks.
#'   \item There are 1 *closed* semilandmark curve of 30 points mapping the right temporalis muscle attachment site.
#'   \item There is 1 *open* semilandmark curve of 11 points mapping the right superior temporal line.
#'   \item Specficially,
#'   \itemize{
#'     \item Landmarks 1-16 are type 1 landmarks.
#'     \item Landmarks 17 - 28 are open semilandmarks.
#'     \item Landmarks 29 - 39 are type 1 landmarks.
#'     \item Landmarks 40 - 70 are closed semilandmarks.
#'     \item Landmarks 71 - 75 are type 1 landmarks.
#'   }
#' }
#'
#'
#' Landmark_matrix <- gen_curves(ncurves = 2,
#'                               nfixland = 32,
#'                               nsemiland = 41,
#'                               def_points = c(1:16, # type 1 landmarks
#'                                              17, 22 ,28 # Keep in mind that these are the anchor points for the open curve
#'                                              29:39, # type 1 landmarks
#'                                              # Since the curve is close we don't include anchor points for the 2nd curve
#'                                              71:75 # type 1 landmarks
#'                                              ),
#'                               def_semi = c(17:28, 40:70) # You must include anchor points as well,
#'                               closed = 1 # There is only 1 closed semilandmark curve in this example. Default is 0.
#'                               )
#'
#'}
gen_curves <- function(ncurves,
                       nfixland,
                       nsemiland,
                       def_points = NULL,
                       def_semi = NULL,
                       closed = FALSE) {

  if (is.null(def_semi) || is.null(def_points)) { # if the curves are equally divided by the number of landmarks
    land.curve = nsemiland/ncurves
    curve.list = vector("list", ncurves)

    for(x in 1:ncurves) {

      curve.list[x] <- list(matrix(NA, land.curve-nfixland, 3))

      if (x <= 1) {
        pattern <- c(1, 2, 3)
        curve <- matrix(NA, land.curve-nfixland, 3)
        for (i in 1:(land.curve - nfixland)) {
          curve[i, ] <- pattern
          pattern <- pattern + 1
          if (pattern[3] > land.curve) {
            pattern <- pattern - 2
          }
        }
        curve.list[[x]] <- curve

      } else {

        pattern2 <- c((land.curve * (x-1)) + 1, (land.curve * (x-1)) + 2, (land.curve * (x-1)) + 3)
        curve2 <- matrix(NA, land.curve-nfixland, 3)

        for (i in 1:(land.curve - nfixland)) {
          curve2[i, ] <- pattern2
          pattern2 <- pattern2 + 1

          if (pattern2[3] > land.curve * x) {
            pattern2 <- pattern2 - 2}
        }
        curve.list[[x]] <- curve2
      }

      curve_matrix <- do.call(rbind, curve.list)

      if (closed == TRUE) {
        curve_matrix <- rbind(curve_matrix, c(pattern2[1]+land.curve, pattern2[1], pattern2[2]))
      }
    }
  } else {

    curve_start_index = def_semi[1] - 1

    n <- length(def_semi) - sum(c(any(seq_along(def_semi)) + 1 & any(seq_along(def_semi)) - 1
                                  %in% any(seq_along(def_points))))

    curve_matrix <- matrix(NA, nrow = n, ncol = 3)

    for(i in 1:n) {

      if (i == 1) {
        curve_matrix[1,] <- c(curve_start_index, def_semi[1], def_semi[1 + 1])

      } else if (i == n + 1) {
        curve_matrix[i,] <- curve_matrix[i, ] <- c(def_semi[i-1], def_semi[i-1] + 1, curve_start_index)

      } else if (closed == TRUE & def_semi[i] > def_semi[i-1] + 2) {
        curve_matrix[i, ] <- c(def_semi[i-1], def_semi[i-1] + 1, curve_start_index)
        curve_start_index <- def_semi[i] - 1

      } else if (def_semi[i] + 2 == any(seq_along(def_points)) & i != length(n)) {
        curve_matrix[i,] <- c(def_semi[i] + 2, def_semi[i] + 3, def_semi[i]+ 4)
        curve_start_index <- def_semi[i]+ 2

      } else if (def_semi[i] - 1 == any(seq_along(def_points)) & i != length(n)) {
        curve_matrix[i, ] <- c(def_semi[i] - 1, def_semi[i], def_semi[i] + 1)

      } else if (def_semi[i] + 1 == any(seq_along(def_points)) & i != length(n)) {
        curve_matrix[i,] <- c(def_semi[1] - 1, def_semi[i], def_semi[i]+1)

      } else if (def_semi[i] + 2 != any(seq_along(def_points))) {
        curve_matrix[i,] <- c(def_semi[i] - 1, def_semi[i], def_semi[i]+1)

      }
    }
    curve_matrix <- na.omit(curve_matrix)
  }
  return(curve_matrix)
}
