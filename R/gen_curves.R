#' @title Generate a curve matrix for landmark data
#'
#' @description
#' This function generates a 'n x 3' curve matrix containing the positions of
#' semilandmarks relative to fixed landmarks. This curve matrix can be
#' subsequently used to calculate equidistant semilandmark or landmark points.
#' This function also satisfies the "gpagen" function's argument 'curves' from
#' the *geomorph* package.
#'
#' @param ncurves Numeric. Total number of distinct curves (e.g., 10).
#' @param nsemiland Numeric. Total number of landmarks in the dataset (e.g., 590).
#' @param def_semi Numeric vector. Indices of the semi-landmarks.
#' @param def_points Numeric vector. Indices of the fixed landmarks (anchors).
#' @param closed Logical. Whether INDIVIDUAL curves are closed loops.
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
#' @author Keeling et al., 2026
#'
#' @seealso \pkg{\link{geomorph}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # You should know the exact positions of your semilandmarks in your landmark configuration before using this function
#' data(semiland_indices)
#' data(anchor_indices)
#' gen_curves(ncurves = 10,
#'            nsemiland = 590,
#'            def_semi = semiland_indices,
#'            def_points = anchor_indices,
#'            closed = FALSE)
#' }

gen_curves <- function(ncurves, nsemiland, def_semi, def_points, closed = FALSE) {

  # 1. Calculate the stride (Length of one curve)
  points_per_curve <- nsemiland / ncurves

  if(floor(points_per_curve) != points_per_curve) {
    stop("Error: Total landmarks must be perfectly divisible by ncurves.")
  }

  # 2. Assign a group id to every landmark to strictly enforce boundaries
  get_group <- function(id) {
    return(ceiling(id / points_per_curve))
  }

  valid_points <- c(def_semi, def_points)
  curve_matrix <- list()
  counter <- 1

  for(s in def_semi) {
    # Theoretical neighbors
    before <- s - 1
    after  <- s + 1
    current_group <- get_group(s)

    # 3. neighbor validation
    #    A neighbor is valid ONLY if:
    #    A) It exists in the dataset (is an anchor or slider)
    #    B) It belongs to the SAME curve group as the current point

    has_before <- (before %in% valid_points) && (get_group(before) == current_group)
    has_after  <- (after  %in% valid_points) && (get_group(after)  == current_group)

    if(has_before && has_after) {
      # Standard internal sliding
      curve_matrix[[counter]] <- c(before, s, after)
      counter <- counter + 1

    } else if (closed == TRUE) {
      # Closed Loop Logic (WITHIN GROUP ONLY)
      # If we are at the start/end of a curve group, wrap around to that SPECIFIC group's limits

      group_start <- ((current_group - 1) * points_per_curve) + 1
      group_end   <- current_group * points_per_curve

      if(!has_before && s == group_start) {
        curve_matrix[[counter]] <- c(group_end, s, after)
        counter <- counter + 1
      } else if (!has_after && s == group_end) {
        curve_matrix[[counter]] <- c(before, s, group_start)
        counter <- counter + 1
      }
    }
  }

  final_mat <- do.call(rbind, curve_matrix)
  colnames(final_mat) <- c("Before", "Slider", "After")
  return(final_mat)
}
