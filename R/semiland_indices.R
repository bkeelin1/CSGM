#' @title semiland_indices: A string of landmark indices which are semilandmark points
#'
#' @description
#' The semiland_indices dataset is a numeric vector identifying the exact landmarks
#' which are semilandmark points within the Corpus_Land dataset. This is useful for
#' using the gen_curves function.
#'
#' @format
#' A numeric string of 560 semilandmark points, 114 for each of the following landmarks
#' within Corpus_Land:
#'
#' \itemize{
#'   \item{LM1M2 (1-118):   Left first-second molar junction}
#'   \item{LP3P4 (119-236):   Left third-fourth premolar junction}
#'   \item{Symphysis (237-354):   Mandibular symphysis}
#'   \item{RP3P4 (355-472):   Right third-fourth premolar junction}
#'   \item{RM1M2 (473-590):   Right first-second molar junction}
#' }
#'
#' @usage data(semiland_indices)
#'
#' @source Keeling et al., 2025
"semiland_indices"
