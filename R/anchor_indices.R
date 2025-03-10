#' @title anchor_indices: A string of landmark indices which are anchor points
#'
#' @description
#' The anchor_indices dataset is a numeric vector identifying the exact landmarks
#' which are anchor points within the Corpus_Land dataset. This is useful for
#' using the gen_curves function.
#'
#' @format
#' A numeric string of 30 anchor points, 6 for each of the following landmarks
#' within Corpus_Land:
#'
#' \itemize{
#'   \item{LM1M2 (1-118)        Left first-second molar junction}
#'   \item{LP3P4 (119-236)      Left third-fourth premolar junction}
#'   \item{Symphysis (237-354)  Mandibular symphysis}
#'   \item{RP3P4 (355-472)      Right third-fourth premolar junction}
#'   \item{RM1M2 (473-590)      Right first-second molar junction}
#' }
#'
#' @usage data(anchor_indices)
#'
#' @source Keeling et al., 2025
"anchor_indices"
