#' @title Corpus_Prop: Biomechanical Properties of five cross-sectional regions along the mandibular corpus
#'
#' @description
#' The Corpus_Land dataset consists of an p x k x n (landmarks, dimensions, individuals)
#' 3D landmark array of 50 recent modern humans and 590 semilandmarks which map
#' the internal and external cortical bone surface of five cross–sections taken
#' along the mandibular corpus: (1) the left side first and second intermolar
#' corpus (abbreviated as LM1-M2), (2) the left side third and fourth interpremolar
#' corpus (abbreviated as LP3-P4), (3) the mandibular symphysis, (4) the right
#' side third and fourth interpremolar corpus (abbreviated as RP3-P4), and (5)
#' the right side first and second intermolar corpus (abbreviated as RM1-M2).
#' For each cross-section, 119 equidistant semilandmarks were collected using
#' six of the semilandmarks as anchor points: (1) external buccal alveolar shelf,
#' (2) external basal corpus, (3) external lingual alveolar shelf, (4) internal
#' buccal alveolar shelf, (5) internal basal corpus, (6) internal lingual alveolar
#' shelf. When including the anchor points and semilandmarks together, the external
#' and internal cross-sectional contours contain 60 semilandmarks each and 30 semilandmarks
#' for the buccal and lingual sides.
#'
#' @format
#' A data frame with 595 rows and 92 variables:
#'
#' \itemize{
#'   \item{ID: Character. Unique specimen identifier}
#'   \item{Collection: Character. Source collection identifier}
#'   \item{Region: Character. Cross-sectional location with 5 levels: LM1M2 (left M1-M2), LP3P4 (left P3-P4), Symphysis, RP3P4 (right P3-P4), RM1M2 (right M1-M2)}
#'   \item{TA: Total area}
#'   \item{CA: Cortical area}
#'   \item{EA: Endosteal area}
#'   \item{RCT: Relative cortical thickness}
#'   \item{MaxCT: Maximum cortical thickness}
#'   \item{ACT: average cortical thickness}
#'   \item{sdCT: standard deviation of the mean cortical thickness}
#'   \item{R1: Maximum radius along the major axis of the cross-section}
#'   \item{R2: Maxmium radius along the minor axis of the cross-section}
#'   \item{Fmin: Feret minimum diameter}
#'   \item{Fmax: Feret maximum diameter}
#'   \item{Fangle: Feret angle}
#'   \item{IP: Internal perimeter of the cross-section}
#'   \item{EP: External perimeter of the cross-section}
#'   \item{PBI: Internal buccal perimeter}
#'   \item{PLI: Internal lingual perimeter}
#'   \item{PBE: External buccal perimeter}
#'   \item{PLE: External lingual perimeter}
#'   \item{Ix: Second moment of area (x-axis)}
#'   \item{Iy: Second moment of area (y-axis)}
#'   \item{Imin: Minimum second moment of area}
#'   \item{Imax: Maxmimum second moment of area}
#'   \item{J: Polar moment of area}
#'   \item{Theta_r: Theta in radians}
#'   \item{Zx: Section modulus of area (x-axis)}
#'   \item{Zy: Sectional modulus of area (y-axis)}
#'   \item{Zmax: Maximum section modulus}
#'   \item{Zmin: Minimum section modulus}
#'   \item{Zp: Polar section modulus}
#'   \item{Point_1 - Point_59: Cortical Thicknesses from the lingual alveolar margin to the buccal alveolar margin}
#' }
#'
#' @details
#' The dataset contains measurements from 50 mandibles, each measured at 5 cross-sections
#' (50 × 5 = 250 total observations). Each cross-section has multiple measurements
#' characterizing its size, shape, and mechanical properties. The data structure
#' facilitates analysis of regional variation in mandibular properties and their
#' relationship to shape variation captured in the Corpus_Land array.
#'
#' @usage data(Corpus_Prop)
#'
#' @source Keeling et al., 2025
"Corpus_Prop"
