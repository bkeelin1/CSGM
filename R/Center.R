#' @title Center Shape of Procrustes aligned coordinates
#'
#' @description
#' The Center dataset consists of a large .ply mesh of a mandible which represents
#' the center shape of all individuals in the study sample based on semilandmarks
#' of five cross-sectional regions across the mandibular corpus.
#'
#' @format
#'
#' A 3D mesh .ply file with the following
#'
#' \itemize{
#'   \item{vb: mesh vertices (connection points)}
#'   \item{it: triangles of the mesh}
#'   \item{primitivetype: The method of mesh construction. In this case: triangles}
#'   \item{material: a material mask. In this case none were available}
#'   \item{normals: Mesh normals or vectors drawn perpendicular to the vertices for lighting reconstruction}
#' }
#'
#' @details
#' This .ply mesh was generated using the mesh.ply function from the geomorph package.
#' Landmarks are organized sequentially by cross-section:
#'
#' @usage data(Center)
#'
#' @source Keeling et al., 2025
"Center"
