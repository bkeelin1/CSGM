#' @title Extract Package Vignette to conduct sample analysis
#'
#' @description This function finds the package's vignette and extracts it as
#' an .Rmd RMarkdown file to be used in the script panel. The file will appear in
#' the project directory.
#'
#' @param vignette_name the name in quotations of the vignette you want to extract
#' @param dest_dir the directory destination you want the file to be stored
#'
#' @return a .Rmd RMarkdown file of the package vignette in your project directory
#'
#' @export

get_vignette <- function(vignette_name = "CSGM-Workflow", dest_dir = getwd()) {
  # Find vignette path
  vignette_path <- system.file("doc", paste0(vignette_name, ".Rmd"), package = "CSGM")
  if (vignette_path == "") {
    vignette_path <- system.file("vignettes", paste0(vignette_name, ".Rmd"), package = "CSGM")
  }

  if (vignette_path == "") {
    stop("Vignette source not found in the installed package")
  }

  # Make sure destination directory exists
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)

  # Copy the file
  dest_file <- file.path(dest_dir, basename(vignette_path))
  file.copy(vignette_path, dest_file, overwrite = TRUE)

  message("Vignette extracted to: ", dest_file)
  return(dest_file)
}
