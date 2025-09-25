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
#' @author Keeling et al., 2025
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

  # Copy the vignette file
  dest_file <- file.path(dest_dir, basename(vignette_path))
  file.copy(vignette_path, dest_file, overwrite = TRUE)

  # Find and copy image files
  # First, look for images in the vignettes/images directory
  img_path <- system.file("vignettes", "images", package = "CSGM")

  # If not found, try the doc/images directory
  if (img_path == "") {
    img_path <- system.file("doc", "images", package = "CSGM")
  }

  # If images exist, copy them and update the vignette paths if needed
  if (img_path != "") {
    # Create images directory in the destination
    img_dest_dir <- file.path(dest_dir, "images")
    if (!dir.exists(img_dest_dir)) dir.create(img_dest_dir)

    # Copy all image files
    img_files <- list.files(img_path, full.names = TRUE)
    if (length(img_files) > 0) {
      for (img in img_files) {
        file.copy(img, file.path(img_dest_dir, basename(img)), overwrite = TRUE)
      }
      message("Copied ", length(img_files), " image files to ", img_dest_dir)
    }

    # Read the vignette content
    vignette_content <- readLines(dest_file)

    # Check if paths need adjustment (if using absolute paths or different structure)
    if (any(grepl("vignettes/images", vignette_content))) {
      vignette_content <- gsub("vignettes/images", "images", vignette_content)
      writeLines(vignette_content, dest_file)
      message("Updated image paths in the vignette file")
    }
  } else {
    message("No image directory found in the package")
  }

  message("Vignette extracted to: ", dest_file)
  return(dest_file)
}
