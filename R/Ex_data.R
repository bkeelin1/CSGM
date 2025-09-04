#' @title Ex_data
#'
#' @description
#' This script takes a 3D landmark array and generates raw landmark files in various formats
#' commonly used in geometric morphometrics. The data structure is preserved across
#' all formats.
#'
#' @param Directory a file list indicating the directory that the neseted folder should be created
#'
#' @param Land a p x k x n (landmarks, dimensions, observations) array of a landmark configuration
#'
#' @return a nested list object in a user-specified directory called "Example"
#' containing folders with their corresponding landmark formats.
#'
#' @importFrom readr write_csv
#' @importFrom writexl write_xlsx
#' @importFrom dplyr select rename %>%
#' @importFrom stats na.omit
#'
#' @author Keeling et al., 2025
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(Corpus_Land)
#'
#' Directory <- "E:/Paper 2 - The CSGM Package/CSGM/data/Example"
#'
#' Ex_data(Directory, Corpus_Land)
#'
#' }
Ex_data <- function(Directory, Land) {

  # Define target directory

  # Define CS positions (assuming these are your dimnames)
  CS_positions <- unique(substr(dimnames(Land)[[3]], 1, 6))

  # Function to convert array slice to data frame
  array_to_df <- function(array_slice, specimen_id) {
    # Get landmark names if they exist
    landmark_names <- dimnames(array_slice)[[1]]
    if (is.null(landmark_names)) {
      landmark_names <- paste0("LM", 1:dim(array_slice)[1])
    }

    # Create data frame
    data.frame(
      Name = landmark_names,
      PosX = array_slice[, 1],
      PosY = array_slice[, 2],
      PosZ = array_slice[, 3],
      ID = specimen_id
    )
  }

  # Create directories for each format
  formats <- c("Dragonfly", "Avizo", "Excel", "TPS", "NTS", "Mimics")
  for(format in formats) {
    for(CS in CS_positions) {
      dir.create(file.path(Directory, format, CS),
                 recursive = TRUE, showWarnings = FALSE)
    }
  }

  # Process each CS position
  for(CS in CS_positions) {
    # Get specimens for this CS
    specimens <- grep(CS, dimnames(Land)[[3]], value = TRUE)

    for(specimen in specimens) {
      # Extract specimen data
      spec_data <- Land[, , specimen]
      spec_id <- sub(paste0(CS, "_"), "", specimen)

      # Convert to data frame
      df <- array_to_df(spec_data, spec_id)

      # Save in Dragonfly format (CSV)
      write_csv(df,
                file.path(Directory, "Dragonfly", CS,
                          paste0(spec_id, ".csv")))

      # Save in Avizo format
      avizo_header <- c(
        "# Avizo Landmark ASCII",
        "# Content: 3D landmark coordinates",
        paste("# Specimen:", spec_id),
        "# Format: X Y Z",
        ""
      )

      coords <- select(df, PosX, PosY, PosZ) %>%
        format(scientific = FALSE, trim = TRUE)

      colnames(coords) <- c("df", "X", "Y","Z")

      writeLines(
        c(avizo_header, apply(coords, 1, paste, collapse = " ")),
        file.path(Directory, "Avizo", CS,
                  paste0(spec_id, ".landmarkAscii"))
      )

      # Save in Excel format
      writexl::write_xlsx(df,
                          file.path(Directory, "Excel", CS,
                                    paste0(spec_id, ".xlsx")))

      # Save in NTS format
      writeLines(
        c(as.character(nrow(df)),
          apply(select(df, X, Y, Z), 1, paste, collapse = "\t")),
        file.path(Directory, "NTS", CS,
                  paste0(spec_id, ".nts"))
      )

      # Save in Mimics format
      write_csv(rename(df, Point = Name, X = X, Y = Y, Z = Z),
                file.path(Directory, "Mimics", CS,
                          paste0(spec_id, ".csv")))
    }

    # Save TPS file for each CS position (all specimens in one file)
    tps_content <- NULL
    for(specimen in specimens) {
      spec_data <- Land[, , specimen]
      spec_id <- sub(paste0(CS, "_"), "", specimen)
      coords <- array_to_df(spec_data, spec_id) %>%
        select(PosX, PosY, PosZ)

      specimen_content <- c(
        paste0("LM3=", nrow(coords)),
        paste(format(as.matrix(coords), scientific = FALSE), collapse = "\n"),
        paste0("ID=", spec_id),
        ""
      )
      tps_content <- c(tps_content, specimen_content)
    }
    writeLines(tps_content,
               file.path(Directory, "TPS", CS, "specimens.tps"))
  }

  # Create README file
  readme_content <- c(
    "Landmark Data in Various Data Formats
    =================================================

      This directory contains mandibular corpus landmark data in various formats:

      Dragonfly/ - Landmarks in ORS Dragonfly format
      Avizo/     - Landmarks in Avizo ASCII format
      Excel/     - Landmarks in Excel format
      TPS/       - Landmarks in TPS format
      NTS/       - Landmarks in NTS (Viewbox) format
      Mimics/    - Landmarks in Mimics format

    Each directory contains data for five corpus cross-section positions, in order:
      - Lm1m2 (Left M1-M2 intermolar corpus)
    - Lp3p4 (Left P3-P4 interpremolar corpus)
    - Symphysis (interincisal septum)
    - Rp3p4 (Right P3-P4 interpremolar corpus)
    - Rm1m2 (Right M1-M2 intermolar corpus)

    These landmarks represent mandibular cross-sections of individuals included in
    the development of the CSGM package and discussed in the following article
    (in prep to Nature Communications - Biological Sciences):

      Keeling, Brian A., Conde-Valverde, M.,  Urciuoli, A., Diez-Valero, J., Martínez, I.,
    and Quam, R., 2025. The Biomechanical Nature of the Human Chin. Prepared Manuscript
    to submit to Nature Communications, Biological Sciences."
  )
  writeLines(readme_content, file.path(Directory, "README.md"))
}
