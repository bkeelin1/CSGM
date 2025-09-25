#' @title Import and Process Landmark Data from Various File Formats
#'
#' @description
#' Imports landmark data from various common formats used in geometric morphometrics
#' and medical imaging. Supports multiple file formats including Dragonfly, TPS,
#' Morphologika, NTS (Viewbox), Avizo, Mimics, and CSV files. Can handle both
#' 2D and 3D landmarks and provides options for processing curves and contours.
#'
#' \itemize{Supported file formats:
#'   \item {dragonfly        Dragonfly landmark export format}
#'   \item {tps:             TPS file format}
#'   \item {morphologika:    Morphologika format}
#'   \item {nts:             Viewbox/Checkpoint NTS format}
#'   \item {avizo:           Avizo landmark ASCII format}
#'   \item {mimics:          Mimics Medical landmark format}
#'   \item {simple_csv:      Simple CSV with coordinates and optional landmark names}
#' }
#'
#' @param directory Path to the directory containing landmark files or folders containing landmark files
#' @param pattern File pattern to match (e.g., "\\.csv", "\\.tps")
#' @param format Character string specifying the file format. One of "dragonfly",
#'        "tps", "morphologika", "nts", "avizo", "mimics", or "simple_csv"
#' @param dims Integer specifying dimensions (2 or 3)
#' @param curve_names Optional for Dragonfly file formats only. Which specific landmark curve names be extracted?
#'
#' @return A list containing:
#' \itemize{
#'   \item landmarks: A single p x k x n array of all landmarks and individuals Wcombined
#' }
#'
#' @details
#' This function provides a unified interface for mass importing landmark data
#' from various directory folders and file formats commonly used in
#' geometric morphometrics. It accepts both 2D and 3D landmark data and
#' provides options for processing separate curves or contours.
#'
#' This function also allows for a simple CSV format, the function expects either
#' pure coordinate data (x,y,z) or coordinates with landmark names listed in the
#' first column. The coordinate columns can be customized using the coord_names parameter.
#'
#' @seealso
#' \code{\link{landmarks_to_array}}
#'
#' @author Keeling et al., 2025
#'
#' @importFrom readr read_csv
#' @importFrom stringr str_extract str_replace
#' @importFrom dplyr filter mutate arrange group_by
#' @importFrom tidyr %>%
#' @importFrom gtools mixedorder
#' @importFrom utils read.table
#' @importFrom geomorph readmulti.tps read.morphologika two.d.array arrayspecs
#' @importFrom purrr is_empty
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Import Dragonfly landmarks with curves
#'
#' landmarks <- import_land(
#'   directory = path = file.path(.libPaths("CSGM"), "CSGM", "Example", "Dragonfly"), # CSGM library path
#'   pattern = "\\.csv$",
#'   format = "dragonfly",
#'   dims = 3,
#'   curve_names = c("LE","BE","LI","BI")
#' )
#'
#' # Import Avizo files
#' landmarks <- import_land(
#'   directory = path = file.path(.libPaths("CSGM"), "CSGM", "Example", "Avizo"),
#'   pattern = "\\.landmarkAscii$",
#'   format = "avizo",
#'   dims = 3
#' )
#'
#' # Import Mimics files
#' landmarks <- import_land(
#'   directory = path = file.path(.libPaths("CSGM"), "CSGM", "Example", "Mimics"),
#'   pattern = "\\.csv$",
#'   format = "mimics",
#'   dims = 3
#' )
#'
#' # Import simple CSV with landmark names
#' landmarks <- import_land(
#'   directory = path = file.path(.libPaths("CSGM"), "CSGM", "Example", "Excel"),
#'   pattern = "\\.csv$",
#'   format = "simple_csv",
#'   dims = 3,
#' )
#'
#' # Import TPS files
#' landmarks <- import_land(
#'   directory = path = file.path(.libPaths("CSGM"), "CSGM", "Example", "TPS"),
#'   pattern = "\\.tps$",
#'   format = "tps",
#'   dims = 3
#' )
#'
#' # Import NTS files
#' landmarks <- import_land(
#'   directory = path = file.path(.libPaths("CSGM"), "CSGM", "Example", "NTS"),
#'   pattern = "\\.nts$",
#'   format = "tps",
#'   dims = 3
#' )
#' }

import_land <- function(directory,
                        pattern = "\\.csv$",
                        format = "dragonfly",
                        dims = 3,
                        curve_names = NULL,
                        folders = NULL) {

  # Validate inputs
  if (!dir.exists(directory)) stop("Directory does not exist")
  if (!dims %in% c(2, 3)) stop("Dimensions must be 2 or 3")
  if (!format %in% c("dragonfly", "tps", "morphologika", "nts", "avizo", "mimics", "simple_csv")) {
    stop("Unsupported file format")
  }

  test <- list.files(file.path(directory,
                               pattern = pattern,
                               full.names = TRUE))

  if(purrr::is_empty(test)) {
    folders = if(is.null(folders)) list.files(directory) else folders
  } else {
    folders = "Land"
  }


  # Initialize results list
  results <- list()
  p_names = vector("character",0)

  for(f in seq_along(folders)) {

    files <- list.files(file.path(directory, folders[f]),
                        pattern = pattern,
                        full.names = TRUE)

    if (length(files) == 0) stop("No files found matching the pattern")


    file_ids <- sapply(files, extract_id)

    # Sort files based on extracted IDs
    # If all IDs are numeric, sort numerically; otherwise, sort alphanumerically
    if (!any(is.na(suppressWarnings(as.numeric(file_ids))))) {
      files <- files[order(as.numeric(file_ids))]
    } else {
      # For mixed or character IDs, use natural sort
      # This handles cases like "Specimen1", "Specimen2", "Specimen10" correctly
      files <- files[gtools::mixedorder(file_ids)]
    }

    # Read files based on format
    landmarks <- switch(format,
                             "dragonfly" = read_dragonfly(files, dims, curve_names),
                             "tps" = readmulti.tps(files,dims),
                             "morphologika" = read.morphologika(files),
                             "nts" = readmulti.nts(files),
                             "avizo" = read_avizo(files, dims),
                             "mimics" = read_mimics(files, dims),
                             "simple_csv" = read_simple_csv(files, dims)
    )
    p_names = c(p_names, paste(rep(folders[[f]],dim(landmarks)[1]), dimnames(landmarks)[[1]],sep = "_"))
    n_names = dimnames(landmarks)[[3]]
    results[[f]] <- data.frame(geomorph::two.d.array(landmarks))
    rm(landmarks)
  }

  results = data.frame(results)
  colnames(results) <- if(dims ==3) {paste(rep(c("X","Y","Z"), times = (ncol(results)/3)), rep(1:(ncol(results)/3), each = 3), sep = ".")
                       } else { paste(rep(c("X","Y"), times = (ncol(results)/2)), rep(1:(ncol(results)/2), each = 2), sep = ".")}

  results = geomorph::arrayspecs(results, p = (ncol(results)/dims), k = dims)
  dimnames(results) <- list(Landmarks = p_names,
                            Dimensions = 1:dims,
                            Observations = n_names)

  return(results)
}


#' Extract IDs from filenames (handles both numeric and alphanumeric IDs)
#' @keywords internal
extract_id <- function(filepath) {
  # Extract everything between the last slash/backslash and the file extension
  base_name <- str_extract(filepath, "(?<=/|\\\\)[^/\\\\]+(?=\\.[a-zA-Z]+$)")
  # Remove any common prefixes if they exist (like "specimen_" or "ID_")
  clean_id <- str_replace(base_name, "^(specimen_|ID_|Sample_)?", "")
  return(clean_id)
}


#' Read Simple CSV format (X,Y,Z coords)
#' @keywords internal
read_simple_csv <- function(files, dims){
  all_landmarks <- list()

  for (file in files) {
    # Read CSV with specific columns based on dimensions
    cols <- c("X", "Y")
    if (dims == 3) cols <- c(cols, "Z")

    land <- read_csv(file = file,
                     col_select = all_of(cols),
                     show_col_types = FALSE)

    land$ID <- str_extract(file, "(?<=/|\\\\)P?[0-9]+(?=\\.[a-zA-Z]+$)")
    all_landmarks[[length(all_landmarks) + 1]] <- land
  }


  combined_data <- do.call(rbind, all_landmarks) %>% filter(!is.na(ID)) %>% group_by(ID) %>% arrange(., .by_group = TRUE)

  p = nrow(combined_data %>% filter(ID == unique(combined_data$ID)[1]))
  k = dims
  n = length(unique(combined_data$ID))
  ID = na.omit(unique(combined_data$ID))

  curves_data = array(dim = c(p,k,n))

  curves_data[,1,] <- matrix(curve_data2$PosX, nrow = n, ncol = p)
  curves_data[,2,] <- matrix(curve_data2$PosY, nrow = n, ncol = p)
  if(k == 3) {curves_data[,3,] <- matrix(curve_data2$PosZ, nrow = n, ncol = p)}

                             names(curves_data) <- c("Landmarks", "Dimensions", "Observations")
                             dimnames(curves_data) <- list(Landmarks = paste(rep("Land", p), rep(1:p), sep = "_"),
                                                           Dimensions = c(1:k),
                                                           Observations = ID)

  return(curves_data)

}


#' Read Dragonfly Format Files
#' @keywords internal
read_dragonfly <- function(files, dims, curve_names) {
  all_landmarks <- list()

  for (file in files) {
    # Read CSV with specific columns based on dimensions
    cols <- c("Name", "PosX", "PosY")
    if (dims == 3) cols <- c(cols, "PosZ")

    land <- read_csv(file = file,
                     col_select = all_of(cols),
                     show_col_types = FALSE)

    land$ID <- str_extract(file, "(?<=/|\\\\)P?[0-9]+(?=\\.[a-zA-Z]+$)")
    all_landmarks[[length(all_landmarks) + 1]] <- land
  }

  combined_data <- do.call(rbind, all_landmarks)

  if (is.null(curve_names)) {
    curve_names = unique(combined_data$Name)
  }
    # Process separate curves
    curve_data <- list()
    for (curve in curve_names) {
      curve_data[[curve]] <- combined_data %>%
        filter(Name == curve) %>%
        arrange(as.numeric(ID))
    }

    curve_data2 = matrix(NA, nrow = 1, ncol = length(colnames(curve_data[[1]])))
    colnames(curve_data2) <- colnames(curve_data[[1]])
    for(i in seq_along(curve_names)) {
      curve_data[[i]]$ID <- sapply(curve_data[[i]]$ID, as.factor)
      curve_data2 = data.frame(rbind(curve_data2, curve_data[[i]]))
    }
    curve_data2$ID <- sapply(curve_data2$ID, as.factor)
    curve_data2 <- curve_data2 %>% filter(!is.na(ID)) %>% group_by(ID) %>% arrange(.,.by_group = TRUE)

    p = nrow(curve_data2 %>% filter(ID == unique(curves_data2$ID)[1]))
    k = dims
    n = length(na.omit(unique(curves_data2$ID)))
    ID = na.omit(unique(curves_data2$ID))

    land_names = vector("character", 0)
    for(curve in curve_names) {
      psamp = nrow(curves_data2 %>% filter(Name == curve, ID == unique(curves_data2$ID)[1]))
      index = paste(rep(curve, each = psamp), rep(1:psamp), sep = "_")
      land_names = c(land_names, index)
    }

    curves_data = array(dim = c(p,k,n))

    curves_data[,1,] <- matrix(curve_data2$PosX, nrow = n, ncol = p)
    curves_data[,2,] <- matrix(curve_data2$PosY, nrow = n, ncol = p)
    if(k == 3) {curves_data[,3,] <- matrix(curve_data2$PosZ, nrow = n, ncol = p)}

    dimnames(curves_data) <- list(Landmarks = land_names,
                                  Dimensions = c(1:k),
                                  Observations = ID)

    return(curves_data)

}

#' Read TPS Format Files
#' @keywords internal
read_tps <- function(files, dims) {
  if (requireNamespace("geomorph", quietly = TRUE)) {
    landmarks <- geomorph::readmulti.tps(files)
    return(landmarks)
  } else {
    stop("geomorph package required for TPS file import")
  }
}

#' Read NTS (Viewbox) Format Files
#' @keywords internal
read_nts <- function(files, dims) {
  all_landmarks <- list()

  for (file in files) {
    # Read NTS file
    lines <- readLines(file)

    # Extract number of landmarks from header
    n_landmarks <- as.numeric(lines[1])

    # Read coordinates
    coords <- read.table(text = lines[-1], header = FALSE,
                         nrows = n_landmarks)

    # Format depends on dimensions
    if (dims == 2) {
      names(coords) <- c("X", "Y")
    } else {
      names(coords) <- c("X", "Y", "Z")
    }

    coords$ID <- extract_id(file)
    all_landmarks[[length(all_landmarks) + 1]] <- coords
  }

  combined_data <- do.call(rbind, all_landmarks) %>% filter(!is.na(ID)) %>% group_by(ID) %>% arrange(., .by_group = TRUE)

  p = nrow(combined_data %>% filter(ID == unique(combined_data$ID)[1]))
  k = dims
  n = length(unique(combined_data$ID))
  ID = na.omit(unique(combined_data$ID))

  curves_data = array(dim = c(p,k,n))

  curves_data[,1,] <- matrix(curve_data2$PosX, nrow = n, ncol = p)
  curves_data[,2,] <- matrix(curve_data2$PosY, nrow = n, ncol = p)
  if(k == 3) {curves_data[,3,] <- matrix(curve_data2$PosZ, nrow = n, ncol = p)}

                             names(curves_data) <- c("Landmarks", "Dimensions", "Observations")
                             dimnames(curves_data) <- list(Landmarks = paste(rep("Land", p), rep(1:p), sep = "_"),
                                                           Dimensions = c(1:k),
                                                           Observations = ID)

                             return(curves_data)
}

#' Read Avizo Landmark Format Files
#' @keywords internal
read_avizo <- function(files, dims) {
  all_landmarks <- list()

  for (file in files) {
    # Read Avizo landmark file (typically ASCII format)
    lines <- readLines(file)

    # Find the start of coordinate data
    data_start <- which(grepl("^\\s*[0-9.-]", lines))[1]

    # Read coordinates
    coords <- read.table(text = lines[data_start:length(lines)],
                         header = FALSE)

    # Format based on dimensions
    if (dims == 2) {
      names(coords) <- c("X", "Y")
    } else {
      names(coords) <- c("X", "Y", "Z")
    }

    coords$ID <- extract_id(file)
    all_landmarks[[length(all_landmarks) + 1]] <- coords
  }

  combined_data <- do.call(rbind, all_landmarks) %>% filter(!is.na(ID)) %>% group_by(ID) %>% arrange(., .by_group = TRUE)

  p = nrow(combined_data %>% filter(ID == unique(combined_data$ID)[1]))
  k = dims
  n = length(unique(combined_data$ID))
  ID = na.omit(unique(combined_data$ID))

  curves_data = array(dim = c(p,k,n))

  curves_data[,1,] <- matrix(curve_data2$PosX, nrow = n, ncol = p)
  curves_data[,2,] <- matrix(curve_data2$PosY, nrow = n, ncol = p)
  if(k == 3) {curves_data[,3,] <- matrix(curve_data2$PosZ, nrow = n, ncol = p)}

                             names(curves_data) <- c("Landmarks", "Dimensions", "Observations")
                             dimnames(curves_data) <- list(Landmarks = paste(rep("Land", p), rep(1:p), sep = "_"),
                                                           Dimensions = c(1:k),
                                                           Observations = ID)

  return(curves_data)
}

#' Read Mimics Medical Format Files
#' @keywords internal
read_mimics <- function(files, dims) {
  all_landmarks <- list()

  for (file in files) {
    # Read Mimics landmark file (CSV format with specific headers)
    data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)

    # Extract coordinates based on Mimics format
    # Mimics typically uses "X", "Y", "Z" or "x", "y", "z" columns
    coords <- data.frame(
      X = data[[grep("^[Xx]$", names(data))[1]]],
      Y = data[[grep("^[Yy]$", names(data))[1]]]
    )

    if (dims == 3) {
      coords$Z <- data[[grep("^[Zz]$", names(data))[1]]]
    }

    # If landmark names/numbers are present in Mimics file
    if ("Name" %in% names(data) || "Point" %in% names(data)) {
      coords$Name <- data[["Name"]] %||% data[["Point"]]
    }

    coords$ID <- extract_id(file)
    all_landmarks[[length(all_landmarks) + 1]] <- coords
  }

  combined_data <- do.call(rbind, all_landmarks) %>% filter(!is.na(ID)) %>% group_by(ID) %>% arrange(., .by_group = TRUE)

  p = nrow(combined_data %>% filter(ID == unique(combined_data$ID)[1]))
  k = dims
  n = length(unique(combined_data$ID))
  ID = na.omit(unique(combined_data$ID))

  curves_data = array(dim = c(p,k,n))

  curves_data[,1,] <- matrix(curve_data2$PosX, nrow = n, ncol = p)
  curves_data[,2,] <- matrix(curve_data2$PosY, nrow = n, ncol = p)
  if(k == 3) {curves_data[,3,] <- matrix(curve_data2$PosZ, nrow = n, ncol = p)}

                             names(curves_data) <- c("Landmarks", "Dimensions", "Observations")
                             dimnames(curves_data) <- list(Landmarks = paste(rep("Land", p), rep(1:p), sep = "_"),
                                                           Dimensions = c(1:k),
                                                           Observations = ID)

  return(curves_data)
}

#' Convert Landmark Data to Array Format
#'
#' @description
#' Converts a landmark data frame to a 3D array format commonly used in geometric
#' morphometrics analyses. The resulting array has dimensions (landmarks, coordinates, specimens).
#'
#' @details
#' This function takes a data frame of landmark coordinates and converts it to a 3D array
#' format where:
#' \itemize{
#'   \item First dimension represents landmarks
#'   \item Second dimension represents coordinates (2 or 3)
#'   \item Third dimension represents specimens
#' }
#' This format is compatible with many geometric morphometrics packages and analyses.
#'
#' @param landmark_data Data frame containing landmark coordinates
#' @param dims Number of dimensions (2 or 3)
#' @param n_landmarks Number of landmarks per specimen
#'
#' @return A 3D array of landmark coordinates with dimensions (n_landmarks, dims, n_specimens)
#'
#' @examples
#' \dontrun{
#' # Convert imported landmarks to array format
#' landmarks <- import_landmarks("path/to/landmarks", format = "simple_csv")
#' landmark_array <- landmarks_to_array(
#'   landmark_data = landmarks$landmarks,
#'   dims = 3,
#'   n_landmarks = 50
#' )
#' }
#'
#' @seealso
#' \code{\link{import_landmarks}} for importing landmark data
#'
#' @importFrom dplyr filter
#' @importFrom tidyr %>%
#'
#' @export
landmarks_to_array <- function(landmark_data, dims, n_landmarks) {
  n_specimens <- length(unique(landmark_data$ID))

  # Create array with appropriate dimensions
  array_dims <- c(n_landmarks, dims, n_specimens)
  landmark_array <- array(dim = array_dims)

  # Fill array with coordinates
  for (i in 1:n_specimens) {
    spec_data <- landmark_data %>%
      filter(ID == unique(ID)[i])

    landmark_array[,1,i] <- spec_data$PosX
    landmark_array[,2,i] <- spec_data$PosY
    if (dims == 3) landmark_array[,3,i] <- spec_data$PosZ
  }

  return(landmark_array)
}
