#' @title Import and Process Landmark Data from Various File Formats
#'
#' @description
#' Imports landmark data from various common formats used in geometric morphometrics
#' and medical imaging. Supports multiple file formats including Dragonfly, TPS,
#' Morphologika, NTS (Viewbox), Avizo, Mimics, fcsv, and CSV files. Can handle both
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
#'   \item {fcsv:            3D Slicer Fudicial CSV file format}
#' }
#'
#' @param directory Path to the directory containing landmark files or folders containing landmark files
#' @param pattern File pattern to match (e.g., "\\.csv", "\\.tps")
#' @param format Character string specifying the file format. One of "dragonfly",
#'        "tps", "morphologika", "nts", "avizo", "mimics", "fcsv", or "simple_csv"
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
#'   format = "nts",
#'   dims = 3
#' )
#'
#' # Import fcsv files
#' landmarks <- import_land(
#'   directory = path = file.path(.libPaths("CSGM"), "CSGM", "Example", "fcsv"),
#'   pattern = "\\.fcsv$",
#'   format = "fcsv",
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
  if (!format %in% c("dragonfly", "tps", "morphologika", "nts", "avizo", "mimics", "simple_csv","fcsv")) {
    stop("Unsupported file format")
  }

  # Check for files in the root directory first
  # pattern is now OUTSIDE file.path so it works correctly
  test <- list.files(directory, pattern = pattern, full.names = TRUE)

  # Logic: If files exist in root, set folders to "." (current dir)
  if(length(test) > 0) {
    folders <- "."
  } else {
    # If no files in root, default to subfolders
    if(is.null(folders)) {
      folders <- list.files(directory)
    }
  }

  # Initialize results list
  results <- list()
  p_names = vector("character",0)
  n_names = vector("character", 0)

  for(f in seq_along(folders)) {

    if(folders[f] == ".") {
      search_path <- directory
    } else {
      search_path <- file.path(directory, folders[f])
    }

    files <- list.files(search_path,
                        pattern = pattern,
                        full.names = TRUE)

    # Skip empty folders smoothly
    if (length(files) == 0) {
      if(folders[f] == ".") next
      warning(paste("No matching files in folder:", folders[f]))
      next
    }

    file_ids <- sapply(files, extract_id)

    # Sort files
    if (!any(is.na(suppressWarnings(as.numeric(file_ids))))) {
      files <- files[order(as.numeric(file_ids))]
    } else {
      files <- files[gtools::mixedorder(file_ids)]
    }

    # Read files
    landmarks <- switch(format,
                        "dragonfly" = read_dragonfly(files, dims, curve_names),
                        "tps" = read_tps(files, dims),
                        "morphologika" = geomorph::read.morphologika(files),
                        "nts" = read_nts(files, dims),
                        "avizo" = read_avizo(files, dims),
                        "mimics" = read_mimics(files, dims),
                        "simple_csv" = read_simple_csv(files, dims),
                        "fcsv" = read_fcsv(files, dims)
    )

    if(is.null(landmarks)) next

    if (folders[f] == ".") {
      current_p_names <- dimnames(landmarks)[[1]]
    } else {
      current_p_names <- paste(folders[f], dimnames(landmarks)[[1]], sep = "_")
    }

    p_names = c(p_names, current_p_names)

    # Handle Observation names safely
    obs_names <- if(length(dim(landmarks)) == 3) dimnames(landmarks)[[3]] else colnames(landmarks)
    n_names = c(n_names, obs_names)

    results[[f]] <- data.frame(geomorph::two.d.array(landmarks))
    rm(landmarks)
  }

  if (length(results) == 0) stop("No landmark data could be imported.")

  results <- do.call(rbind, results[!sapply(results, is.null)])

  colnames(results) <- if(dims == 3) {
    paste(rep(c("X","Y","Z"), times = (ncol(results)/3)), rep(1:(ncol(results)/3), each = 3), sep = ".")
  } else {
    paste(rep(c("X","Y"), times = (ncol(results)/2)), rep(1:(ncol(results)/2), each = 2), sep = ".")
  }

  results = geomorph::arrayspecs(results, p = (ncol(results)/dims), k = dims)

  dimnames(results) <- list(Landmarks = p_names,
                            Dimensions = 1:dims,
                            Observations = unique(n_names))

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

#' Read 3D Slicer Fiducial CSV (.fcsv) Format Files
#' @keywords internal
read_fcsv <- function(files, dims) {
  all_landmarks <- list()

  for (file in files) {
    # Read fcsv file.
    # The format has 3 header lines marked with '#'. The data starts on line 4.
    # Column definitions on line 3 specify x, y, z are at indices 2, 3, and 4.
    data <- tryCatch({
      utils::read.csv(file, skip = 3, header = FALSE, stringsAsFactors = FALSE)
    }, error = function(e) {
      warning(paste("Could not read file:", file, "-", e$message))
      return(NULL)
    })

    if (is.null(data)) next

    # Extract coordinates based on dimensions
    if (dims == 3) {
      if (ncol(data) < 4) {
        warning(paste("File", file, "does not have enough columns for 3D data. Skipping."))
        next
      }
      coords <- data[, 2:4]
      colnames(coords) <- c("X", "Y", "Z")
    } else {
      if (ncol(data) < 3) {
        warning(paste("File", file, "does not have enough columns for 2D data. Skipping."))
        next
      }
      coords <- data[, 2:3]
      colnames(coords) <- c("X", "Y")
    }

    coords$ID <- extract_id(file)
    all_landmarks[[length(all_landmarks) + 1]] <- coords
  }

  if (length(all_landmarks) == 0) return(NULL)

  combined_data <- do.call(rbind, all_landmarks) %>%
    dplyr::filter(!is.na(ID)) %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(., .by_group = TRUE)

  # Get dimensions
  ID_list = na.omit(unique(combined_data$ID))
  n = length(ID_list)
  if (n == 0) stop("No valid landmark data found after processing.")
  p = nrow(combined_data %>% dplyr::filter(ID == ID_list[1]))
  k = dims

  # Check for consistency across specimens
  if (nrow(combined_data) != n * p) {
    stop("Landmark data is inconsistent: not all specimens have the same number of landmarks.")
  }

  # Create the 3D array (p x k x n)
  land_data_array = array(dim = c(p, k, n))

  # Fill array from the sorted combined data.
  # matrix(..., nrow = p, ncol = n) creates a p x n matrix where each column corresponds to a specimen.
  land_data_array[, 1, ] <- matrix(combined_data$X, nrow = p, ncol = n)
  land_data_array[, 2, ] <- matrix(combined_data$Y, nrow = p, ncol = n)
  if (k == 3) {
    land_data_array[, 3, ] <- matrix(combined_data$Z, nrow = p, ncol = n)
  }

  # Set dimension names
  dimnames(land_data_array) <- list(
    Landmarks = paste(rep("Land", p), 1:p, sep = "_"),
    Dimensions = c(1:k),
    Observations = ID_list
  )

  return(land_data_array)
}

