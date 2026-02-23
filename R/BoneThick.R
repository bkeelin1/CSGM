#' @title BoneThick - Calculate Bone Thicknesses with Landmark Coordinates
#'
#' @description
#' BoneThick calculates the cortical bone thicknesses of any bone, provided that
#' the landmarks are in a configuration where the external and internal landmarks
#' are sequentially paired in some manner.
#'
#' @param External A landmark object as a "p x k x n" array (landmarks, dimensions, observations), "n x p" 2D matrix or data frame of the external cortical surface.
#' @param Internal A landmark object as a "p x k x n" array (landmarks, dimensions, observations), "n x p" 2D matrix or data frame of the internal cortical surface.
#' @param Points A numeric value representing the total number of landmark points in  you have in a single individual specimen landmark configuration.
#' @param k A numeric value indicating the number of landmark dimensions. Only accepts 2 or 3.
#' @param Spec_ID An optional vector containing the identifying labels for each observation.
#' @param Row_ID An optional logical value indicating whether the External or Internal arguments is a data frame with the first row (not column names) consisting of data column names.
#'
#' @returns An n x p dataframe with the row names identified by Spec_ID which contains landmark euclidean distances between the external and internal surfaces.
#'
#' @details
#' This function takes two inputs: (1) an External landmark configuration as either
#' an "p x k x n" array (landmarks, dimensions, observations), "n x p" 2D matrix
#' or data frame, and (2) an Internal landmark configuraiton in the same configuration.
#' If the landmark shape data is in a dataframe, there are three options for its
#' organization: (1) a dataframe where columns represent a sequence of landmarks
#' and rows represent the individual specimens, (2) a dataframe where the rows
#' represent all of the landmarks (x1,y1,z1,x2,y2,z2) and each column a specimen,
#' or a 2D or 3D array, or (4) the dataframe consists of a .morphologika, .tps, o
#' r .dat format where three columns designates the X, Y, and Z landmark position
#' and consists of all landmarks in a configuration grouping all individuals in
#' these three columns. The function will take both of these configurations and
#' sequentially calculate the euclidean distances between the external and internal
#' landmarks to calculate bone thicknesses. Landmarks should be converted either
#' before the analysis or after to a millimeter scale to ensure meaningful results.
#'
#' @author Keeling et al., 2025
#'
#' @importFrom geomorph two.d.array
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Scenario: Calculate the Cortical bone thicknesses of the mandibular symphysis
#'
#' # Load in the Corpus_Land dataset
#' data(Corpus_Land)
#'
#' # Isolate the External and Internal Symphysis
#' # Symphysis landmarks in Corpus_Land = 241:360
#'
#' External_Symphysis = Corpus_Land[241:300,,]
#' Internal_Symphysis[301:360,,]
#'
#' # Calculate the Bone Thicknesses
#' Thicknesses = BoneThick(External_symphysis, # External landmarks
#'                         Internal_symphysis, # Internal landmarks
#'                         Points = 60, # There are 60 semilandmarks along the external and internal symphyseal bone contour
#'                         k = 3, # three dimensions
#'                         Spec_ID = dimnames(Corpus_Land)[[3]]
#'                         Row_ID = FALSE # The landmark configurations are arrays not data tables
#'                         )
#'
#' # View the Thicknesses as a thickness plot
#'
#' Mean_Thicness = colMeans(Thicknesses) %>%
#'                   pivot_longer(cols = everything(), names_to = "Point",
#'                                 values_to = "Thickness") %>%
#'                   mutate(ID = rep(dimnames(Corpus_Land)[[3]], each = 60))
#'
#' Thickness_long = Thicknesses %>%
#'                   pivot_longer(cols = everything(), names_to = "Point",
#'                                 values_to = "Thickness") %>%
#'                   mutate(ID = rep(dimnames(Corpus_Land)[[3]], each = 60)) %>%
#'
#' ggplot(Thicknesses, aes(x = Point, y = Thickness)) +
#'   geom_line(aes(color = ID)) +
#'   geom_line(data = Mean_Thickness, aes(x = Point, y = Thickness, color = "Black"), size = 5) +
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +  # Rotate x-axis labels
#'   labs(title = "Bone Thickness Across Different Measurement Points", x = "Semilandmark Point", y = "Thickness (mm)") +
#'   ggpubr::theme_pubclean()
#'
#' }

BoneThick <- function(External,
                      Internal,
                      Points,
                      k,
                      Spec_ID = NULL,
                      Row_ID = FALSE) {

  if (k == 3) { # for 3D landmark configurations

    if (is.array(External) || is.array(Internal)) {
      External = data.frame(geomorph::two.d.array(External))
      Internal = data.frame(geomorph::two.d.array(Internal))
    } else {}

    if (nrow(External) == Points*3) {
      External <- t(External)
    } else {}

    if (isTRUE(Row_ID)) {
      External <- External[-1,]
      Internal <- Internal[-1,]
    } else {}

    if (ncol(External) <= 3) {
      ID = nrows(External) / Points
      n <- c(rep(1:ID, each = length(Points)))
      p <- c(rep(1:Points, times = ID))
      matrix1 <- matrix(Procrustes_data[,1], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      matrix2 <- matrix(Procrustes_data[,2], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      matrix3 <- matrix(Procrustes_data[,3], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      External_array <- array(dim = c(length(ID), k, nrow(External)/length(ID)))


      External_array[,1,] <- matrix1
      External_array[,2,] <- matrix2
      External_array[,3,] <- matrix3
      External = data.frame(geomorph::two.d.array(External_array))


      ID = nrows(Internal) / Points
      n <- c(rep(1:ID, each = length(Points)))
      p <- c(rep(1:Points, times = ID))
      matrix1 <- matrix(Internal[,1], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      matrix2 <- matrix(Internal[,2], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      matrix3 <- matrix(Internal[,3], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      External_array <- array(dim = c(length(ID), k, nrow(External)/length(ID)))


      External_array[,1,] <- matrix1
      External_array[,2,] <- matrix2
      External_array[,3,] <- matrix3
      Internal = data.frame(geomorph::two.d.array(Internal_array))

    } else {}

    Obs <- nrow(External)

    Ext_Points <- vector("list", length = ncol(External)/3)

    for (i in seq(1, ncol(External), by = 3)) {
      subset_1 <- External[, i:(i + 2), drop = FALSE]
      Ext_Points[[i/3 + 1]] <- subset_1
    }

    for (i in 1:length(Ext_Points)) {
      for (j in 1:length(Ext_Points[[i]])) {
        names(Ext_Points[[i]]) <- c("X", "Y", "Z")
        Ext_Points[[i]][[j]] <- as.numeric(Ext_Points[[i]][[j]])
      }
    }


    Int_Points <- vector("list", length = ncol(Internal)/3)

    for (i in seq(1, ncol(Internal), by = 3)) {
      subset_2 <- Internal[, i:(i + 2), drop = FALSE]
      Int_Points[[i/3 + 1]] <- subset_2
    }


    for (i in 1:length(Int_Points)) {
      for (j in 1:length(Int_Points[[i]])) {
        names(Int_Points[[i]]) <- c("X", "Y", "Z")
        Int_Points[[i]][[j]] <- as.numeric(Int_Points[[i]][[j]])
      }
    }

    Thickness <- vector("list", Points)
    Distance <- data.frame(rep(NA, Obs))

    for (i in 1:Points) {
      Distance[i] <-  sqrt((as.numeric(Ext_Points[[i]]$'X') - as.numeric(Int_Points[[i]]$'X'))^2 +
                             (as.numeric(Ext_Points[[i]]$'Y') - as.numeric(Int_Points[[i]]$'Y'))^2 +
                             (as.numeric(Ext_Points[[i]]$'Z') - as.numeric(Int_Points[[i]]$'Z'))^2)
      Thickness[[i]] <- Distance[i]
    }

    Thickness <- unlist(Thickness)

    Thickness <- data.frame(matrix(Thickness, nrow = Obs, ncol = Points))

    colnames(Thickness) <- paste("Point", 1:Points, sep = "")

    if (is.null(Spec_ID)) {
    } else {
      Spec_ID <- unlist(Spec_ID)
      rownames(Thickness) <- Spec_ID
    }

    return(Thickness)
  }

  else { # for 2D landmark configurations

    if (is.array(External) || is.array(Internal)) {
      External <- matrix(External, nrow = dim(External)[3], ncol = dim(External)[1])
      External <- data.frame(External)
      Internal <- matrix(Internal, nrow = dim(Internal)[3], ncol = dim(Internal)[1])
      Internal <- data.frame(Internal)
    } else {}

    if (nrow(External) == Points*2) {
      External <- t(External)
    } else {}

    if (isTRUE(Row_ID)) {
      External <- External[-1,]
      Internal <- Internal[-1,]
    } else {}

    if (ncol(External) == 2) {
      ID = nrows(External) / Points
      n <- c(rep(1:ID, each = length(Points)))
      p <- c(rep(1:Points, times = ID))
      matrix1 <- matrix(Procrustes_data[,1], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      matrix2 <- matrix(Procrustes_data[,2], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      External_array <- array(dim = c(ID, k, Points))

      # Fill the array
      External_array[,1,] <- matrix1
      External_array[,2,] <- matrix2
      External = data.frame(geomorph::two.d.array(External_array))


      ID = nrows(Internal) / Points
      n <- c(rep(1:ID, each = length(Points)))
      p <- c(rep(1:Points, times = ID))
      matrix1 <- matrix(Internal[,1], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      matrix2 <- matrix(Internal[,2], nrow = ID, ncol = Points,
                        dimnames = list(unique(n), unique(p)))
      Internal_array <- array(dim = c(ID, k, Points))

      # Fill the array
      Internal_array[,1,] <- matrix1
      Internal_array[,2,] <- matrix2
      Internal = data.frame(geomorph::two.d.array(Internal_array))
    } else {}

    Obs <- nrow(External)

    Ext_Points <- vector("list", length = ncol(External)/2)

    for (i in seq(1, ncol(External), by = 2)) {
      subset_1 <- External[, i:(i + 1), drop = FALSE]
      Ext_Points[[i/2 + 1]] <- subset_1
    }

    for (i in 1:length(Ext_Points)) {
      for (j in 1:length(Ext_Points[[i]])) {
        names(Ext_Points[[i]]) <- c("X", "Y")
        Ext_Points[[i]][[j]] <- as.numeric(Ext_Points[[i]][[j]])
      }
    }


    Int_Points <- vector("list", length = ncol(Internal)/2)

    for (i in seq(1, ncol(Internal), by = 2)) {
      subset_2 <- Internal[, i:(i + 1), drop = FALSE]
      Int_Points[[i/2 + 1]] <- subset_2
    }


    for (i in 1:length(Int_Points)) {
      for (j in 1:length(Int_Points[[i]])) {
        names(Int_Points[[i]]) <- c("X", "Y")
        Int_Points[[i]][[j]] <- as.numeric(Int_Points[[i]][[j]])
      }
    }

    Thickness <- vector("list", Points)
    Distance <- data.frame(rep(NA, Obs))

    for (i in 1:Points) {
      Distance[i] <-  sqrt((as.numeric(Ext_Points[[i]]$'X') - as.numeric(Int_Points[[i]]$'X'))^2 +
                             (as.numeric(Ext_Points[[i]]$'Y') - as.numeric(Int_Points[[i]]$'Y'))^2)
      Thickness[[i]] <- Distance[i]
    }

    Thickness <- unlist(Thickness)

    Thickness <- data.frame(matrix(Thickness, nrow = Obs, ncol = Points))

    colnames(Thickness) <- paste("Point", 1:Points, sep = "")

    if (is.null(Spec_ID)) {
    } else {
      Spec_ID <- unlist(Spec_ID)
      rownames(Thickness) <- Spec_ID
    }

    return(Thickness)
  }
}
