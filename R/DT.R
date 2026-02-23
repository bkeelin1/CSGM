#' @title Data Transformation (DT) - A function to transform data
#'
#' @description This function takes an array, matrix, or data table object and
#' transforms the data into either a distance matrix, principal components,
#' geometric mean, or 2D data table. This function can also normalize the data
#' and allows for subsetting of data.
#'
#' @param Response an array, matrix, or data table for the response(dependent) variables.
#'
#' @param Predictor an array, matrix, or data table for the predictor(independent) variables.
#'
#' #' @param Res_transform Character string indicating which data transformation of the response
#' (i.e., landmark array) variable is desired. There are several acceptable options for data
#' transformation of the response variable:
#' \itemize{
#'   \item procdist # Procrustes or other distances matrix
#'   \itemize{
#'     \item Calculates Procrustes or other distances as sqrt(sum(squared differences))
#'     \item Accepts array, matrix, or data table
#'   }
#'   \item pca # Principal Component Analysis
#'   \itemize{
#'     \item Direct PCA on array, matrix, or data table
#'   }
#'   \item procdist+pca # Sequential Procrustes or other distances + PCA
#'   \itemize{
#'     \item Calculates Procrustes or other distances
#'     \item PCA on resulting distance matrix
#'   }
#'   \item GM # Geometric Mean
#'   \itemize{
#'     \item Calculates geometric mean of all variables as exp(mean(log(x + 1e-8)))
#'     \item For arrays: calculates landmark configuration centroid
#'   }
#'   \item procdist+GM # Sequential Procrustes or other distances + Geometric Mean
#'   \itemize{
#'     \item Calculates Procrustes or other distances
#'     \item Geometric mean on resulting distances
#'   }
#'   \item 2D # Dimension reduction
#'   \itemize{
#'     \item Transforms 3D array to data frame
#'   }
#'   \item NULL # No transformation
#'   \itemize{
#'     \item Data remains unchanged
#'   }
#' }
#'
#' @param Pred_transform Character string indicating which data transformation of the predictor
#' (i.e., landmark array) variable is desired. There are several acceptable options for data
#' transformation of the response variable:
#' \itemize{
#'   \item procdist # Procrustes or other distances matrix
#'   \itemize{
#'     \item Calculates Procrustes or other distances as sqrt(sum(squared differences))
#'     \item Accepts array, matrix, or data table
#'   }
#'   \item pca # Principal Component Analysis
#'   \itemize{
#'     \item Direct PCA on array, matrix, or data table
#'   }
#'   \item procdist+pca # Sequential Procrustes or other distances + PCA
#'   \itemize{
#'     \item Calculates Procrustes or other distances
#'     \item PCA on resulting distance matrix
#'   }
#'   \item GM # Geometric Mean
#'   \itemize{
#'     \item Calculates geometric mean of all variables as exp(mean(log(x + 1e-8)))
#'     \item For arrays: calculates landmark configuration centroid
#'   }
#'   \item procdist+GM # Sequential Procrustes or other distances + Geometric Mean
#'   \itemize{
#'     \item Calculates Procrustes or others distances
#'     \item Geometric mean on resulting distances
#'   }
#'   \item 2D # Dimension reduction
#'   \itemize{
#'     \item Transforms 3D array to data frame
#'   }
#'   \item NULL # No transformation
#'   \itemize{
#'     \item Data remains unchanged
#'   }
#' }
#'
#' @param Res_ncomp Optional integer. Number of components to retain for response variables.
#' Only applicable when Res_transform string includes PCA transformation. Default is NULL.
#'
#' @param Pred_ncomp Optional integer. Number of components to retain for predictor variables.
#' Only applicable when Pred_transform string includes PCA transformation. Default is NULL.
#'
#' @param subset an optional logical indicating whether only a subset of the variables
#' in either the response and/or predictor variables should be considered in the
#' data transformation.
#'
#' @param Res_point_set an optional vector containing the numeric representation of variables
#' in the data table, matrix, or array to subset the response variables for this
#' data tranformation. If a subset is provided, the function will only perform
#' the data transformation on the selected variables.
#'
#' @param Pred_point_set an optional vector containing the numeric representation of variables
#' in the data table, matrix, or array to subset the predictor variables for this
#' data tranformation. If a subset is provided, the function will only perform
#' the data transformation on the selected variables.
#'
#' @param dt_method an optional character string used only when data transformation
#' for the response and/or predictor variable is set to an option with "procdist".
#' The default is "procdist". If default is changed, the function will overwrite the
#' procdist iteration and perform a distance matrix from the following methods:
#' "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard",
#' "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao",
#' "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or
#' "robust.aitchison". Please see the *vegdist* "method" argument from the
#' **vegan** package for more details on distance matrix calculations.
#'
#' @details This is an internal function was designed to transform variables prior to performing
#' a regression. It accepts either a p x k x n (variables, dimensions, observations)
#' array, a matrix, or data table object as the Response and/or Predictor variables.
#' This function allows for data transformation of subsetted data by switching
#' the argument "subset" to TRUE and providing a numeric vector of the numeric
#' representations (i.e., column number) of the variables of interest in the data
#' object. This function provides a number of data transformation types such as:
#' "procdist" (procrustes or another distance matrix), "procdist+pca" (Procrustes
#' or another distance matrix then subjected to principal component analysis),
#' "procdist+GM" (Procrustes or another distance matrix then the Geometric Mean
#' is calculated), "pca" (Principal Component Analysis), "GM" (Geometric Mean),
#' 2D (Change a 3D array to a 2D data table), NULL (do not transform data). Each
#' of these options are common methods to transform geometric morphometric or
#' biological variables for downstream analysis. If 'procdist', or any option with
#' 'procdist', is selected the method argument will be used to select which
#' distance matrix calculation will be used to trasnform the data. The option
#' 'procdist' or Procrustes distances will be selected by default unless changed.
#' Procrustes distances are calculated as the square root sum of all squared distances.
#' All other distance matrix calculations use the *vegdist* function from the
#' **vegan** package. Principal Component Analysis uses the *gm.prcomp* function
#' from the **geomorph** package for 3D arrays. For 2D data tables, the *PCA*
#' function from **factoextra** is applied.
#'
#' @returns a nested list object with the following structure:
#'
#' \itemize{
#'    \item{Response:   Optional list object with transformed data (Response). If a variant of "pca" is selected, several objects will be included: Res_dim (original data dimensions), Res_p (number of original variables), Res_pca (principal component analysis object), Res_rot (the rotation matrix of the PCA), Res_pca_scores (the principal component scores). Used in other CSGM functions}
#'
#'    \item{Predictor:  Optional list object with transformed data (Predictor).}
#' }
#'
#' @author Keeling et al., 2025
#'
#' @seealso \code{gm.prcomp} \code{PCA}
#'
#' @importFrom geomorph arrayspecs two.d.array gm.prcomp
#' @importFrom stats mahalanobis cov var
#' @importFrom vegan vegdist
#' @importFrom FactoMineR PCA
#' @importFrom factoextra get_pca
#' @importFrom ade4 dudi.pca bca
#' @importFrom dplyr filter mutate
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Scenario: Create a Procrustes distance matrix of a landmark configuration
#' # and find the geometric mean for all corresponding biomechanical variables.
#'
#' # Import data sample of mandibular corpus: Corpus_Land
#' # Import data sample of mandibular corpus data table: Biomechanics
#'
#' data(Corpus_Land)
#' data(Biomechanics)
#'
#' # Conduct the data transformation
#'
#' # NOTE: only numeric values can be taken for geometric means for the biomechanics dataset
#' # NOTE: Use subset for Biomechanic data table to select only numeric values
#'
#' output <- DT(Response = Corpus_Land,
#'              Predictor = Biomechanics,
#'              Res_transformation = "procdist",
#'              Pred_transformation = "GM",
#'              subset = TRUE,
#'              point_set_res = NULL, # The response variable is not subset.
#'              point_set_pred = c(5:93) # only numeric variables
#'              ) # all other arguments leave as default
#'
#' Land_matrix = output$Response # our Procrustes distances
#'
#' GM_Bio = output$Predictor # our Geometric Means of all biomechanical variables
#'
#' # Add identifying information back into GM_Bio
#'
#' GM_Bio = data.frame(Biomechanics[,1:4], GM_Bio)
#'
#' # Since the Biomechanic data has 5 cross-sections and duplicated IDs, we need to
#' # restructure the data so each cross-section is a variable of a geometric mean
#' # for each individual.
#'
#' GM_Bio = GM_Bio %>% summarise('ID' = unique(GM_Bio[,1]),
#'                               'LM1M2' = GM_Bio[,5] %>% filter(Region == "LM1M2"),
#'                               'LP3P4' = GM_Bio[,5] %>% filter(Region == "LP3P4"),
#'                               'Symphysis' = GM_Bio[,5] %>% filter(Region == "Symphysis"),
#'                               'RP3P4' = GM_Bio[,5] %>% filter(Region == "RP3P4"),
#'                               'RM1M2' = GM_Bio[,5] %>% filter(Region == "RM1M2")
#'                              )
#' # view the data
#'
#' view(Land_matrix)
#'
#' view(GM_Bio)
#' }
DT <- function(Response = NULL,
               Predictor = NULL,
               Res_transform = "NULL",
               Pred_transform = "NULL",
               Res_ncomp = NULL,
               Pred_ncomp = NULL,
               Res_point_set = NULL,
               Pred_point_set = NULL,
               dt_method = "procdist") {

  Output = list()

  IDT <- function(Dataset,
                  transform = "NULL",
                  ncomp = NULL,
                  point_set = NULL,
                  token = FALSE,
                  method = "procdist") {

        DT_Output = list()

        if(length(dim(Dataset)) == 3) {
            if(!is.null(point_set)) {
              Dataset = Dataset[point_set,,]
            }

          if(startsWith(transform, "procdist")) {
            #print(paste("Dataset transformation via Procrustes Distances of landmarks selected."))

            if(method == "mahalanobis") {
              Dataset = data.frame(geomorph::two.d.array(Dataset))
              if (ncol(Dataset) >= nrow(Dataset)) {
                PCA = gm.prcomp(Dataset)
                PCA$A = Dataset

                if(is.null(ncomp)){
                  npc = PCA_variances(PCA)
                  npc = npc$PCs
                  Dataset = PCA$x[,1:npc]
                } else {
                  Dataset = PCA$x[,1:ncomp]
                }
                Mean = colMeans(Dataset)
                Cov = cov(Dataset)
                Euc_dist = mahalanobis(Dataset, center = Mean, cov = Cov)
                Dataset_dist = dist(Euc_dist)
                Dataset = as.dist(Dataset_dist)

              } else {
                Euc_dist = vegdist(Dataset, method = method)
                Dataset_dist = Euc_dist
                Dataset = as.dist(Dataset_dist)
              }
            } # end of mahalanobis distances

            if(method != "procdist" & method != "mahalanobis") {
              Dataset = data.frame(geomorph::two.d.array(Dataset))
              Euc_dist = vegdist(Dataset, method = method)
              Dataset_dist = Euc_dist
              Dataset = as.dist(Dataset_dist)
              colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")
            }

            if(method == "procdist") {

                place = matrix(NA, ncol = dim(Dataset)[1], nrow = dim(Dataset)[3])
                place = data.frame(place)

                mean_shape = array(dim = c(dim(Dataset)[1],dim(Dataset)[2], 1))
                Res = Dataset

                mean_shape[,1,] <- mean(Res[,1,])
                mean_shape[,2,] <- mean(Res[,2,])
                mean_shape[,3,] <- mean(Res[,3,])

                Res = data.frame(two.d.array(Res))
                D1 <- vector("list", length = ncol(Res)/3)
                for (i in seq(1, ncol(Res), by = 3)) {
                  subset_2 <- Res[, i:(i + 2), drop = FALSE]
                  D1[[i/3 + 1]] <- subset_2
                }
                for (i in 1:length(D1)) {
                  for (j in 1:length(D1[[i]])) {
                    names(D1[[i]]) <- c("X", "Y", "Z")
                    D1[[i]][[j]] <- as.numeric(D1[[i]][[j]])
                  }
                }

                mean = data.frame(two.d.array(mean_shape))
                D2 <- vector("list", length = ncol(mean)/3)
                for (i in seq(1, ncol(mean), by = 3)) {
                  subset_2 <- mean[, i:(i + 2), drop = FALSE]
                  D2[[i/3 + 1]] <- subset_2
                }

                for (i in 1:length(D2)) {
                  for (j in 1:length(D2[[i]])) {
                    names(D2[[i]]) <- c("X", "Y", "Z")
                    D2[[i]][[j]] <- as.numeric(D2[[i]][[j]])
                  } # end of j for loop
                } # end of i for loop

                for (j in 1:length(D1)) {
                  place[,j] <-  sqrt((as.numeric(D1[[j]]$'X') - as.numeric(D2[[j]]$'X'))^2 +
                                       (as.numeric(D1[[j]]$'Y') - as.numeric(D2[[j]]$'Y'))^2 +
                                       (as.numeric(D1[[j]]$'Z') - as.numeric(D2[[j]]$'Z'))^2)
                } # end of j for loop

                Dataset <- place
                colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")
                scale  = FALSE
            } # end of procdist method
          } else if (startsWith(transform, "pca")) {
            #print(paste("Dataset transformation of Dataset data via PCA of data selected."))

            Res_o = Dataset
            Res_dim = dim(Dataset)[2]
            Res_p = dim(Dataset)[1]
            Res_pca <- gm.prcomp(Dataset)
            Res_rot = Res_pca$rotation
            Res_pca_scores = data.frame(Res_pca$x)
            Res_pca_center = Res_pca$center
            Dataset <- data.frame(Res_pca$x)
            if(!is.null(ncomp)) {
              if(length(1:ncomp) > ncol(Dataset)) {

              } else {
                Dataset <- Dataset[,1:as.numeric(ncomp)]
              }
            }
            scale = FALSE

          } else if (startsWith(transform, "bgpca")) {
            #print(paste("Dataset transformation of 3D array data via BGPCA of data selected."))

            Dataset = data.frame(geomorph::two.d.array(Dataset))
            Res_o = Dataset
            Res_dim = dim(Dataset)[2]
            Res_p = dim(Dataset)[1]
            Res_pca = ade4::dudi.pca(Dataset, scannf = FALSE, nf = 1000)
            Res_rot = Res_pca$c1
            Res_pca_center = as.matrix(bca_output$cent)
            scale = FALSE
            # conduct bgpca
            if(!is.null(ncomp)) {
              if(length(1:ncomp) < ncol(Dataset)) {
                Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = max(ncomp))
              } else {
                Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = 1000)
              }
            } else {
              Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = 1000)
            }
            Res_pca_scores = as.matrix(Dataset$ls) %*% t(Dataset$c1)

          } else if (startsWith(transform, "GM")) {
            calculate_centroid <- function(landmarks) {
              apply(landmarks, 2, mean)
            }
            Dataset <- apply(Dataset, 3, calculate_centroid)
            scale = FALSE

          } else if (startsWith(transform, "2D")) {
            #print(paste("No Dataset transformation of a 3D Dataset array selected.","To avoid rescaling, make columns names contain: Dist, .X, Dim, CS, or Land."))
            Dataset <- data.frame(geomorph::two.d.array(Dataset))
            if(any(grepl("Comp", colnames(Dataset))) | any(grepl("Dist", colnames(Dataset))) |
               any(grepl(".X", colnames(Dataset))) | any(grepl("Dim", colnames(Dataset))) |
               any(grepl("Land", colnames(Dataset))) | any(grepl("CS", colnames(Dataset))) |
               any(grepl("*\\.1", colnames(Dataset)))) {
              scale = FALSE
            } else {
              scale = TRUE
            }
            scale = FALSE

          } else if (transform == "NULL") {
            print(paste("No Dataset transformation of a 3D array selected."))
            Dataset <- Dataset
            scale = FALSE
          }
          #_________________________________________________________________________________________
          # Now perform possible subsequent combinations

          if(endsWith(transform, "+procdist")) {
            #print(paste("Dataset transformation via Procrustes Distances of landmarks selected."))

            if(method == "mahalanobis") {
              Dataset = data.frame(Dataset)
              if (ncol(Dataset) >= nrow(Dataset)) {
                PCA = geomorph::gm.prcomp(Dataset)
                PCA$A = Dataset
                if(is.null(ncomp)){
                  npc = PCA_variances(PCA)
                  npc = npc$PCs
                  Dataset = PCA$x[,1:npc]
                } else {
                  Dataset = PCA$x[,1:ncomp]
                }
                Mean = colMeans(Dataset)
                Cov = cov(Dataset)
                Euc_dist = mahalanobis(Dataset, center = Mean, cov = Cov)
                Dataset_dist = dist(Euc_dist)
                Dataset = as.dist(Dataset_dist)

              } else {
                Euc_dist = vegan::vegdist(Dataset, method = method)
                Dataset_dist = Euc_dist
                Dataset = as.dist(Dataset_dist)
              }
            } # end of mahalanobis distances

            if(method != "procdist" & method != "mahalanobis") {
              Euc_dist = vegan::vegdist(Dataset, method = method)
              Dataset_dist = Euc_dist
              Dataset = as.dist(Dataset_dist)
              #colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")
            }

            if(method == "procdist") {
              Eucdist = matrix(nrow = nrow(Dataset), ncol = ncol(Dataset))
              m_pred = data.frame(mean = colMeans(Dataset))
              Eucdist = data.frame(sweep(Dataset, 2, as.numeric(m_pred[[1]]), "-"))
              Dataset = Eucdist
              Res_o = Eucdist
              colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")
            } # end of procdist method
            scale = FALSE
            #colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")

          } else if (endsWith(transform, "+pca")) {
            #print(paste("Dataset transformation via PCA of data selected."))

            Res_o = Dataset
            Res_dim = dim(Dataset)[2]
            Res_p = dim(Dataset)[1]
            Res_pca <- geomorph::gm.prcomp(Dataset)
            Res_rot = Res_pca$rotation
            Res_pca_center = Res_pca$center
            Res_pca_scores = data.frame(Res_pca$x)
            Dataset <- data.frame(Res_pca$x)
            if(!is.null(ncomp)) {
              if(length(1:ncomp) > ncol(Dataset)) {

              } else {
                Dataset <- Dataset[,as.numeric(1:ncomp)]
              }
            }
            scale = FALSE

          } else if (endsWith(transform, "+bgpca")) {
            #print(paste("Dataset transformation via BGPCA of data selected."))

            Res_o = Dataset
            Res_dim = dim(Dataset)[2]
            Res_p = dim(Dataset)[1]
            Res_pca <- ade4::dudi.pca(Dataset, scannf = FALSE, nf = 1000)
            Res_rot = Res_pca$c1
            Res_pca_center = as.matrix(bca_output$cent)
            scale = FALSE
            # conduct bgpca
            if(!is.null(ncomp)) {
              if(length(1:ncomp) < ncol(Dataset)) {
                Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = max(ncomp))
              } else {
                Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = 1000)
              }
            } else {
              Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = 1000)
            }
            Res_pca_scores = as.matrix(Dataset$ls) %*% t(Dataset$c1)

          } else if (endsWith(transform, "+GM")) {
            #print(paste("Dataset transformation via Geometric Mean of variables selected"))
            Dataset <- apply(Dataset, 1, function(x) exp(mean(log(x+1e-8))))
            Dataset <- data.frame(GM = Dataset)
            scale = FALSE
          }

          #_________________________________________________________________________________________
          # Now perform starting combinations of 2D data tables or matrices

        } else {

          if(!is.null(point_set)) {
            Dataset = Dataset[,point_set]
          }

          if(startsWith(transform, "procdist")) {
            #print(paste("Dataset transformation via Procrustes Distances of landmarks selected."))

            if(method == "mahalanobis") {
              Dataset = data.frame(Dataset)
              if (ncol(Dataset) >= nrow(Dataset)) {
                PCA = geomorph::gm.prcomp(Dataset)
                PCA$A = Dataset
                if(is.null(ncomp)){
                  npc = PCA_variances(PCA)
                  npc = npc$PCs
                  Dataset = PCA$x[,1:npc]
                } else {
                  Dataset = PCA$x[,1:ncomp]
                }
                Mean = colMeans(Dataset)
                Cov = cov(Dataset)
                Euc_dist = mahalanobis(Dataset, center = Mean, cov = Cov)
                Dataset_dist = dist(Euc_dist)
                Dataset = as.dist(Dataset_dist)

              } else {
                Euc_dist = vegan::vegdist(Dataset, method = "mahalanobis")
                Dataset_dist = Euc_dist
                Dataset = as.dist(Dataset_dist)
              }
            } # end of mahalanobis distances

            if(method != "procdist" & method != "mahalanobis") {
              Euc_dist = vegan::vegdist(Dataset, method = method)
              Dataset_dist = Euc_dist
              Dataset = as.dist(Dataset_dist)
              #colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")
            }

            if(method == "procdist") {
              Eucdist = matrix(nrow = nrow(Dataset), ncol = ncol(Dataset))
              m_pred = data.frame(mean = colMeans(Dataset))
              Eucdist = data.frame(sweep(Dataset, 2, as.numeric(m_pred[[1]]), "-"))
              Dataset = Eucdist
              colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")
            } # end of procdist method
            scale = FALSE

          } else if (startsWith(transform, "pca")) {
            #print(paste("Dataset transformation via PCA of data selected."))

            Res_o = data.frame(Dataset)
            Res_dim = dim(Dataset)[2]
            Res_p = dim(Dataset)[1]
            Res_pca <- geomorph::gm.prcomp(Dataset)
            Res_rot = Res_pca$rotation
            Res_pca_scores = data.frame(Res_pca$x)
            Res_pca_center = Res_pca$center
            Dataset <- data.frame(Res_pca$x)
            if(!is.null(ncomp)) {
              if(length(1:ncomp) > ncol(Dataset)) {

              } else {
                Dataset <- Dataset[,1:as.numeric(ncomp)]
              }
            }
            scale = FALSE

          } else if (startsWith(transform, "bgpca")) {
            #print(paste("Dataset transformation via BGPCA of data selected."))

            Res_o = data.frame(Dataset)
            Res_dim = dim(Dataset)[2]
            Res_p = dim(Dataset)[1]
            Res_pca <- ade4::dudi.pca(Dataset, scannf = FALSE, nf = 1000)
            Res_rot = Res_pca$c1
            Res_pca_center = as.matrix(bca_output$cent)
            scale = FALSE
            # conduct bgpca
            if(!is.null(ncomp)) {
              if(length(1:ncomp) < ncol(Dataset)) {
                Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = max(ncomp))
              } else {
                Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = 1000)
              }
            } else {
              Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = 1000)
            }
            Res_pca_scores = as.matrix(Dataset$ls) %*% t(Dataset$c1)

          } else if (startsWith(transform, "GM")) {
            #print(paste("Dataset transformation via Geometric Mean of variables selected"))
            Dataset <- apply(Dataset, 1, function(x) exp(mean(log(x+1e-8))))
            Dataset <- data.frame(GM = Dataset)
            scale = FALSE

          } else if (transform == "NULL") {
            #print(paste("No Dataset transformation of a 3D array selected."))
            Dataset <- Dataset
            scale = FALSE
          }
          #_________________________________________________________________________________________
          # Now perform possible subsequent combinations

          if(endsWith(transform, "+procdist")) {
            #print(paste("Dataset transformation via Procrustes Distances of landmarks selected."))

            if(method == "mahalanobis") {
              Dataset = data.frame(Dataset)
              if (ncol(Dataset) >= nrow(Dataset)) {
                PCA = geomorph::gm.prcomp(Dataset)
                PCA$A = Dataset
                if(is.null(ncomp)){
                  npc = PCA_variances(PCA)
                  npc = npc$PCs
                  Dataset = PCA$x[,1:npc]
                } else {
                  Dataset = PCA$x[,1:ncomp]
                }
                Mean = colMeans(Dataset)
                Cov = cov(Dataset)
                Euc_dist = mahalanobis(Dataset, center = Mean, cov = Cov)
                Dataset_dist = dist(Euc_dist)
                Dataset = as.dist(Dataset_dist)

              } else {
                Euc_dist = vegdist(Dataset, method = method)
                Dataset_dist = Euc_dist
                Dataset = as.dist(Dataset_dist)
              }
            } # end of mahalanobis distances

            if(method != "procdist" & method != "mahalanobis") {
              Euc_dist = vegan::vegdist(Dataset, method = method)
              Dataset_dist = Euc_dist
              Dataset = as.dist(Dataset_dist)
              #colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")
            }

            if(method == "procdist") {
              Eucdist = matrix(nrow = nrow(Dataset), ncol = ncol(Dataset))
              m_pred = data.frame(mean = colMeans(Dataset))
              Eucdist = data.frame(sweep(Dataset, 2, as.numeric(m_pred[[1]]), "-"))
              Dataset = Eucdist
              Res_o = Eucdist
              colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")
            } # end of procdist method
            scale = FALSE
            #colnames(Dataset) <- paste(rep("Dist",ncol(Dataset)), rep(1:ncol(Dataset)), sep = ".")

          } else if (endsWith(transform, "+pca")) {
            #print(paste("Dataset transformation via PCA of data selected."))

            Res_o = Dataset
            Res_dim = dim(Dataset)[2]
            Res_p = dim(Dataset)[1]
            Res_pca <- geomorph::gm.prcomp(Dataset)
            Res_rot = Res_pca$rotation
            Res_pca_scores = data.frame(Res_pca$x)
            Res_pca_center = Res_pca$center
            Dataset <- data.frame(Res_pca$x)
            if(!is.null(ncomp)) {
              if(length(1:ncomp) > ncol(Dataset)) {

              } else {
                Dataset <- Dataset[,1:as.numeric(ncomp)]
              }
            }
            scale = FALSE

          } else if (endsWith(transform, "+bgpca")) {
            #print(paste("Dataset transformation via BGPCA of data selected."))

            Res_o = Dataset
            Res_dim = dim(Dataset)[2]
            Res_p = dim(Dataset)[1]
            Res_pca <- ade4::dudi.pca(Dataset, scannf = FALSE, nf = 1000)
            Res_rot = Res_pca$c1
            Res_pca_center = as.matrix(bca_output$cent)
            scale = FALSE
            # conduct bgpca
            if(!is.null(ncomp)) {
              if(length(1:ncomp) < ncol(Dataset)) {
                Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = max(ncomp))
              } else {
                Dataset = ade4::bca(Res_pca, as.factor(Response), scannf = FALSE, nf = 1000)
              }
            }
            Res_pca_scores = as.matrix(Dataset$ls) %*% t(Dataset$c1)

          } else if (endsWith(transform, "+GM")) {
            #print(paste("Dataset transformation via Geometric Mean of variables selected"))
            Dataset <- apply(Dataset, 1, function(x) exp(mean(log(x+1e-8))))
            Dataset <- data.frame(GM = Dataset)
            scale = FALSE
          }

          if(isTRUE(scale)) {
            Dataset = data.frame(scale(Dataset))
          }
        } # end of transformation

        if(isTRUE(token)) {
          DT_Output$Response <- Dataset
          if(transform == "pca" || transform == "procdist+pca") {
            DT_Output$Res_dim <- Res_dim
            DT_Output$Res_p = Res_p
            DT_Output$Res_pca <- Res_pca
            DT_Output$Res_rot = Res_rot
            DT_Output$Res_pca_scores = Res_pca_scores
            DT_Output$Res_pca_center = Res_pca_center
          }
        } else {
          DT_Output$Predictor <- Dataset
        }
        return(DT_Output)
      }

    if(!is.null(Response)) {
      Output <- IDT(Dataset = Response,
                     transform = Res_transform,
                     ncomp = Res_ncomp,
                     point_set = Res_point_set,
                     token = TRUE,
                     method = dt_method)

    }

    if(!is.null(Predictor)) {
      Output$Predictor <- IDT(Dataset = Predictor,
                               transform = Pred_transform,
                               ncomp = Pred_ncomp,
                               point_set = Pred_point_set,
                               token = FALSE,
                               method = dt_method)$Predictor

    }

  return(Output)
} # end of data transformation function

