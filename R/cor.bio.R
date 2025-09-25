#' @title Correlation Analysis for sets of Multivariate data
#'
#' @description
#' cor.bio is a function designed to investigate, in parallel, relationships between sets
#' of multivariate data through several permuted correlation tests. This function uses a
#' model-based hypothesis testing structure where a \code{Models} nested list object
#' contains hypothesis tests which usually involve one or more correlation test sets
#' of data. Data can be both subsetted and transformed using the \code{DT} function
#' (see details). The function calculates correlation coefficients and generates
#' correlograms to visualize dataset relationships as well as conducts either
#' Mantel tests (for same data structures) or canonical correlation tests (for
#' different data structures). This function also output correlation summary statistics
#' which can be useful for contexts of integration.
#'
#' @param Models Nested list object containing regression hypotheses to test. Each list element should
#' contain Response and Predictor data, structured as:
#'
#' \itemize{
#'   \item 'Hypothesis Test Models' # Argument Models object name
#'   \itemize{
#'     \item 'Symphyseal shape ~ Posterior Corpus Properties' --- 1st hypothesis model
#'     \itemize{
#'       \item 'Symphyseal shape ~ LM1-M2 Corpus Properties ' --- 1st hypothesis test
#'       \itemize{
#'         \item Symphyseal Landmarks --- Response Data within hypothesis test
#'         \item LM1-M2 Corpus Properties --- Predictor data within hypothesis test
#'         \item dt_parameters --- Optional data transformation parameters
#'       }
#'       \item 'Symphyseal shape ~ LP3-P4 Corpus Properties ' --- 2nd hypothesis test
#'       \itemize{
#'         \item Symphyseal Landmarks --- Response Data within hypothesis test
#'         \item LP3-P4 Corpus Properties --- Predictor data within hypothesis test
#'         \item dt_parameters --- Optional data transformation parameters
#'       }
#'       \item 'Symphyseal shape ~ RP3-P4 Corpus Properties ' --- 3rd hypothesis test
#'       \itemize{
#'         \item Symphyseal Landmarks --- Response Data within hypothesis test
#'         \item RP3-P4 Corpus Properties --- Predictor data within hypothesis test
#'         \item dt_parameters --- Optional data transformation parameters
#'       }
#'       \item 'Symphyseal shape ~ RM1-M2 Corpus Properties ' --- 4th hypothesis test
#'       \itemize{
#'         \item Symphyseal Landmarks --- Response Data within hypothesis test
#'         \item RM1-M2 Corpus Properties --- Predictor data within hypothesis test
#'         \item dt_parameters --- Optional data transformation parameters
#'       }
#'     }
#'     \item 'Symphyseal shape ~ Posterior Corpus Properties' --- 2nd hypothesis model
#'     \itemize{
#'       \item 'Symphyseal shape ~ All Corpus Properties ' --- 1st hypothesis test
#'       \itemize{
#'         \item Symphyseal Landmarks --- Response Data within hypothesis test
#'         \item Posterior Corpus Properties --- Predictor data within hypothesis test
#'         \item dt_parameters --- Optional data transformation parameters
#'       }
#'     }
#'     \item ..... --- 3rd hypothesis model
#'   }
#' }
#'
#' @param point_set Optional nested list specifying variable subsets to analyze. Structure must
#' match Models argument with added subset definitions:
#' \itemize{
#'   \item point_set
#'   \itemize{
#'     \item Model_1 # Argument object name
#'     \itemize{
#'       \item Response               # Response list subsets (e.g., landmarks)
#'       \itemize{
#'         \item Subset_1 = c(1:20)    # Landmarks 1-20
#'         \item Subset_2 = c(21:40)   # Landmarks 21-40
#'         \item Subset_3 = c(41:60)   # Landmarks 41-60
#'       }
#'       \item Predictor              # Predictor subsets (e.g., data table)
#'       \itemize{
#'         \item Subset_1 = c(1:5)     # Variables 1-5
#'         \item Subset_2 = c(6:10)    # Variables 6-10
#'         \item Subset_3 = c(11:15)   # Variables 11-15
#'       }
#'     }
#'   }
#' }
#'
#' @param parallel Logical value indicating whether to use parallel processing
#'
#' @param ... Arguments to pass for the \code{\link{DT}} function.
#'
#' @returns A list containing correlation analysis results for each model:
#' \itemize{
#'   \item{Cor.all.Mantels:   Nested list structure Mantel test results of unsubsetted data including:
#'     \itemize{
#'        \item{r:    Mantel/Canonical correlation statistic}
#'        \item{p:    Statistical significance}
#'        \item{I:    Mean squared correlation (overall correlation strength)}
#'        \item{I_x:    Within-set mean squared correlation for predictor variables}
#'        \item{I_y:    Within-set mean squared correlation for response variables}
#'     }
#'   }
#'   \item{Cor.subset.Mantels:   Nested list structure of Mantel test results of subsetted data including:
#'     \itemize{
#'        \item{r:    Mantel/Canonical correlation statistic}
#'        \item{p:    Statistical significance}
#'        \item{I:    Mean squared correlation (overall correlation strength)}
#'        \item{I_x:    Within-set mean squared correlation for predictor variables}
#'        \item{I_y:    Within-set mean squared correlation for response variables}
#'     }
#'   }
#' }
#'
#' @details
#' cor.bio conducts correlation analyses in parallel for each test within a hypothesis
#' model using both correlation coefficient matrices. For each test within a
#' hypothesis model, the function first performs user-specified data transformations
#' via the \code{DT} function parameters. If no parameters additional DT related
#' parameters are provided through argument (...), the function proceeds with
#' the raw data. When geometric morphometric data is supplied in "p x k x n"
#' (landmarks, dimensions, observations) array format, the function automatically
#' converts the data to a data frame while preserving spatial relationships.
#'
#' The function implements a two-stage correlation analysis approach. First, it calculates
#' standard Pearson correlation coefficients between all variables, generating a complete
#' correlation matrix with associated p-values for significance testing. These correlations
#' are visualized through correlograms generated using the \pkg{corrplot} package, with
#' significant correlations (p < 0.05) distinctly marked. Second, the function conducts
#' either a permuted Mantel test for same structured data or a permuted canonical correlation tests
#' with 1000 permutations to assess matrix-wide correlations. The function outputs
#' a correlation statistic (r) and significance value (p), along with overall (I) and
#' intradata (I_x, I_y) correlation measures.
#'
#' When point_set is provided, the function performs separate correlation analyses
#' sequentially for each defined subset. The function automatically performs two separate
#' tests in subsetted cases: (1) a complete analysis without subsets (Cor.All.Mantels),
#' (2) a series of paired subset analysis (Cor.subset.Mantels). Additionally,
#' the joint data correlation matrices are also output by the function for the
#' same subset comparison (same.subsets) and a trimmed matrix of variables where r >= 0.7
#' (trim.same.subsets) as well as a paired subset comparison matrix (pair.subsets)
#' and a trimmed matrix of variables where r >= 0.7 (trim.dif.subsets).
#'
#' The function generates several output files stored in an organized directory structure.
#' Correlograms are automatically saved as high-resolution PNG files, named according
#' to the model and test structure. For subset analyses, separate plots are generated
#' for each subset or subset combination. Statistical results are stored in the returned
#' list object, organized hierarchically to match the input model structure.
#'
#' For geometric morphometric data, cor.bio first transforms the data. If no data
#' transformation is selected, the data will be transformed as a "n x p" matrix and
#' subject to principal component analysis. This ensures that shape information is
#' appropriately considered in the correlation calculations. All other data types
#' will first be subject to a user-specified data transformation using the \code{DT}
#' function parameters via the argument (...) or as the third object of each hypothesis
#' model test in the Models argument. If no \code{DT} function parameters are
#' provided, no data transformation will occur prior to correlation analysis.
#' The function automatically detects landmark data based on array dimensionality
#' and adjusts its procedures accordingly.
#'
#' @returns A list of class "cor.bio" containing correlation analysis results organized by model:
#'   \itemize{
#'     \item{Cor.all.Mantels:   Mantel test results
#'         \itemize{
#'              \item{r:    Mantel/Canonical correlation statistic}
#'              \item{p:    Statistical significance}
#'              \item{I:    Mean squared correlation (overall correlation strength)}
#'              \item{I_x:    Within-set mean squared correlation for predictor variables}
#'              \item{I_y:    Within-set mean squared correlation for response variables}
#'        }
#'     }
#'     \item{Cor.subset.Mantels:    Optional when point_set is provided
#'         \itemize{
#'              \item{r:    Mantel/Canonical correlation statistic}
#'              \item{p:    Statistical significance}
#'              \item{I:    Mean squared correlation (overall correlation strength)}
#'              \item{I_x:    Within-set mean squared correlation for predictor variables}
#'              \item{I_y:    Within-set mean squared correlation for response variables}
#'        }
#'     }
#'   }
#'
#' @seealso \code{\link{AABA}} \code{\link{Bio.VIP}} \code{\link{DT}}
#'
#' @author Keeling et al., 2025
#'
#' @importFrom stats cor cov cov2cor var as.dist
#' @importFrom Hmisc rcorr
#' @importFrom vegan mantel
#' @importFrom PMA CCA.permute
#' @importFrom corrplot corrplot
#' @importFrom grDevices png dev.off
#' @importFrom dplyr filter
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#'
#' @references
#'
#' Wei, T., & Simko, V. (2021). R package 'corrplot': Visualization of a Correlation
#' Matrix (Version 0.92).
#'
#' Adams, D. C., & Otárola‐Castillo, E. (2013). geomorph: an R package for the
#' collection and analysis of geometric morphometric shape data. Methods in Ecology
#' and Evolution, 4(4), 393-399.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Test for correlations between mandibular symphysis shape and cross-sectional properties
#' # at different cross-sections
#'
#' # Create model structure
#' Models <- list(
#'   'Symphyseal_Shape ~ Posterior_Corpus_Biomechanics' = list(
#'     'Symphysis_LM1M2' = list(
#'       'Response' = Corpus_Land[241:360,,],  # Symphysis landmarks
#'       'Predictor' = Biomechanics %>%
#'         filter(Region == "LM1M2") %>%
#'         select(5:93)  # Biomechanical variables
#'     ),
#'     'Symphysis_LP3P4' = list(
#'       'Response' = Corpus_Land[241:360,,],
#'       'Predictor' = Biomechanics %>%
#'         filter(Region == "LP3P4") %>%
#'         select(5:93)
#'     ),
#'     'Symphysis_RP3P4' = list(
#'       'Response' = Corpus_Land[241:360,,],
#'       'Predictor' = Biomechanics %>%
#'         filter(Region == "RP3P4") %>%
#'         select(5:93)
#'     ),
#'     'Symphysis_LP3P4' = list(
#'       'Response' = Corpus_Land[241:360,,],
#'       'Predictor' = Biomechanics %>%
#'         filter(Region == "RM1M2") %>%
#'         select(5:93)
#'     ),
#'   )
#' )
#'
#' # Define variable subsets for testing which cross-sectional properties are important in symphyseal shape
#' point_set <- list('Symphyseal Shape ~ Posterior Corpus Properties' =
#'                      list(
#'                        'Symphysis_Shape ~ LM1M2' = list(
#'                          'Symphysis_Shape' = list(
#'                               'Subset_1' = c(1:120), # symphysis has 120 total landmarks
#'                               'Subset_2' = c(1:120),
#'                               'Subset_3' = c(1:120),
#'                               'Subset_4' = c(1:120),
#'                               'Subset_5' = c(1:120)
#'                               ),
#'
#'                          'LM1M2' = list(
#'                               'Dimensions' = c(5,7:12,14:19),
#'                               'Orientation' = c(13, 28),
#'                               'Cortical Thickness' = c(6, 20:22, 34:93),
#'                               'Bending Resistance' = c(23:27),
#'                               'Breaking Strength' = c(29:33)
#'                               )
#'                          ),
#'                        'Symphysis_Shape ~ LP3P4' = list(
#'                          'Symphysis_Shape' = list(
#'                               'Subset_1' = c(1:120), # symphysis has 120 total landmarks
#'                               'Subset_2' = c(1:120),
#'                               'Subset_3' = c(1:120),
#'                               'Subset_4' = c(1:120),
#'                               'Subset_5' = c(1:120)
#'                               ),
#'
#'                          'LP3P4' = list(
#'                               'Dimensions' = c(5,7:12,14:19),
#'                               'Orientation' = c(13, 28),
#'                               'Cortical Thickness' = c(6, 20:22, 34:93),
#'                               'Bending Resistance' = c(23:27),
#'                               'Breaking Strength' = c(29:33)
#'                               )
#'                          ),
#'                         'Symphysis_Shape ~ RP3P4' = list(
#'                          'Symphysis_Shape' = list(
#'                               'Subset_1' = c(1:120), # symphysis has 120 total landmarks
#'                               'Subset_2' = c(1:120),
#'                               'Subset_3' = c(1:120),
#'                               'Subset_4' = c(1:120),
#'                               'Subset_5' = c(1:120)
#'                               ),
#'
#'                          'RP3P4' = list(
#'                               'Dimensions' = c(5,7:12,14:19),
#'                               'Orientation' = c(13, 28),
#'                               'Cortical Thickness' = c(6, 20:22, 34:93),
#'                               'Bending Resistance' = c(23:27),
#'                               'Breaking Strength' = c(29:33)
#'                               )
#'                          ),
#'                         'Symphysis_Shape ~ RM1M2' = list(
#'                          'Symphysis_Shape' = list(
#'                               'Subset_1' = c(1:120), # symphysis has 120 total landmarks
#'                               'Subset_2' = c(1:120),
#'                               'Subset_3' = c(1:120),
#'                               'Subset_4' = c(1:120),
#'                               'Subset_5' = c(1:120)
#'                               ),
#'
#'                          'RM1M2' = list(
#'                               'Dimensions' = c(5,7:12,14:19),
#'                               'Orientation' = c(13, 28),
#'                               'Cortical Thickness' = c(6, 20:22, 34:93),
#'                               'Bending Resistance' = c(23:27),
#'                               'Breaking Strength' = c(29:33)
#'                               )
#'                          )
#'                       )
#'                   )
#'
#' # Run correlation analysis
#' results <- cor.bio(
#'   Models = Models, # add in hypothesis models
#'   point_set = point_set, # add in subsets
#'   paired_subsets = FALSE # since the Response data stays the same, its unnecessary to do cross comparisons
#' ) # keep defaults
#'
#' # Examine results
#' # View the trimmed correlation matrix based off of the highest correlations
#' print(results$Symphysis_Correlations$Symphysis_LM1M2$Trim.subset)
#'
#' # Check the overall Mantel test results
#' print(results$Symphysis_Correlations$Symphysis_LM1M2$Cor.all.Mantels)
#'
#' # Check the first subset test (Corpus Dimensions) results for 'Symphysis_Shape ~ LM1M2'
#' print(results[[1]][[1]][[1]])
#'
#' # Reading the result structure
#' # results     This is the output
#' # [[1]]       This is the first hypothesis model 'Symphyseal Shape ~ Posterior Corpus Properties'
#' # [[1]]       This is the first test in the first hypothesis model 'Symphysis_Shape ~ LM1M2'
#' # [[1]]       (Optional if point_set is provided) First subset test for Symphysis_Shape ~ LM1M2'
#'
#' # Check the second subset test (Corpus Orientation) results now
#' print(results[[1]][[1]][[2]])
#'
#' # Check the first subset test (Corpus dimensions) results for 'Symphysis_Shape ~ RM1M2'
#' print(results[[1]][[4]][[1]])
#' }

cor.bio <- function(Models,
                    point_set = NULL,
                    paired_subsets = FALSE,
                    parallel = TRUE,
                    core_choice = "detect",
                    ...) {

  set.seed(1)
  #___________________________________________________________________________________________________________________________
  Dir = getwd()
  if(!dir.exists("Correlations_Bio")) {
    dir.create("Correlations_Bio")
  }
  ct.dir = "Correlations_Bio"
  path_cor = paste(Dir, "/", ct.dir, sep = "")
  setwd(path_cor)
  #___________________________________________________________________________________________________________________________
  cor.bio = vector("list", length(Models))

  subset = if(is.null(point_set)) FALSE else TRUE
  #___________________________________________________________________________________________________________________________
  for(m in seq_along(Models)) {

    #cor.all = vector("list", length(Models[[m]]))
    #cor.subset = vector("list", length(Models[[m]]))
    #trim.subset = vector("list", length(Models[[m]]))
    Cor.subset.Mantels = vector("list", length(Models[[m]]))
    Cor.all.Mantels = vector("list", length(Models[[m]]))
    cam = list()
    csm = list()
    cdm = list()

    if(isTRUE(subset)) {
      for(i in seq_along(Models[[m]])) {
        #cor.subset[[i]] <- vector("list", length(point_set[[m]][[i]][[2]]))
        #trim.subset[[i]] <- vector("list", length(point_set[[m]][[i]][[2]]))
        Cor.subset.Mantels[[i]] <- vector("list", length(point_set[[m]][[i]][[2]]))
        if(length(point_set[[m]][[i]][[2]]) == 1) {
          p_subset <- 1
          #cor.subset[[i]] <- vector("list", p_subset)
          #trim.subset[[i]] <- vector("list", p_subset)
          Cor.subset.Mantels[[i]] <- vector("list", p_subset)
        } else {
          if(isTRUE(paired_subsets)) {
            p_subset = expand.grid(seq_along(point_set[[m]][[i]][[1]]), seq_along(point_set[[m]][[i]][[2]]))
          } else {
            p_subset = data.frame(Response = rep(1:length(point_set[[m]][[i]][[1]])),
                                  Predictor = rep(1:length(point_set[[m]][[i]][[2]]))
            )
          }
         # cor.subset[[i]] <- vector("list", nrow(p_subset))
         # trim.subset[[i]] <- vector("list", nrow(p_subset))
          Cor.subset.Mantels[[i]] <- vector("list", nrow(p_subset))
        }
      }
    } # end of if subsets are true

    if(!dir.exists(paste(names(Models)[m]))) {
      dir.create(paste(names(Models)[m]))
    }
    model.dir = paste(names(Models)[m])
    if(!dir.exists(paste(model.dir))) {
      dir.create(paste(model.dir))
    }
    path = file.path(Dir, ct.dir, model.dir)
    #________________________________________________________________________________________________________________________________

    if(core_choice == "detect") {
      core_choice = parallel::detectCores()- 2
    }
    if(core_choice == "high") {
      core_choice = 28
    }
    if(core_choice == "medium-high") {
      core_choice = 20
    }
    if(core_choice == "medium") {
      core_choice = 14
    }
    if(core_choice == "medium-low") {
      core_choice = 10
    }
    if(core_choice == "low") {
      core_choice = 8
    }
    core_free = (detectCores()- 2) - core_choice

    #_______________________________________________________________________________________________________________________________

    args = list(...)

    dt_params = names(formals(DT))

    dt_parameters <- list(Response = NULL,
                          Predictor = NULL,
                          Res_transform = "NULL",
                          Pred_transform = "NULL",
                          Res_ncomp = NULL,
                          Pred_ncomp = NULL,
                          point_set_res = NULL,
                          point_set_pred = NULL)

    dt_parameters <- modifyList(dt_parameters, args) %>% .[intersect(names(.), dt_params)]

    #________________________________________________________________________________________________________________________________


    if (isTRUE(parallel)) {
      doParallel::registerDoParallel(parallel::detectCores() - 2)
    } else {
      registerDoSEQ()
    }
    parallels <- foreach(i = 1:length(Models[[m]]), .combine = "rbind", .packages =
                           c("corrplot", "factoextra", "FactoMineR", "geomorph","grDevices", "graphics", "caret", "vegan", "stats", "PMA", "purrr",
                             "Hmisc", "tidyverse")) %dopar% {

     dt_parameters = tryCatch({modifyList(dt_parameters,Models[[m]][[i]][[3]])}, error = function(e){dt_parameters})

     dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                    Predictor = Models[[m]][[i]][[2]],
                                                    Res_transform = if(length(dim(Models[[m]][[i]][[1]])) == 3 & dt_parameters$Res_transform == "NULL") "2D" else dt_parameters$Res_transform,
                                                    Pred_transform = if(length(dim(Models[[m]][[i]][[2]])) == 3 & dt_parameters$Pred_transform == "NULL") "2D" else dt_parameters$Pred_transform)
     ) %>% .[intersect(names(.), dt_params)]

       Response = data.frame(do.call(DT, dt_parameters)$Response)
       Predictor = data.frame(do.call(DT, dt_parameters)$Predictor)

       if(length(dim(Response)) == 3) {
         y = as.matrix.data.frame(geomorph::two.d.array(Response))
         I_y = mean(cor(y)^2)
         y = as.matrix.data.frame(geomorph::gm.prcomp(y)$x %>% if(ncol(.) < 10) . else .[,1:10])
       } else {
         y = as.matrix.data.frame(Response)
         I_y = mean(cor(y)^2)
       }
       if(length(dim(Predictor)) == 3) {
         x = as.matrix.data.frame(geomorph::two.d.array(Predictor))
         I_x = mean(cor(x)^2)
         x = as.matrix.data.frame(geomorph::gm.prcomp(x)$x %>% if(ncol(.) < 10) . else .[,1:10])
       } else {
         x = as.matrix.data.frame(Predictor)
         I_x = mean(cor(x)^2)
       }
       corr_results <- Hmisc::rcorr(x,y)
       cor_matrix <- corr_results$r
       I = mean(cor_matrix^2, na.rm = TRUE)
       p_matrix <- corr_results$P
       p_matrix[is.na(p_matrix)] <- 0.99
       corr_mask <- p_matrix <= 0.05
       #cor.all <- corr_results
       if(ncol(x) == ncol(y)) {
         cam <- vegan::mantel(cor(x), cor(y))
         Cor.all.Mantels <- data.frame(r = cam$statistic, p = cam$signif, I = I, I_x = I_x, I_y = I_y)
       } else {
           cam = PMA::CCA.permute(x, y, "standard", "standard")
           cam = PMA::CCA.permute(x, y, "standard", "standard",cam$bestpenaltyx, cam$bestpenaltyz)
           cam = data.frame(r = mean(cam$corperms, na.rm = TRUE), p = mean(cam$pvals, na.rm = TRUE),I = I, I_x = I_x, I_y = I_y)
           Cor.all.Mantels <- cam
       }

       if(ncol(x) > 100){
         nsize = 0.1
         tsize = 0.25
         lsize = 0.3
       }else if(ncol(x) > 60) {
         nsize = 0.2
         tsize = 0.3
         lsize = 0.4
       }else if(ncol(x) > 50) {
         nsize = 0.25
         tsize = 0.35
         lsize = 0.5
       } else if(ncol(x) < 50) {
         nsize = 0.4
         tsize = 0.4
         lsize = 0.6
       }

       # Call an image then save it to the created directory
       grDevices::png(paste(path,"/", names(Models[[m]][[i]])[1], "-", names(Models[[m]][[i]])[2], ".png", sep = ""),
                      res = 600, width = 1920, height = 1080)
       # Generate the correlation plot using corrplot


       if(identical(x,y)) {
         corrplot::corrplot(cor_matrix[-c(1:ncol(x)), c(1:(ncol(x)))], method = "circle", type = "upper",
                            order="original", tl.col="black", tl.srt=45, tl.cex = tsize, cl.cex = lsize, number.cex = nsize,
                            p.mat = p_matrix[-c(1:ncol(x)), c(1:(ncol(x)))], sig.level = 0.05, insig = "blank",
                            addCoef.col = "black")

       } else {

         corrplot::corrplot(cor_matrix[-c(1:ncol(x)), c(1:(ncol(x)))], method = "circle", type = "full",
                            order="original", tl.col="black", tl.srt=45, tl.cex = tsize, cl.cex = lsize, number.cex = nsize,
                            p.mat = p_matrix[-c(1:ncol(x)), c(1:(ncol(x)))], sig.level = 0.05, insig = "blank",
                            addCoef.col = "black")
       }

       dev.off()

       if(isTRUE(subset)) {
         if(isTRUE(paired_subsets)) {
           p_subset = expand.grid(seq_along(point_set[[m]][[i]][[1]]), seq_along(point_set[[m]][[i]][[2]]))
         } else {
           p_subset = data.frame(Response = rep(1:length(point_set[[m]][[i]][[1]])),
                                 Predictor = rep(1:length(point_set[[m]][[i]][[2]]))
           )
         }



         for(j in 1:nrow(p_subset)) {

           dt_parameters = tryCatch({modifyList(dt_parameters,Models[[m]][[i]][[3]]) %>% .[intersect(names(.), dt_params)]},
                    error = function(e){dt_parameters})

           dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                          Predictor = Models[[m]][[i]][[2]],
                                                          Res_point_set = point_set[[m]][[i]][[1]][[p_subset[j,1]]],
                                                          Pred_point_set = point_set[[m]][[i]][[2]][[p_subset[j,2]]]
                                                          )
                                      ) %>% .[intersect(names(.), dt_params)]

           Response = data.frame(do.call(DT, dt_parameters)$Response)
           Predictor = data.frame(do.call(DT, dt_parameters)$Predictor)

           if(length(dim(Response)) == 3) {
             y = as.matrix.data.frame(two.d.array(Response))
             I_y = mean(cor(y)^2)
             y = as.matrix.data.frame(gm.prcomp(y)$x %>% if(ncol(.) < 10) . else .[,1:10])
           } else {
             y = as.matrix.data.frame(Response)
             I_y = mean(cor(y)^2)
           }
           if(length(dim(Predictor)) == 3) {
             x = as.matrix.data.frame(two.d.array(Predictor))
             I_x = mean(cor(x)^2)
             x = as.matrix.data.frame(gm.prcomp(x)$x %>% if(ncol(.) < 10) . else .[,1:10])
           } else {
             x = as.matrix.data.frame(Predictor)
             I_x = mean(cor(x)^2)
           }

           corr_results <- Hmisc::rcorr(x,y)
           cor_matrix <- corr_results$r
           I = mean(cor_matrix^2, na.rm = TRUE)
           p_matrix <- corr_results$P
           p_matrix[is.na(p_matrix)] <- 0.99
           corr_mask <- p_matrix <= 0.05
           #cor.subset[[j]] <- corr_results
           #remove_var <- findCorrelation(cor_matrix, cutoff = 0.7, names = TRUE)
           #trim.subset[[j]] <- cor_matrix[(rownames(cor_matrix) %in% remove_var), (colnames(cor_matrix) %in% remove_var)]

           if(identical(x,y)) {
             cdm <- mantel(cor(x), cor(y))
             Cor.subset.Mantels[[i]][[j]] <- data.frame(r = cdm$statistic, p = cdm$signif, I = I, I_x = I_x, I_y = I_y)
           } else {
               cam = PMA::CCA.permute(x,y, "standard", "standard")
               cam = PMA::CCA.permute(x,y, "standard", "standard",cam$bestpenaltyx, cam$bestpenaltyz)
               cam = data.frame(r = mean(cam$corperms, na.rm = TRUE), p = mean(cam$pvals, na.rm = TRUE),I = I, I_x = I_x, I_y = I_y)
               Cor.subset.Mantels[[i]][[j]] <- cam
           }
           cdm = list()

           if(ncol(x) > 100){
             nsize = 0.1
             tsize = 0.25
             lsize = 0.3
           }else if(ncol(x) > 60) {
             nsize = 0.2
             tsize = 0.3
             lsize = 0.4
           }else if(ncol(x) > 50) {
             nsize = 0.25
             tsize = 0.35
             lsize = 0.5
           } else if(ncol(x) < 50) {
             nsize = 0.4
             tsize = 0.4
             lsize = 0.6
           }

           # Call an image then save it to the created directory
           grDevices::png(paste(path,"/", names(Models[[m]])[i], "_", "Subsets", "_",
                                names(point_set[[m]][[i]][[1]])[[p_subset[j,1]]], "_",
                                names(point_set[[m]][[i]][[2]])[[p_subset[j,2]]], ".png", sep = ""),
                          res = 600, width = 1920, height = 1080)
           # Generate the correlation plot using corrplot

           if(identical(x,y)) {
             corrplot::corrplot(cor_matrix[-c(1:ncol(x)), c(1:(ncol(x)))], method = "circle", type = "upper",
                                order="original", tl.col="black", tl.srt=45, tl.cex = tsize, cl.cex = lsize, number.cex = nsize,
                                p.mat = p_matrix[-c(1:ncol(x)), c(1:(ncol(x)))], sig.level = 0.05, insig = "blank",
                                addCoef.col = "black")

           } else {

             corrplot::corrplot(cor_matrix[-c(1:ncol(x)), c(1:(ncol(x)))], method = "circle", type = "full",
                                order="original", tl.col="black", tl.srt=45, tl.cex = tsize, cl.cex = lsize, number.cex = nsize,
                                p.mat = p_matrix[-c(1:ncol(x)), c(1:(ncol(x)))], sig.level = 0.05, insig = "blank",
                                addCoef.col = "black")
           }
           dev.off()
         } # end of j for loop
       } # end of subset

       if(isTRUE(subset)) {
         return(list(#cor.all = cor.all,
                     Cor.all.Mantels = Cor.all.Mantels,
                     #cor.subset = cor.subset,
                     Cor.subset.Mantels = Cor.subset.Mantels
                     #trim.subset = trim.subset)
                ))
       } else {
         return(list(#cor.all = cor.all,
                     Cor.all.Mantels = Cor.all.Mantels))
       }
     } # end of i for loop

    if(isTRUE(parallel)) {
      closeAllConnections()
    }

    if(length(Models[[m]]) == 1) {
      list2env(parallels, envir = environment())
    } else {

      if(isTRUE(subset)) {
        var_boot = c("Cor.all.Mantels", "cor.subset.Mantels")
        booted.list = vector("list", length(var_boot))
        for(v in 1:length(var_boot)) {
          boot_temp = vector("list", nrow(parallels))
          for(b in 1:nrow(parallels)) {
            boot_temp[[b]] <- unlist(parallels[b,v], recursive = FALSE)
            boot_temp[[b]] <- Filter(function(x) !is.null(x), boot_temp[[b]])
            boot_temp[[b]] <- unlist(boot_temp[[b]], recursive = FALSE)
            boot_temp[[b]] <- Filter(function(x) !is.null(x), boot_temp[[b]])
          }
          booted.list[[v]] <- boot_temp
        }
        names(booted.list) <- var_boot
        list2env(booted.list, envir = environment())
      } else {
        var_boot = c("Cor.all.Mantels")
        booted.list = vector("list", length(var_boot))
        boot_temp = vector("list", nrow(parallels))
        for(b in 1:nrow(parallels)) {
          boot_temp[[b]] <- unlist(parallels[b,v], recursive = FALSE)
          boot_temp[[b]] <- Filter(function(x) !is.null(x), boot_temp[[b]])
          boot_temp[[b]] <- unlist(boot_temp[[b]], recursive = FALSE)

        }
        booted.list[[v]] <- boot_temp
      }
      names(booted.list) <- var_boot
      list2env(booted.list, envir = environment())
    }

    #___________________________________________________________________________________________________________________________

    #cor.bio[[m]]$Correlations <- cor.all
    cor.bio[[m]]$Cor.all.Mantels <- Cor.all.Mantels

    if(isTRUE(subset)) {
      cor.bio[[m]]$Cor.subset.Mantels <- cor.subset.Mantels
      #cor.bio[[m]]$Cor.subsets <- cor.subset
      #cor.bio[[m]]$Trim.subset <- trim.subset
    }

    names(cor.bio) <- names(Models)

    if(isTRUE(subset)) {
      cor.bio[[m]] <- lapply(cor.bio[[m]], function(comparisons) {
        names(comparisons) <- names(Models[[m]])
        comparisons
        lapply(comparisons, function(subsets) {
          names(subsets) <- paste(rep("Subset"), rep(1:length(subsets)))
          subsets
        })
      })
    } else {
      cor.bio[[m]] <- lapply(cor.bio[[m]], function(comparisons) {
        names(comparisons) <- names(Models[[m]])
        comparisons
      })
    }
  } # end of m for loop
  #___________________________________________________________________________________________________________________________
  setwd(Dir)
  class(cor.bio) <- "cor.bio"
  return(cor.bio)
} #end of Correlation function for Biological variables
