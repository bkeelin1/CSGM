#' @title Procrustes Shape Linear Regression (pslr.bio)
#'
#' @description This function conducts, in parallel, a Procrustes linear regression
#' on multiple models between landmark arrays (Response variable) and predictor
#' variables (univariate or multivariate). The function was designed to aid in
#' studies where multiple regression hypotheses need to be tested on different
#' and/or the same predictor variables. This function permits conducting regressions
#' on subsets or combinations of subsetted landmarks and/or variables.
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
#' @param paired_subsets a logical value indicating whether or not the function
#' should iterate through regressions for every possible combination of subsets.
#' If argument paired_subsets = TRUE and argument point_set is not NULL, then
#' the function will calculate all of the possible subset combinations between
#' the response and predictor variables and conduct separate regressions for
#' each subsetted combination.
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
#' Each subset should be a numeric vector specifying which landmarks/variables to include.
#' Subsets can be of different sizes and overlap. When paired_subsets=TRUE, all
#' possible combinations of response and predictor subsets will be analyzed.
#'
#' @param all_vips Optional nested list object containing variables of interest to trim the predictor dataset.
#' The structure should follow a similar structure to the  point_set argument where
#' for each m hypothesis model, there are i hypothesis tests. Within each hypothesis
#' test should be a vector called "vip" of variables to include in the analysis.
#' \itemize{
#'   \item all_vips
#'   \itemize{
#'     \item Model_1 (1st hypothesis model)
#'     \itemize{
#'       \item Test_1 = c("var1", "var2", "var3", ...) # Variables for 1st test
#'       \item Test_2 = c("var1", "var4", "var8", ...) # Variables for 2nd test
#'       \item ... # More tests
#'     }
#'     \item Model_2 (2nd hypothesis model)
#'     \itemize{
#'       \item Test_1 = c("var2", "var5", "var9", ...) # Variables for 1st test
#'       \item ... # More tests
#'     }
#'     \item ... # More models
#'   }
#' }
#'
#' @param subset_vips Optional nested list object containing variables for each subset within hypothesis tests.
#' The structure should follow a similar structure to the  point_set argument where
#' for each m hypothesis model, there are i hypothesis tests. Within each hypothesis
#' test should be j list objects of character or numeric vectors called "vip" to
#' include in the analysis for each subset.
#' \itemize{
#'   \item subset_vips
#'   \itemize{
#'     \item Model_1 (1st hypothesis model)
#'     \itemize{
#'       \item Test_1 (1st hypothesis test)
#'       \itemize{
#'         \item Subset_1 = c("var1", "var2", ...) # Variables for subset 1
#'         \item Subset_2 = c("var3", "var5", ...) # Variables for subset 2
#'         \item ... # More subsets
#'       }
#'       \item Test_2 (2nd hypothesis test)
#'       \itemize{
#'         \item Subset_1 = c("var2", "var4", ...) # Variables for subset 1
#'         \item ... # More subsets
#'       }
#'       \item ... # More tests
#'     }
#'     \item Model_2 (2nd hypothesis model)
#'     \itemize{
#'       \item ... # Similar structure
#'     }
#'     \item ... # More models
#'   }
#' }
#'
#' @param key Character string. Identifier for output files and directories. Default is NULL.
#'
#' @details The pslr.bio function performs Procrustes shape linear regression analysis by:
#'
#' 1. Creating working directories for results storage
#' 2. Processing multiple regression models in parallel
#' 3. Handling data transformations (PCA, Procrustes distances, etc.)
#' 4. Supporting subset analysis of landmarks and variables
#' 5. Implementing paired comparisons when requested
#' 6. Managing error handling and alternative approaches
#' 7. Writing results to Excel in summarized format
#'
#' The function uses geomorph::procD.lm() for the regression analysis with:
#' - Parallel processing enabled
#' - 1000 iterations
#' - Randomized residual permutation procedure (RRPP)
#'
#' @returns Returns a nested list containing:
#' \itemize{
#'   \item{Model_1:   First Hypothesis Model
#'     \itemize{ First Test within first Hypothesis Model
#'       \item{ANOVA tables         Contains SS, MS, Rsq, F, Z, and P-values}
#'       \item{Effect size          R-squared values for each model}
#'       \item{Model coefficients   When applicable}
#'     }
#'   }
#' }
#'
#' If write.results=TRUE, also generates Excel files containing:
#' - Separate sheets for each regression analysis
#' - Summary statistics and test results
#' - Organized by model structure
#'
#' Directory structure created:
#'  - pslr.bio_Results_(key)/Model_(model_name).xlsx
#'
#' @seealso \code{procD.lm} \code{DT} \code{restructure_results} \code{AABA}
#'
#' @author Keeling et al., 2025
#'
#' @importFrom geomorph procD.lm arrayspecs two.d.array
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom writexl write_xlsx
#' @importFrom dplyr filter mutate %>% intersect
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Scenario Question: How much of symphyseal shape variance can be explained
#' #   by a linear regression of the cross-sectional properties of the posterior corpus?
#'
#' # Posterior corpus in this example refers to the following cross-sections:
#'
#' # LM1M2 - interdental junction of the left side first and second molars.
#' # LP3P4 - interdental junction of the left side third and fourth premolars.
#' # RP3P4 - interdental junction of the right side third and fourth premolars.
#' # RM1M2 - interdental junction of the right side first and second molars.
#'
#' # Import data sample of mandibular corpus: Corpus_Land
#' # Import data sample of mandibular corpus data table: Biomechanics
#'
#' data(Corpus_Land)
#' data(Biomechanics)
#'
#'
#' # Create a data analysis model
#'
#' Models <- list('Regression_Model' =
#'                   list(
#'                        'Symphysis_Shape ~ LM1M2' = list(
#'                          'Symphysis_Shape' = Corpus_Land[241:360,,], # symphyseal landmarks
#'                          'LM1M2' = Biomechanics %>% filter(Region = "LM1M2") # LM1M2 filtered data
#'                          ),
#'                        'Symphysis_Shape ~ LP3P4' = list(
#'                          'Symphysis_Shape' = Corpus_Land[241:360,,],
#'                          'LM1M2' = Biomechanics %>% filter(Region = "LM1M2")
#'                          ),
#'                        'Symphysis_Shape ~ RP3P4' = list(
#'                          'Symphysis_Shape' = Corpus_Land[241:360,,],
#'                          'LM1M2' = Biomechanics %>% filter(Region = "LM1M2")
#'                          ),
#'                        'Symphysis_Shape ~ RM1M2' = list(
#'                          'Symphysis_Shape' = Corpus_Land[241:360,,],
#'                          'LM1M2' = Biomechanics %>% filter(Region = "LM1M2")
#'                          ),
#'                        )
#'                )
#'
#'
#' # Let's subset the filtered data tables by cross-sectional property types:
#'
#'   Dimensions: Variables 1:3, 5:8, 10:15
#'   Orientation: Variables 9, 24
#'   Cortical Thickness: Variables 4, 16:18, 30:60
#'   Bending Resistance:Variables 19:23
#'   Breaking Strength: Variables 25:29
#'
#'
#' # Create the point_set argument model for subsetted variables.
#' # IMPORTANT: If no subsets are desired for either the response or predictor data
#' # all variables must be inputted in point_set subsets. For example:
#'
#' point_set <- list('Regression' =
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
#'                               'Dimensions' = c(1:3, 5:8, 10:15),
#'                               'Orientation' = c(9, 24),
#'                               'Cortical Thickness' = c(4, 16:18, 30:60),
#'                               'Bending Resistance' = c(19:23),
#'                               'Breaking Strength' = c(25:29)
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
#'                               'Dimensions' = c(1:3, 5:8, 10:15),
#'                               'Orientation' = c(9, 24),
#'                               'Cortical Thickness' = c(4, 16:18, 30:60),
#'                               'Bending Resistance' = c(19:23),
#'                               'Breaking Strength' = c(25:29)
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
#'                               'Dimensions' = c(1:3, 5:8, 10:15),
#'                               'Orientation' = c(9, 24),
#'                               'Cortical Thickness' = c(4, 16:18, 30:60),
#'                               'Bending Resistance' = c(19:23),
#'                               'Breaking Strength' = c(25:29)
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
#'                               'Dimensions' = c(1:3, 5:8, 10:15),
#'                               'Orientation' = c(9, 24),
#'                               'Cortical Thickness' = c(4, 16:18, 30:60),
#'                               'Bending Resistance' = c(19:23),
#'                               'Breaking Strength' = c(25:29)
#'                               )
#'                          )
#'                       )
#'                   )
#'
#' # Run the pslr.bio function to conduct parallel regressions.
#'
#' output <- pslr.bio(Models, # the regression model created above
#'                subset = TRUE, # there are subsets in this regression test
#'                paired_subsets = FALSE, # Since response variable contains no unique subsets, this is not necessary
#'                point_set, # this is the subset model that was created above
#'                key = "Corpus", # Name for the output file folder.
#'                ) # keep everything else as default
#'
#'
#'
#'
#' # View the results of each regression test separately in the output object
#'
#' # View the Symphysis ~ LM1M2 results
#'
#' output[[1]][[1]][[1]] # structure is output, model, regression, subset
#' # brackets indicate first model "Regression" and the first regression test
#' # "Symphysis_Shape~LM1M2" and the first variable subset "Dimensions"
#'
#'
#' # view the results via the excel file listed in the folder "pslr.bio Results-Corpus":
#' # File name is "pslr.bio_Model_Regression.xlsx"
#'
#' pslr.bio_Results = read_xlsx("pslr.bio_Model_Regression.xlsx")
#'
#' view(pslr.bio_Results)
#' }

pslr.bio <- function(Models,
                 paired_subsets = FALSE,
                 point_set = NULL,
                 all_vips = NULL,
                 subset_vips = NULL,
                 key = NULL,
                 ...) {

  set.seed(1)

  Dir = getwd()
  if(!is.null(key)) {
    path = paste("pslr.bio_Results-", key)
    path = file.path(Dir, path)
    if(!dir.exists(path)) {
      dir.create(path)
    }

  } else {
    path = file.path(Dir, "pslr.bio_Results")
    if(!dir.exists("pslr.bio_Results")) {
      dir.create("pslr.bio_Results")
    }
  }
  setwd(path)

  args = list(...)
  subset = if(!is.null(point_set)) TRUE else FALSE

  dt_params = names(formals(DT))

  dt_parameters = args %>% .[intersect(names(.), dt_params)]
  if(length(dt_parameters) == 0){
    dt_parameters = list(dt_method = "procdist")
  }

  Output = vector("list", length(Models))
  names(Output) = names(Models)

  #___________________________________________________________________________________________________________

  for(m in seq_along(Models)) {

    pslr.results.all = vector("list", length(Models[[m]]))
    names(pslr.results.all) <- names(Models[[m]])
    if(isTRUE(subset)) {
      pslr.results = vector("list", length(Models[[m]]))
      names(pslr.results) <- names(Models[[m]])
      for(i in seq_along(point_set[[m]])) {
        if(isTRUE(paired_subsets)) {
          p_subset = expand.grid(seq_along(point_set[[m]][[i]][[1]]), seq_along(point_set[[m]][[i]][[2]]))
          pslr.results[[i]] <- vector("list", nrow(p_subset))
          names(pslr.results[[i]]) <- paste(rep("Subset",nrow(p_subset)),rep(1:nrow(p_subset)), sep = "_")
        } else {
          p_subset = data.frame(Response = rep(1:length(point_set[[m]][[i]][[1]])),
                                Predictor = rep(1:length(point_set[[m]][[i]][[2]])))
          pslr.results[[i]] <- vector("list", nrow(p_subset))
          names(pslr.results[[i]]) <- paste(rep("Subset",nrow(p_subset)),rep(1:nrow(p_subset)), sep = "_")
        }
      }
    }

    for(i in seq_along(Models[[m]])) {

      if(isTRUE(subset)) {
        if(isTRUE(paired_subsets)) {
          p_subset = expand.grid(seq_along(point_set[[m]][[i]][[1]]), seq_along(point_set[[m]][[i]][[2]]))
        } else {
          p_subset = data.frame(Response = rep(1:length(point_set[[m]][[i]][[1]])),
                                Predictor = rep(1:length(point_set[[m]][[i]][[2]])))
        }

        for(j in 1:nrow(p_subset)) {

          dt_parameters = tryCatch({modifyList(dt_parameters,Models[[m]][[i]][[3]]) %>% .[intersect(names(.),dt_params)]},
                                   error = function(e){dt_parameters})

          dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                         Predictor = Models[[m]][[i]][[2]],
                                                         Res_point_set = point_set[[m]][[i]][[1]][[p_subset[j,1]]],
                                                         Pred_point_set = point_set[[m]][[i]][[2]][[p_subset[j,2]]]))

          Response <<- do.call(DT, dt_parameters)$Response
          if(length(dim(Response)) != 3){
            Response <<- as.matrix.data.frame(Response)
          }

          Predictor <<- do.call(DT, dt_parameters)$Predictor
          if(length(dim(Predictor)) != 3){
            Predictor <<- as.matrix.data.frame(Predictor)
          }

          if(!is.null(subset_vips)) {
            Predictor <<- as.matrix.data.frame(if(length(dim(Predictor)) == 3) Predictor[subset_vips[[m]][[i]][[j]],,] else Predictor[,subset_vips[[m]][[i]][[j]]])
          }

          pslr.results[[i]][[j]]$summary <- data.frame(suppressWarnings(tryCatch({base::summary(geomorph::procD.lm(Response ~ Predictor,
                                                                                                        Parallel = TRUE, iter = 1000, RRPP = TRUE))$table},
                                                                                      error = function(e) {
                                                                                        Predictor <<- as.matrix.data.frame(geomorph::gm.prcomp(Predictor)$x)

                                                                                        base::summary(geomorph::procD.lm(Response ~ Predictor, Parallel = TRUE, iter = 1000,
                                                                                                                             RRPP = TRUE))$table
                                                                                      })) %>% cbind(.[1,], .[2,1:4], .[3,1:2]) %>% .[1,])


        }
      } # end of subsets

      dt_parameters = tryCatch({modifyList(dt_parameters,Models[[m]][[i]][[3]]) %>% .[intersect(names(.),dt_params)]},
               error = function(e){dt_paramters})

      dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                     Predictor = Models[[m]][[i]][[2]],
                                                     Res_point_set = NULL,
                                                     Pred_point_set = NULL))

      Response <<- do.call(DT, dt_parameters)$Response
      if(length(dim(Response)) != 3){
        Response <<- as.matrix.data.frame(Response)
      }

      Predictor <<- do.call(DT, dt_parameters)$Predictor
      if(length(dim(Predictor)) != 3){
        Predictor <<- as.matrix.data.frame(Predictor)
      }
      if(!is.null(all_vips)) {
        Predictor <<- as.matrix.data.frame(if(length(dim(Predictor)) == 3) Predictor[all_vips[[m]][[i]],,] else Predictor[,all_vips[[m]][[i]]])
      }

      pslr.results.all[[i]]$summary <- data.frame(suppressWarnings(tryCatch({base::summary(geomorph::procD.lm(Response ~ Predictor,
                                                                                                   Parallel = TRUE, iter = 1000, RRPP = TRUE))$table},
                                                                                 error = function(e) {
                                                                                   Predictor <<- as.matrix.data.frame(geomorph::gm.prcomp(Predictor)$x)

                                                                                   base::summary(geomorph::procD.lm(Response ~ Predictor, Parallel = TRUE, iter = 1000,
                                                                                                    RRPP = TRUE))$table
                                                                                 })) %>% cbind(.[1,], .[2,1:4], .[3,1:2]) %>% .[1,])


    } # end of i comparisons for loop

    Output[[m]]$pslr.results.all <- pslr.results.all
    if(isTRUE(subset)) {
      Output[[m]]$pslr.results <- pslr.results
    }
  } # end of m loop

  setwd(Dir)

  class(Output) <- "pslr.bio"
  return(Output)
} # end of function
