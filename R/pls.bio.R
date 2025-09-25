#' @title Partial Least Squares Analysis between multivariate sets of Biological Variables
#'
#' @description
#' pls.bio is a function designed to investigate covariance relationships
#' between sets of multivariate data. This function uses a model-based hypothesis
#' testing structure where a \code{Models} nested list object contains hypothesis tests
#' which usually involve one or more covariation test sets of data. Data can be
#' both subsetted and transformed using the \code{DT} function (see details).
#' This function allows multiple data types such as geometric morphometric, genetic,
#' and biomechancial data sets to be compared in a partial least squares test.
#' Covariation tests between these sets of data within a hypothesis is ran in
#' parallel using a user-selected Partial Least Squares calculation. For each
#' covariation test within a hypothesis model, the function will generate summary
#' statistics along with an optional output for variable importance in projection
#' plots to evaluate important predictor variables.
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
#' @param paired_subsets Logical value indicating whether to analyze all possible combinations
#' of response and predictor subsets. By default, paired_subsets is FALSE. When FALSE,
#' subsets are paired (Response subset 1, Predictor subset 1; Response subset 2; Predictor subset 2;...).
#'
#' @param vips Logical value indicating whether to estimate and plot
#' variable importance in projection scores for which variables in each model's
#' Predictor data best explains the variation in the complete response data.
#'
#' @param vip_method Character string specifying PLS method to use:
#'  \itemize{
#'      \item{"spls"}   {sparce PLS via \pkg{mixOmics}}
#'      \item{"pls"}    {traditional PLS via \pkg{pls}}
#'      \item{"pcr"}    {traditional principal component regression via \pkg{pls}}
#'      \item{"plsr2"}  {two-block PLS via \pkg{plsdepot}}
#'  }
#' @param lv_method Character specifying method for selecting important variables:
#' \itemize{
#'      \item{"trim"}   {Select Predictors with VIP scores > 1 across all latent components.}
#'      \item{"median"} {Select Predictors with VIP scores greater than the median VIP scores across all latent components.}
#'      \item{"mean"}   {Select Predictors with VIP scores greater than the meean VIP scores across all latent components.}
#' }
#' @param cv_type Character specifying cross-validation type
#' @param covs Logical value indicating whether to perform covariance analysis for
#' shape datasets. This mode can handle 2D or 3D datasets and will provide covariation
#' statistics for each hypothesis in a model
#' @param key Optional character string for naming output directories
#' @param parallel Logical value indicating whether to use parallel processing
#' @param print_progress Logical value indicating whether to display progress bars
#' @param core_choice Character value specifying number of cores to use.
#' The following options are available:
#' \itemize{
#'    \item{"detect"}       {:  total available cores - 2}
#'    \item{"high"}         {:  28 cores}
#'    \item{"medium-high"}  {:   20 cores}
#'    \item{"medium"}       {:  14 cores}
#'    \item{"medium-low"}   {:  10 cores}
#'    \item{"low"}          {:   8 cores}
#' } ("detect" , "high", "medium-high", "medium", "medium-low", "low", or specific number)
#'
#' @returns A list of class "pls.bio" containing:
#' \itemize{
#'   \item{scores}          {Nested list with data tables of VIP scores for each test}
#'   \item{plsr.results}    {PLS model results including:
#'     \itemize{summary
#'       \item{R2}          {Model fit statistics}
#'       \item{Q2}          {Cross-validation prediction statistics}
#'       \item{RMSEP}       {Root mean square error of prediction}
#'       \item{p_value}     {cross-validated statistical significance value}
#'       \item{reg}         {Nested list with model coefficients and loadings}
#'     }
#'   }
#'   \item{Cov.All}         {Nested list with two-block PLS shape covariation test statistics
#'     \itemize{
#'       \item{pls2b}       {Object output from the \code{pls2B} function from \pkg{Morpho}
#'          \itemize{
#'            \item{lollipop.LV.1.Pred_and_Res}   {Optional \pkg{plotly} object when both data sets are shape arrays}
#'            \item{lollipop.LV.1.Response}       {Optional \pkg{plotly} object when Response data is shape array but not Predictor}
#'          }
#'        }
#'      }
#'    }
#'    \item{Cov.VIP}         {Nested list with two-block PLS shape covariation test statistics trimmed using VIP analysis
#'     \itemize{
#'       \item{pls2b}       {Object output from the \code{pls2B} function from \pkg{Morpho} with VIP trimming
#'          \itemize{
#'            \item{lollipop.LV.1.Pred_and_Res}   {Optional \pkg{plotly} object when both data sets are shape arrays}
#'            \item{lollipop.LV.1.Response}       {Optional \pkg{plotly} object when Response data is shape array but not Predictor}
#'          }
#'        }
#'      }
#'    }
#' }
#' @details
#' pls.bio conducts partial least squares (PLS) analyses in parallel for each test within
#' a hypothesis model. The following methods in PLS are accepted:
#'
#' \itemize{
#'   \item{Sparse PLS ("spls")} {Via mixOmics package, optimal for high-dimensional,
#'   linearly related data which requires automatic feature selection (e.g., variable
#'   trimming) using LASSO (Least Absolute Shrinkage and Selection Operator)
#'   based regularization parameters to include only the influential predictor
#'   variables in predicting the response data.}
#'   \item{Orthogonal PLS ("oscorepls")} {via pls package, optimal for
#'   high-dimensional, linearly related data which requires automatic feature
#'   selection using a singular value decomposition (i.e., PCA) based regularization
#'   parameters to only include predictor variables which are orthogonally related
#'   to the response data.}
#'   \item{Traditional Kernel PLS ("kernelpls")} {Via pls package, optimal for
#'   non-linear, multidimensional data.}
#'   \item{Wide Kernel PLS ("widekernelpls")} {via pls package, optimal for
#'   non-linear, multidimensional data with a larger variable to observation ratio.}
#'   \item{Canonically Powered PLS ("cppls")} {via pls package, optimal for
#'   maximizing overall correlations between linearly related, multivariate data.}
#'   \item{Principal Correlation Regression ("svdpc")} {via pls package, optimal
#'   for maximizing correlations between data sets that need to be subject to
#'   multidimensionality reduction (principal component analysis).}
#'   \item{Two-Block PLS ("plsr2")} {Via plsdepot package, optimal for analyzing
#'   relationships between two blocks of linearly related data.}
#'   \item{Two-Block PLS (covs = TRUE)} {via Morpho package, optimal for analyzing
#'   relationships between two blocks of multidimensional shape data.}
#' }
#'
#' The function sequentially evaluates each hypothesis model in the Models object
#' and runs, in parallel, a user-selected PLS for each test within the hypothesis
#' model. If point_set is not NULL, data will be subsetted prior to analysis. All
#' data is subjected to user-specified transformation using the \code{DT} function
#' parameters via the argument (...). If no \code{DT} function parameters are
#' provided, no data transformation will occur prior to PLS analysis.
#'
#' Partial Least Squares requires a specification for the number of latent components
#' to explore to maximize the covariance/correlation relationships between sets
#' of data. For method "spls", the \code{tune.spls} function is carried out to
#' identify the number of latent components by maximizing the squared correlation
#' coefficient (R^2) and minimizing the Root Square Mean Error of Prediction (RMSEP).
#' Additionally, only with method "spls", an autotuning procedure is taken based
#' on the number of predictor and response variables in the test to identify the
#' best number of meaningful variables in both the predictor and response that
#' produces the highest squared correlation coefficient (R^2) as well as the
#' lowest modelerror (RMSEP). Autotuning is a feature specific to the SPLS methodology
#' in feature reduction, but can be turned off by setting autotune = FALSE.
#' Similarly, for all other methods, an internal \code{NComp} function is run to
#' identify the number of meaningful latent components following the same criteria
#' but without feature autotuning.
#'
#' This function also supports an optional two-block PLS method through the
#' \pkg{Morpho} package particularly for analyzes of integration tests. If geometric
#' morphometric data is supplied in a "p x k x n" (landmarks, dimensions, observations)
#' array, the function will also store interactive .html plots using the \pkg{plotly}
#' package of the covariance relationships of the response data. Additionally,
#' summary statistics will be calculated for each hypothesis test.
#' The following statistics are provided:
#'
#' \itemize{
#'    \item{df}                 {List of all variables in the predictor data}
#'    \item{PLS2B}              {Data table of the latent variable weights and correlations of the test}
#'    \item{All RV Coefficient} {The RV integration statistic for the test based on all predictor data}
#'    \item{All RV P-Value}     {The p-value associated with the test on all prdictor data}
#'    \item{VIP Vars}           {List of variables after VIP selection}
#'    \item{VIP PLS2B}          {Data table of the latent variable weights and correlations of the test after VIP trimming}
#'    \item{VIP RV Coefficient} {The RV integration statistic for the test based on VIP trimmed predictor data}
#'    \item{VIP RV P-Value}     {The p-value associated with the test on VIP trimmed prdictor data}
#' }
#'
#' As each test is performed, an internal function \code{plsr_stats} is conducted
#' to produce cross-validated summary statistics including an overall squared (R^2)
#' and cross-validated (Q2) correlation coefficient as well as the Root Mean Squared
#' Error of Prediction (RMSEP) and a statistical significance value (p-value).
#'
#' Variable importance in pls.bio is assessed through Variable Importance in
#' Projection scores (VIPs), calculated by the weight that each variable of the
#' predictor data which explains the most variance across the selected latent variables.
#' Variables can be selected using the lv_method argument based on different
#' thresholds (VIP > 1, median VIP, or mean VIP) across all selected latent variables.
#' VIP plots are subsequently generated and stored in a directory folder called
#' "PLS Regression - VIP Scores" with an optionally specified key argument string
#' to avoid data overwriting when running more than once. To stop VIP plot generation,
#' set vips = FALSE.
#'
#' This function has a progress bar built in and can be removed by setting
#' print.progress = FALSE.
#'
#' Parallel processing is used to maxmize computational efficiency. As such, the
#' number of cores available for analysis can be automatically detected by using
#' "detect" string in argument core_choice or as one of the options defined in the
#' core_choice argument.
#'
#'
#' @seealso \code{\link{AABA}} \code{\link{NComp}} \code{\link{plsr_stats}}
#'
#' @author Brian Anthony Keeling
#'
#' @importFrom mixOmics spls tune.spls
#' @importFrom pls mvr plsr
#' @importFrom plsdepot plsreg2
#' @importFrom caret createDataPartition
#' @importFrom Morpho pls2B plsCoVar plsCoVarCommonShape
#' @importFrom dplyr filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom grDevices png dev.off graphics.off
#' @importFrom ggplot2 ggplot aes geom_boxplot stat_summary geom_hline ggtitle theme element_text
#'
#' @references
#' Rohart F, Gautier B, Singh A, Lê Cao K-A (2017). "mixOmics: An R package for 'omics feature
#' selection and multiple data integration." PLoS Computational Biology, 13(11), e1005752.
#'
#' Lê Cao K-A, Rossouw D, Robert-Granié C, Besse P (2008). "A sparse PLS for variable selection
#' when integrating Omics data." Statistical Applications in Genetics and Molecular Biology, 7(1), 35.
#'
#' Mevik, B.H., Wehrens, R. and Liland, K.H. (2019). pls: Partial Least Squares and Principal
#' Component Regression. R package version 2.7-3.
#'
#' Sanchez, G. (2012). plsdepot: Partial Least Squares (PLS) Data Analysis Methods.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Scenario Example:
#' # Question: Which variables and to what degree is mandibular symphyseal shape
#' # influenced by the cross-sectional properties of masticatory regions of the
#' # posterior corpus?
#'
#' # Create hypothesis model structure
#' Models <- list(
#'   'Symphyseal_Shape ~ Posterior_Corpus_Properties' = list( # This is your hypothesis model
#'     'Symphysis_LM1M2' = list( # This is your first hypothesis test
#'       'Response' = Corpus_Land[241:360,,],  # Symphysis landmarks
#'       'Predictor' = Biomechanics %>%
#'         filter(Region == "LM1M2") %>%
#'         select(5:93)  # Biomechanical variables,
#'        'dt_parameters' =  list(Res_transform = "NULL", # No data transformation. Default is "NULL"
#'                                Pred_transform = "NULL", # No data transformation. Default is "NULL"
#'                                point_set_res = NULL, # Not necessary since we will add a point_set but useful for single subset cases
#'                                point_set_pred = NULL
#'                                )
#'     ),
#'     'Symphysis_LP3P4' = list(
#'       'Response' = Corpus_Land[241:360,,],
#'       'Predictor' = Biomechanics %>%
#'         filter(Region == "LP3P4") %>%
#'         select(5:93),
#'       'dt_parameters' = list() # adding this is not necessary if not transforming data prior to analysis.
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
#' point_set <- list('Symphyseal Shape ~ Posterior_Corpus_Properties' =
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
#' # Run covariation analysis
#' results <- pls.bio(
#'   Models = Models, # add in hypothesis models
#'   point_set = point_set, # add in subsets
#'   paired_subsets = FALSE, # since the Response data stays the same, its unnecessary to do cross comparisons
#'   covs = TRUE # since we have 3D shape data, let's view the covarying shape changes by the predictor variables
#' ) # keep defaults
#'
#' # Identify the VIP scores for the first subsetted test ('Symphysis_Shape ~ LM1M2') in the first hypothesis model ('Symphysis_Shape ~ Posterior_Corpus_Properties')
#'
#' results[[1]]$scores[[1]][[1]]
#'
#' # Explanation
#' # results is the output
#' # [[1]] after results refers to the first hypothesis model ('Symphysis_Shape ~ Posterior_Corpus_Properties')
#' # $scores refers to the VIP scores
#' # [[1]] after $scores refers to the first test in the first hypothesis model ('Symphysis_Shape ~ LM1M2')
#' # [[1]] after [[1]] refers to the first subset test in the first test in the first hypothesis model ('Symphysis_Shape ~ LM1M2' -- "Corpus Dimensions")
#'
#' # You can view the VIP plot generated for this test as well!
#' library(magick)
#' VIP_plot <- read.image("/PLS VIP Scores/Symphysis_Shape ~ LM1M2_Subset_1_Dimensions.png")
#'
#' # View covariance results between symphyseal shape and corpus dimensions of the LM1M2 corpus region.
#'
#' results[[1]]$COV.VIP[[1]][[1]]$pls2b
#'
#' # View interactive plot for the first latent variable
#' # NOTE: all plots are stored in file directory
#'
#' results[[1]]$COV.VIP[[1]][[1]]$pls2b$lollipop.LV.1.Response
#'
#' # View covariance results from the excel output
#'
#' view(
#'      read_xlsx("/PLS VIP Scores/2BPLS_Subset_1_Dimensions_Model_1.xlsx")
#'      )
#' }


pls.bio <- function(Models = NULL,
                    point_set = NULL,
                    paired_subsets = FALSE,
                    vips = TRUE,
                    vip_method = "spls",
                    lv_method = "mean",
                    cv_type = "CV",
                    covs = FALSE,
                    key = NULL,
                    parallel = TRUE,
                    print_progress = TRUE,
                    core_choice = "detect",
                    ...) {

  set.seed(1)
  output = vector("list", length(Models))

  Dir = getwd()

  if(!is.null(key)) {
    if(!dir.exists(paste("PLS VIP Scores -", key))) {
      dir.create(paste("PLS VIP Scores -", key))
    }
    pls.dir = paste("PLS VIP Scores -", key)
  } else {
    if(!dir.exists("PLS VIP Scores")) {
      dir.create("PLS VIP Scores")
    }
    pls.dir = "PLS VIP Scores"
  }

  path_pls = paste(Dir, "/", pls.dir, sep = "")
  path_pls = file.path(path_pls)
  setwd(path_pls)

  #______________________________________________________________________________________________________________________
  if (isTRUE(print_progress)) {
    pb <- progress::progress_bar$new(
      format = "Processing :what [:bar] :percent ETA: :eta",
      clear = FALSE, total = 100, width = 80)
  } else {}

  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = "Starting pls.bio Modeling"))
  }
  #______________________________________________________________________________________________________________________

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

  #_______________________________________________________________________________________________________________________

  subset = if(!is.null(point_set)) TRUE else FALSE

  #_______________________________________________________________________________________________________________________

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


  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = "pls.bio paramters set"))
  }

  #_______________________________________________________________________________________________________________________

  for(m in seq_along(Models)) {

    if(!dir.exists(paste(names(Models)[m]))) {
      dir.create(paste(names(Models)[m]))
    }
    model.dir = paste(names(Models)[m])
    if(!dir.exists(paste(model.dir))) {
      dir.create(paste(model.dir))
    }
    path = file.path(Dir, pls.dir, model.dir)

    if(isTRUE(print_progress)) {
      pb$tick(10/length(Models), tokens = list(what = paste("Model", m, "PLS and VIP Calculations:")))
    }

    plsr.results = vector("list", length(Models[[m]]))
    plsr.results.all = vector("list", length(Models[[m]]))
    plots.vip = vector("list", length(Models[[m]]))
    plots.vip.all = vector("list", length(Models[[m]]))

    if(isTRUE(subset)) {
      for(i in seq_along(point_set[[m]])) {
        if(isTRUE(paired_subsets)) {
          p_subset = expand.grid(seq_along(point_set[[m]][[i]][[1]]), seq_along(point_set[[m]][[i]][[2]]))
          plots.vip[[i]] <- vector("list", nrow(p_subset))
          plsr.results[[i]] <- vector("list", nrow(p_subset))
        } else {
          plots.vip[[i]] <- vector("list", length(point_set[[m]][[i]][[2]]))
          plsr.results[[i]] <- vector("list", length(point_set[[m]][[i]][[2]]))
        }
      }
    }

    if (isTRUE(parallel)) {
      clust <- parallel::makeCluster(core_choice)
      doParallel::registerDoParallel(clust)
      #PPARAM <- BiocParallel::DoparParam()
    } else {
      foreach::registerDoSEQ()
    }

    parallel.vip <- foreach(i = 1:length(Models[[m]]), .packages =
                              c("factoextra", "FactoMineR", "geomorph","grDevices", "graphics", "caret", "doParallel", "foreach", "BiocParallel", "parallel", "stats", "BiocManager",
                                "doParallel", "Hmisc", "Morpho", "pls", "mixOmics", "plsdepot", "tidyverse", "corrplot", "CSGM")) %dopar% {

      if(isFALSE(isTRUE(file.access(path, 2) == 0))) {
        setwd(path)
      }

      if(isTRUE(subset)) {

        if(isTRUE(paired_subsets)) {
          p_subset = expand.grid(seq_along(point_set[[m]][[i]][[1]]), seq_along(point_set[[m]][[i]][[2]]))
        } else {
          p_subset = data.frame(Response = rep(1:length(point_set[[m]][[i]][[1]])),
                                Predictor = rep(1:length(point_set[[m]][[i]][[2]]))
          )
        }

        for(j in 1:nrow(p_subset)) {

          dt_parameters = tryCatch({modifyList(dt_parameters,Models[[m]][[i]][[3]])},
                   error = function(e){dt_parameters})

          dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                         Predictor = Models[[m]][[i]][[2]],
                                                         Res_point_set = point_set[[m]][[i]][[1]][[p_subset[j,1]]],
                                                         Pred_point_set = point_set[[m]][[i]][[2]][[p_subset[j,2]]],
                                                         Res_transform = if(length(dim(Models[[m]][[i]][[1]])) == 3 & dt_parameters$Res_transform == "NULL") "2D" else dt_parameters$Res_transform,
                                                         Pred_transform = if(length(dim(Models[[m]][[i]][[2]])) == 3 & dt_parameters$Pred_transform == "NULL") "2D" else dt_parameters$Pred_transform
                                                        )) %>% .[intersect(names(.), dt_params)]

          Response <- data.frame(do.call(DT, dt_parameters)$Response)
          Predictor <- data.frame(do.call(DT, dt_parameters)$Predictor)

          if(!vip_method == "plsr2") {

            if(vip_method == "spls") {

              ncomp <- NComp(model = list(Response = Response, Predictor = Predictor), vip_method = vip_method)
              reg <- mixOmics::spls(Predictor, Response, ncomp = ncomp, mode = "regression", near.zero.var = TRUE)
              reg$model <- list(Y = Response, X = Predictor)
              keepX = rep(ncol(VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)),ncomp)
              reg <- mixOmics::spls(Response, Predictor, ncomp = ncomp, mode = "regression", near.zero.var = TRUE)
              reg$model <- list(Y = Predictor, X = Response)
              keepY = rep(ncol(VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)),ncomp)

              reg <- mixOmics::spls(Predictor, Response, ncomp = ncomp, keepX = keepX[1:ncomp], keepY = keepY[1:ncomp], mode = "regression", near.zero.var = TRUE)
              reg$model <- list(Y = Response, X = Predictor)
              scores <- data.frame(VIP(reg, Scores = TRUE, ncomp = ncomp, vip_threshold = "NULL"))
              plots.vip[[i]][[j]]$scores <- scores
              vip <- VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)
              plsr.results[[i]][[j]]$summary <- data.frame(plsr_stats(reg, ncomp, vip, cv_type = cv_type, parallel.me = FALSE))
              plsr.results[[i]][[j]]$vip <- vip
              plsr.results[[i]][[j]]$ncomp <- ncomp

              rm(vip, ncomp, reg)
              gc()

            } else {
              ncomp <- NComp(model = list(Response = Response, Predictor = Predictor), vip_method = vip_method)
              reg <- pls::mvr(as.matrix.data.frame(Response) ~ as.matrix.data.frame(Predictor), ncomp = ncomp, validation = "CV", scale = FALSE, method = vip_method)
              reg$model <- list(Y = as.matrix.data.frame(Response), X = as.matrix.data.frame(Predictor))
              scores <- data.frame(VIP(reg, Scores = TRUE, ncomp = ncomp, vip_threshold = "NULL"))
              plots.vip[[i]][[j]]$scores <- scores
              vip <- VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)
              plsr.results[[i]][[j]]$summary <- data.frame(plsr_stats(reg, ncomp, vip, cv_type = cv_type, parallel.me = FALSE))
              plsr.results[[i]][[j]]$vip <- vip
              plsr.results[[i]][[j]]$ncomp <- ncomp
              rm(vip, ncomp, reg)
              gc()
            }
          }

          if(vip_method == "plsr2") {

            ncomp <- NComp(model = list(Response = Response, Predictor = Predictor), vip_method = vip_method)
            reg <- plsreg2(responses = as.matrix.data.frame(Response), predictors = as.matrix.data.frame(Predictor), comps = ncomp, crosval = "CV")
            reg$model <- list(Y = as.matrix.data.frame(Response), X = as.matrix.data.frame(Predictor))
            scores <- data.frame(VIP(reg, Scores = TRUE, ncomp = ncomp, vip_threshold = "NULL"))
            plots.vip[[i]][[j]]$scores <- scores
            vip <- VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)
            plsr.results[[i]][[j]]$summary <- data.frame(plsr_stats(reg, ncomp, vip, cv_type = cv_type, parallel.me = FALSE))
            plsr.results[[i]][[j]]$vip <- vip
            plsr.results[[i]][[j]]$ncomp <- ncomp
            rm(ncomp, reg,vip)
            gc()
          } # end of plsr2 calculation

          #Plot generation________________________________________________________________________________________

          if(isTRUE(vips)){

            if(length(unique(scores$Predictor)) > 100){
              nsize = 1
            }else if(length(unique(scores$Predictor)) > 80) {
              nsize = 2
            }else if(length(unique(scores$Predictor)) > 50) {
              nsize = 3
            } else if(length(unique(scores$Predictor)) < 50) {
              nsize = 4
            }

            placeholder <- scores %>%
              mutate(Predictor = factor(Predictor, levels = unique(Predictor))) %>%
              ggplot(aes(x = Predictor, y = Score)) +
              geom_boxplot(outlier.shape = NA) +
              stat_summary(fun = median, geom = "line", aes(group = 1),
                           colour = "blue", linetype = "solid", size = 0.8) +
              geom_hline(aes(yintercept = 1), color = "red", linewidth = 1.5) +
              ggtitle("Boxplot of Predictor Variables by VIP Scores") +
              ggpubr::theme_pubclean() +
              theme(axis.text.x = element_text(angle = 90, size = nsize),
                    axis.text.y = element_text(size = nsize),
                    plot.title = element_text(face = "bold", hjust = 0.5, size = 7),
                    axis.title.x = element_text(size = 7),
                    axis.title.y = element_text(size = 7))

            filename <- file.path(path, paste("Model", m, "Test", i, "Subset", j, ".png", sep = "_"))
            grDevices::png(filename = filename, res = 600, width = 1920, height = 1080)
            print(placeholder)
            dev.off()

            place.group = vector("list", length(plots.vip.all[[i]]$Response))

            for(r in seq_along(place.group)) {
              place.group[[r]] <- scores %>%
                filter(Response == paste("t", r, sep = "")) %>%
                mutate(Predictor = factor(Predictor, levels = unique(Predictor))) %>%
                ggplot(aes(x = Predictor, y = Score)) +
                geom_boxplot(outlier.shape = NA) +
                stat_summary(fun = median, geom = "line", aes(group = 1),
                             colour = "blue", linetype = "solid", size = 0.8) +
                geom_hline(aes(yintercept = 1), color = "red", linewidth = 1) +
                ggtitle(paste("Latent Variable", r, sep = "_")) +
                ggpubr::theme_pubclean() +
                theme(axis.text.x = element_text(angle = 90, size = nsize),
                      axis.text.y = element_text(size = nsize),
                      plot.title = element_text(face = "bold", hjust = 0.5, size = 5),
                      axis.title.x = element_text(size = 5),
                      axis.title.y = element_text(size = 5))
            } # end r for loop

            arrange = NULL

            if(length(place.group) == 2) {

              arrange <- ggpubr::ggarrange(place.group[[1]], place.group[[2]], nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")
              LV = 2

            } else if(length(place.group) == 3) {

              arrange <- ggpubr::ggarrange(place.group[[1]], place.group[[2]], place.group[[3]], nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
              LV = 3

            } else if(length(place.group) >= 4) {

              arrange <- ggpubr::ggarrange(place.group[[1]], place.group[[2]], place.group[[3]],
                                   place.group[[4]], nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
              LV = 4
            }

            if(!is.null(arrange)) {
              arrange <- ggpubr::annotate_figure(arrange, top = ggpubr::text_grob(paste("VIP Scores of: ", names(Models[[m]][[i]])[[2]],
                                                                        "for LVs 1-",LV, sep = ""), color = "black", face = "bold", size = 7))

              filename <- file.path(path, paste("Model", m, "Test", i, "Subset", j, "of LVs", ".png", sep = "_"))
              grDevices::png(filename = filename, res = 600, width = 1920, height = 1080)

              print(arrange)
              graphics.off()
            }
          }
        #end of plots_________________________________________________________________________________________________
        } # end of j for loop
      }

   #all____________________________________________________________________________________________________________
      dt_parameters = tryCatch({modifyList(dt_parameters,Models[[m]][[i]][[3]])},
                      error = function(e){dt_parameters})

        dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                       Predictor = Models[[m]][[i]][[2]],
                                                       Res_transform = if(length(dim(Models[[m]][[i]][[1]])) == 3 & dt_parameters$Res_transform == "NULL") "2D" else dt_parameters$Res_transform,
                                                       Pred_transform = if(length(dim(Models[[m]][[i]][[2]])) == 3 & dt_parameters$Pred_transform == "NULL") "2D" else dt_parameters$Pred_transform,
                                                       Pred_point_set = NULL,
                                                       Res_point_set = NULL)) %>% .[intersect(names(.), dt_params)]

        Response <- data.frame(do.call(DT, dt_parameters)$Response)
        Predictor <- data.frame(do.call(DT, dt_parameters)$Predictor)

        if(!vip_method == "plsr2") {

          if(vip_method == "spls") {

            ncomp <- NComp(model = list(Response = Response, Predictor = Predictor), vip_method = vip_method)
            reg <- mixOmics::spls(Predictor, Response, ncomp = ncomp, mode = "regression", near.zero.var = TRUE)
            reg$model <- list(Y = Response, X = Predictor)
            keepX = rep(ncol(VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)),ncomp)
            reg <- mixOmics::spls(Response, Predictor, ncomp = ncomp, mode = "regression", near.zero.var = TRUE)
            reg$model <- list(Y = Predictor, X = Response)
            keepY = rep(ncol(VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)),ncomp)

            reg <- mixOmics::spls(Predictor, Response, ncomp = ncomp, keepX = keepX[1:ncomp], keepY = keepY[1:ncomp], mode = "regression", near.zero.var = TRUE)
            reg$model <- list(Y = Response, X = Predictor)
            scores <- data.frame(VIP(reg, Scores = TRUE, ncomp = ncomp, vip_threshold = "NULL"))
            plots.vip.all[[i]]$scores <- scores
            vip <- VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)
            plsr.results.all[[i]]$summary <- data.frame(plsr_stats(reg, ncomp, vip, cv_type = cv_type, parallel.me = FALSE))
            plsr.results.all[[i]]$vip <- vip
            plsr.results.all[[i]]$ncomp <- ncomp

            rm(vip, ncomp, reg)
            gc()

          } else {
            ncomp <- NComp(model = list(model = list(Response = Response, Predictor = Predictor)), vip_method = vip_method)
            reg <- pls::mvr(as.matrix.data.frame(Response) ~ as.matrix.data.frame(Predictor), ncomp = ncomp, validation = "CV", scale = FALSE, method = vip_method)
            reg$model <- list(Y = as.matrix.data.frame(Response), X = as.matrix.data.frame(Predictor))
            scores <- data.frame(VIP(reg, Scores = TRUE, ncomp = ncomp, vip_threshold = "NULL"))
            plots.vip.all[[i]]$scores <- scores
            vip <- VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)
            plsr.results.all[[i]]$summary <- data.frame(plsr_stats(reg, ncomp, vip, cv_type = cv_type, parallel.me = FALSE))
            plsr.results.all[[i]]$vip <- vip
            plsr.results.all[[i]]$ncomp <- ncomp
            rm(vip, ncomp, reg)
            gc()
          }
        }

        if(vip_method == "plsr2") {

          ncomp <- NComp(model = list(model = list(Response = Response, Predictor = Predictor)), vip_method = vip_method, repeats = 20)
          reg <- plsreg2(responses = as.matrix.data.frame(Response), predictors = as.matrix.data.frame(Predictor), comps = ncomp, crosval = "CV")
          reg$model <- list(Y = as.matrix.data.frame(Response), X = as.matrix.data.frame(Predictor))
          scores <- data.frame(VIP(reg, Scores = TRUE, ncomp = ncomp, vip_threshold = "NULL"))
          plots.vip.all[[i]]$scores <- scores
          vip <- VIP(reg, Scores = FALSE, ncomp = ncomp, vip_threshold = lv_method)
          plsr.results.all[[i]]$summary <- data.frame(plsr_stats(reg, ncomp, vip, cv_type = cv_type, parallel.me = FALSE))
          plsr.results.all[[i]]$vip <- vip
          plsr.results.all[[i]]$ncomp <- ncomp
          rm(ncomp, reg, vip)
          gc()
        } # end of plsr2 calculation
      #_Plot generation___________________________________________________________________________________________
      if(isTRUE(vips)){

        if(length(unique(scores$Predictor)) > 100){
          nsize = 1
        }else if(length(unique(scores$Predictor)) > 80) {
          nsize = 2
        }else if(length(unique(scores$Predictor)) > 50) {
          nsize = 3
        } else if(length(unique(scores$Predictor)) < 50) {
          nsize = 4
        }

        placeholder.all <- scores %>%
          mutate(Predictor = factor(Predictor, levels = unique(Predictor))) %>%
          ggplot(aes(x = Predictor, y = Score)) +
          geom_boxplot(outlier.shape = NA) +
          stat_summary(fun = median, geom = "line", aes(group = 1),
                       colour = "blue", linetype = "solid", size = 0.8) +
          geom_hline(aes(yintercept = 1), color = "red", linewidth = 1.5) +
          ggtitle("Boxplot of Predictor Variables by VIP Scores") +
          ggpubr::theme_pubclean() +
          theme(axis.text.x = element_text(angle = 90, size = nsize),
                axis.text.y = element_text(size = nsize),
                plot.title = element_text(face = "bold", hjust = 0.5, size = 7),
                axis.title.x = element_text(size = 7),
                axis.title.y = element_text(size = 7))

        filename <- file.path(path, paste("Model", m, "Test", i, ".png", sep = "_"))
        grDevices::png(filename = filename, res = 600, width = 1920, height = 1080)
        print(placeholder.all)
        dev.off()

        place.group.all = vector("list", length(plots.vip.all[[i]]$Response))

        for(r in seq_along(place.group.all)) {
          place.group.all[[r]] <- scores %>%
            filter(Response == paste("t", r, sep = "")) %>%
            mutate(Predictor = factor(Predictor, levels = unique(Predictor))) %>%
            ggplot(aes(x = Predictor, y = Score)) +
            geom_boxplot(outlier.shape = NA) +
            stat_summary(fun = median, geom = "line", aes(group = 1),
                         colour = "blue", linetype = "solid", size = 0.8) +
            geom_hline(aes(yintercept = 1), color = "red", linewidth = 1) +
            ggtitle(paste("Latent Variable", r, sep = "_")) +
            ggpubr::theme_pubclean() +
            theme(axis.text.x = element_text(angle = 90, size = nsize),
                  axis.text.y = element_text(size = nsize),
                  plot.title = element_text(face = "bold", hjust = 0.5, size = 5),
                  axis.title.x = element_text(size = 5),
                  axis.title.y = element_text(size = 5))
        } # end r for loop

        arrange = NULL

        if(length(place.group.all) == 2) {

          arrange <- ggpubr::ggarrange(place.group.all[[1]], place.group.all[[2]], nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")
          LV = 2

        } else if(length(place.group.all) == 3) {

          arrange <- ggpubr::ggarrange(place.group.all[[1]], place.group.all[[2]], place.group.all[[3]], nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
          LV = 3

        } else if(length(place.group.all) >= 4) {

          arrange <- ggpubr::ggarrange(place.group.all[[1]], place.group.all[[2]], place.group.all[[3]],
                                place.group.all[[4]], nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
          LV = 4
        }

        if(!is.null(arrange)) {
          arrange <- ggpubr::annotate_figure(arrange, top = ggpubr::text_grob(paste("VIP Scores of: ", names(Models[[m]][[i]])[[2]],
                                                                    "for LVs 1-",LV, sep = ""), color = "black", face = "bold", size = 7))

          filename <- file.path(path, paste("Model", m, "Test", i, "of LVs", ".png", sep = "_"))
          grDevices::png(filename = filename, res = 600, width = 1920, height = 1080)

          print(arrange)
          graphics.off()
        }
      }
      #__end of plots generation_______________________________________________________________________________________________________


      if(isTRUE(subset)){
        return(list(plots.vip.all = plots.vip.all,
                    plots.vip = plots.vip,
                    plsr.results.all = plsr.results.all,
                    plsr.results = plsr.results))
      } else {
        return(list(plots.vip.all = plots.vip.all,
                    plsr.results.all = plsr.results.all))
      }
    } # end of i parallel for loop

    if(isTRUE(parallel)) {
      stopCluster(cl = clust)
    }

    clean_and_flatten <- function(x) {
      # Remove NULL or empty lists
      x <- Filter(function(x) !is.null(x) && length(x) != 0, x)

      x <- lapply(x, function(z) Filter(function(x) !is.null(x), z))

      x <- Filter(function(d) !is_empty(d), x)

      # If after filtering, the list contains only one item and it's a list, return it directly
      if (length(x) == 1 && is.list(x[[1]])) {
        return(clean_and_flatten(x[[1]]))
      } else {
        # Otherwise, recurse through any lists within this list and apply the same cleaning
        return(lapply(x, function(item) {
          if (is.list(item)) {
            clean_and_flatten(item)
          } else {
            item
          }
        }))
      }
    }

    if(isTRUE(print_progress)) {
      pb$tick(40/length(Models), tokens = list(what = paste("Model", m, "Finished PLS and VIP Calculations:")))
    }

    # Applying the cleaning function to the overall list
    cleaned_list <- lapply(parallel.vip, clean_and_flatten)

    # Example to assign cleaned data to the global environment for each statistic name
    var_boot = if(isTRUE(subset)) c("plots.vip.all", "plots.vip", "plsr.results.all" , "plsr.results") else c("plots.vip.all", "plsr.results.all")
    for (stat_name in var_boot) {
      assign(stat_name, lapply(cleaned_list, `[[`, stat_name))
    }

    rm(parallel.vip, var_boot, cleaned_list)
    gc()

    # Analysis of Covariance for Landmark Data______________________________________________________________________________________________________________________
    if (isTRUE(print_progress)) {
      pb$tick(10/length(Models), tokens = list(what = paste("Model", m, "Optional Shape Covariances:")))
    }
    if(isTRUE(covs)) {

      LV.VIP.all <- vector("list", length(plots.vip.all))
      LV.all <- vector("list", length(plots.vip.all))

      if(isTRUE(subset)) {
        LV.VIP <- vector("list", length(plots.vip))
        LV <- vector("list", length(plots.vip))
        for(i in 1:length(plots.vip)) {
          LV.VIP[[i]] <- vector("list", length(plots.vip[[i]]))
          LV[[i]] <- vector("list", length(plots.vip[[i]]))
        }
      }

      for(i in seq_along(plots.vip)) {

        LV.all[[i]] <- vector("list", 4)
        names(LV.all[[i]]) <- c("summary", "scores", "shapes", "vars")

        LV.VIP.all[[i]] <- vector("list", 4)
        names(LV.VIP.all[[i]]) <- c("summary", "scores", "shapes", "vars")

        if(isTRUE(subset)) {
          for(j in seq_along(plots.vip[[i]])) {
            LV[[i]][[j]] <- vector("list", 4)
            names(LV[[i]][[j]]) <- c("summary", "scores", "shapes", "vars")

            LV.VIP[[i]][[j]] <- vector("list", 4)
            names(LV.VIP[[i]][[j]]) <- c("summary", "scores", "shapes", "vars")
          } # end of j for loop
        }
      }  # end of i for loop

      if (isTRUE(parallel)) {
        doParallel::registerDoParallel(cores = core_choice)
      }

      i = 1

      parallel.cov <- foreach(i = 1:length(plots.vip), .combine = "rbind", .packages =
                                c("factoextra", "FactoMineR", "geomorph","grDevices", "graphics", "vegan",
                                  "Hmisc", "Morpho","tidyverse", "RColorBrewer", "plotly", "htmlwidgets", "CSGM")) %dopar% {

        if(isTRUE(subset)) {

          if(isTRUE(paired_subsets)) {
            p_subset = expand.grid(seq_along(point_set[[m]][[i]][[1]]), seq_along(point_set[[m]][[i]][[2]]))
          } else {
            p_subset = data.frame(Response = rep(1:length(point_set[[m]][[i]][[1]])),
                                  Predictor = rep(1:length(point_set[[m]][[i]][[2]]))
            )
          }
            for(j in 1:nrow(p_subset)) {

              tryCatch({dt_parameters = modifyList(dt_parameters,Models[[m]][[i]][[3]])},
                       error = function(e){})

              if(length(dim(Models[[m]][[i]][[1]])) == 3){
                dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                               Predictor = Models[[m]][[i]][[2]],
                                                               Res_point_set = point_set[[m]][[i]][[1]][[p_subset[j,1]]],
                                                               Pred_point_set = point_set[[m]][[i]][[2]][[p_subset[j,2]]],
                                                               Res_transform = "NULL")
                                            ) %>% .[intersect(names(.), dt_params)]


              }
              if(length(dim(Models[[m]][[i]][[2]])) == 3) {
                dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                               Predictor = Models[[m]][[i]][[2]],
                                                               Res_point_set = point_set[[m]][[i]][[1]][[p_subset[j,1]]],
                                                               Pred_point_set = point_set[[m]][[i]][[2]][[p_subset[j,2]]],
                                                               Pred_transform = "NULL")
                                            ) %>% .[intersect(names(.), dt_params)]
              }
              if(length(dim(Models[[m]][[i]][[1]])) != 3 & length(dim(Models[[m]][[i]][[2]])) != 3) {
                dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                               Predictor = Models[[m]][[i]][[2]],
                                                               Res_point_set = point_set[[m]][[i]][[1]][[p_subset[j,1]]],
                                                               Pred_point_set = point_set[[m]][[i]][[2]][[p_subset[j,2]]])
                                          ) %>% .[intersect(names(.), dt_params)]
              }

              Response = do.call(DT, dt_parameters)$Response
              Predictor = do.call(DT, dt_parameters)$Predictor

              if(length(dim(Response)) == 3 & length(dim(Predictor)) == 3) {
                config = TRUE
              } else {
                config = FALSE
              }
              placeholder = pls2B(y = Response, x = Predictor, same.config = config, cv = TRUE, rounds = 100)

              LV[[i]][[j]]$summary <- data.frame(df = ncol(Predictor), RV = placeholder$rv, p_value = placeholder$p.value.RV)

              if(isTRUE(config)) {

                shape = plsCoVarCommonShape(placeholder, 1)

                LV[[i]][[j]]$shapes$lollipop.LV.1.Pred_and_Res <- lolshape(shape.list = list(X_neg = array(shape[,,1], dim = c(dim(shape)[1], 3, 1)),
                                                                                                     X_pos = array(shape[,,2], dim = c(dim(shape)[1], 3, 1))),
                                                                                   ID = c("Pred_&_Res_Neg", "Pred_&_Res_Pos"), title = "Covariance between Predictor and Response Variables for LV 1")
                saveWidget(LV[[i]][[j]]$shapes$lollipop.LV.1.Pred_and_Res, paste("Cov.LV.1", "_","Model", m, "_", names(Models[[m]])[i], "_", "sub", "_", j, ".html", sep = ""))

                shape = plsCoVarCommonShape(placeholder, 2)
                LV[[i]][[j]]$shapes$lollipop.LV.2.Pred_and_Res <- lolshape(shape.list = list(X_neg = array(shape[,,1], dim = c(dim(shape)[1], 3, 1)),
                                                                                                     X_pos = array(shape[,,2], dim = c(dim(shape)[1], 3, 1))),
                                                                                   ID = c("Predictor_Neg", "Response_Neg"), title = "Covariance between Predictor and Response Variables for LV 2")
                saveWidget(LV[[i]][[j]]$shapes$lollipop.LV.2.Pred_and_Res, paste("Cov.LV.2", "_","Model", m, "_", names(Models[[m]])[i], "_", "sub", "_", j, ".html", sep = ""))

              } else {

                if(length(dim(Response)) == 3) {

                  shape = plsCoVar(placeholder, i = 1, sdx = 2, sdy = 2)
                  LV[[i]][[j]]$shapes$lollipop.LV.1.Response <- lolshape(shape.list = list(Y_neg = array(shape$y[,,1], dim = c(dim(shape$y)[1],3,1)),
                                                                                                   Y_pos = array(shape$y[,,2], dim = c(dim(shape$y)[1],3,1))),
                                                                                 ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Response Variable for LV 1")
                  shape = plsCoVar(placeholder, i = 2, sdx = 2, sdy = 2)
                  LV[[i]][[j]]$shapes$lollipop.LV.2.Response <- lolshape(shape.list = list(Y_neg = array(shape$y[,,1], dim = c(dim(shape$y)[1],3,1)),
                                                                                                   Y_pos = array(shape$y[,,2], dim = c(dim(shape$y)[1],3,1))),
                                                                                 ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Response Variable for LV 2")
                }
                if(length(dim(Predictor)) == 3) {
                  shape = plsCoVar(placeholder, i = 1, sdx = 2, sdy = 2)
                  LV[[i]][[j]]$shapes$lollipop.LV.1.Predictor <- lolshape(shape.list = list(X_neg = array(shape$x[,,1], dim = c(dim(shape$x)[1],3,1)),
                                                                                                    X_pos = array(shape$x[,,2], dim = c(dim(shape$x)[1],3,1))),
                                                                                  ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Predictor Variable for LV 1")
                  shape = plsCoVar(placeholder, i = 2, sdx = 2, sdy = 2)
                  LV[[i]][[j]]$shapes$lollipop.LV.2.Predictor <- lolshape(shape.list = list(X_neg = array(shape$x[,,1], dim = c(dim(shape$x)[1],3,1)),
                                                                                                    X_pos = array(shape$x[,,2], dim = c(dim(shape$x)[1],3,1))),
                                                                                  ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Predictor Variable for LV 2")

                }
              } # end of lollipop generation for subsets
              #______________________________________________________________________________________________________

              LV.VIP[[i]][[j]]$vars <- colnames(data.frame(plsr.results[[i]][[j]]$vip))
              Predictor = data.frame(plsr.results[[i]][[j]]$vip)
              Predictor = if(ncol(Predictor) == 0) NULL else Predictor
              #if there are selected variables ______________________________________________________________________________________________________
              if(!is.null(Predictor)) {

                if(length(dim(Response)) == 3 & length(dim(Predictor)) == 3) {
                  config = TRUE
                } else {
                  config = FALSE
                }
                placeholder = pls2B(y = Response, x = Predictor, same.config = config, cv = TRUE, rounds = 100)

                LV.VIP[[i]][[j]]$summary <- data.frame(VIP_df = ncol(Predictor), VIP_RV = placeholder$rv, VIP_p_value = placeholder$p.value.RV)

                if(isTRUE(config)) {

                  shape = plsCoVarCommonShape(placeholder, 1)

                  LV.VIP[[i]][[j]]$shapes$lollipop.LV.1.Pred_and_Res <- lolshape(shape.list = list(X_neg = array(shape[,,1], dim = c(dim(shape)[1], 3, 1)),
                                                                                                   X_pos = array(shape[,,2], dim = c(dim(shape)[1], 3, 1))),
                                                                                 ID = c("Pred_&_Res_Neg", "Pred_&_Res_Pos"), title = "Covariance between Predictor and Response Variables for LV 1")
                  saveWidget(LV.VIP[[i]][[j]]$shapes$lollipop.LV.1.Pred_and_Res, paste("Cov.LV.1", "_","Model", m, "_", names(Models[[m]])[i], "_", "sub", "_", j, ".html", sep = ""))

                  shape = plsCoVarCommonShape(placeholder, 2)
                  LV.VIP[[i]][[j]]$shapes$lollipop.LV.2.Pred_and_Res <- lolshape(shape.list = list(X_neg = array(shape[,,1], dim = c(dim(shape)[1], 3, 1)),
                                                                                                   X_pos = array(shape[,,2], dim = c(dim(shape)[1], 3, 1))),
                                                                                 ID = c("Predictor_Neg", "Response_Neg"), title = "Covariance between Predictor and Response Variables for LV 2")
                  saveWidget(LV.VIP[[i]][[j]]$shapes$lollipop.LV.2.Pred_and_Res, paste("Cov.LV.2", "_","Model", m, "_", names(Models[[m]])[i], "_", "sub", "_", j, ".html", sep = ""))

                } else {
                  if(length(dim(Response)) == 3) {

                    shape = plsCoVar(placeholder, i = 1, sdx = 2, sdy = 2)
                    LV.VIP[[i]][[j]]$shapes$lollipop.LV.1.Response <- lolshape(shape.list = list(Y_neg = array(shape$y[,,1], dim = c(dim(shape$y)[1],3,1)),
                                                                                                 Y_pos = array(shape$y[,,2], dim = c(dim(shape$y)[1],3,1))),
                                                                               ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Response Variable for LV 1")
                    shape = plsCoVar(placeholder, i = 2, sdx = 2, sdy = 2)
                    LV.VIP[[i]][[j]]$shapes$lollipop.LV.2.Response <- lolshape(shape.list = list(Y_neg = array(shape$y[,,1], dim = c(dim(shape$y)[1],3,1)),
                                                                                                 Y_pos = array(shape$y[,,2], dim = c(dim(shape$y)[1],3,1))),
                                                                               ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Response Variable for LV 2")
                  }
                  if(length(dim(Predictor)) == 3) {
                    shape = plsCoVar(placeholder, i = 1, sdx = 2, sdy = 2)
                    LV.VIP[[i]][[j]]$shapes$lollipop.LV.1.Predictor <- lolshape(shape.list = list(X_neg = array(shape$x[,,1], dim = c(dim(shape$x)[1],3,1)),
                                                                                                  X_pos = array(shape$x[,,2], dim = c(dim(shape$x)[1],3,1))),
                                                                                ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Predictor Variable for LV 1")
                    shape = plsCoVar(placeholder, i = 2, sdx = 2, sdy = 2)
                    LV.VIP[[i]][[j]]$shapes$lollipop.LV.2.Predictor <- lolshape(shape.list = list(X_neg = array(shape$x[,,1], dim = c(dim(shape$x)[1],3,1)),
                                                                                                  X_pos = array(shape$x[,,2], dim = c(dim(shape$x)[1],3,1))),
                                                                                ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Predictor Variable for LV 2")
                  }
                } # end of lollipop gen
              } else {
                LV.VIP[[i]][[j]]$summary <- data.frame(VIP_df = "No Data", VIP_RV = "No Data", VIP_p_value = "No Data")
              } # end of variable trimming
            } # end of j for loop paired subsets
          } # end of subsets
          # Unsubsetted data ___________________________________________________________________________________

          tryCatch({dt_parameters = modifyList(dt_parameters,Models[[m]][[i]][[3]])},
                   error = function(e){})

          if(length(dim(Models[[m]][[i]][[1]])) == 3){
            dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                           Predictor = Models[[m]][[i]][[2]],
                                                           Res_transform = "NULL")
                                        ) %>% .[intersect(names(.), dt_params)]


          }
          if(length(dim(Models[[m]][[i]][[2]])) == 3) {
            dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                           Predictor = Models[[m]][[i]][[2]],
                                                           Pred_transform = "NULL")
                                        ) %>% .[intersect(names(.), dt_params)]
          }
          if(length(dim(Models[[m]][[i]][[1]])) != 3 & length(dim(Models[[m]][[i]][[2]])) != 3) {
            dt_parameters = modifyList(dt_parameters, list(Response = Models[[m]][[i]][[1]],
                                                           Predictor = Models[[m]][[i]][[2]])
                                        ) %>% .[intersect(names(.), dt_params)]
          }

          Response = do.call(DT, dt_parameters)$Response
          Predictor = do.call(DT, dt_parameters)$Predictor

          if(length(dim(Response)) == 3 & length(dim(Predictor)) == 3) {
            config = TRUE
          } else {
            config = FALSE
          }

          placeholder = pls2B(y = Response, x = Predictor, same.config = config, cv = TRUE, rounds = 100)
          LV.all[[i]]$summary <- data.frame(df = ncol(Predictor), RV = placeholder$rv, p_value = placeholder$p.value.RV)

          if(isTRUE(config)) {

            shape = plsCoVarCommonShape(placeholder, 1)

            LV.all[[i]]$shapes$lollipop.LV.all.1.Pred_and_Res <- lolshape(shape.list = list(X_neg = array(shape[,,1], dim = c(dim(shape)[1], 3, 1)),
                                                                                            X_pos = array(shape[,,2], dim = c(dim(shape)[1], 3, 1))),
                                                                          ID = c("Pred_&_Res_Neg", "Pred_&_Res_Pos"), title = "Covariance between Predictor and Response Variables for LV.all 1")
            saveWidget(LV.all[[i]]$shapes$lollipop.LV.all.1.Pred_and_Res, paste("Cov.LV.all.1", "_","Model", m, "_", names(Models[[m]])[i], ".html", sep = ""))

            shape = plsCoVarCommonShape(placeholder, i = 2)
            LV.all[[i]]$shapes$lollipop.LV.all.2.Pred_and_Res <- lolshape(shape.list = list(X_neg = array(shape[,,1], dim = c(dim(shape)[1], 3, 1)),
                                                                                            X_pos = array(shape[,,2], dim = c(dim(shape)[1], 3, 1))),
                                                                          ID = c("Predictor_Neg", "Response_Neg"), title = "Covariance between Predictor and Response Variables for LV.all 2")
            saveWidget(LV.all[[i]]$shapes$lollipop.LV.all.2.Pred_and_Res, paste("Cov.LV.all.2", "_","Model", m, "_", names(Models[[m]])[i], ".html", sep = ""))

          } else {

            if(length(dim(Response)) == 3) {

              shape = plsCoVar(placeholder, i = 1, sdx = 2, sdy = 2)
              LV.all[[i]]$shapes$lollipop.LV.all.1.Response <- lolshape(shape.list = list(Y_neg = array(shape$y[,,1], dim = c(dim(shape$y)[1],3,1)),
                                                                                          Y_pos = array(shape$y[,,2], dim = c(dim(shape$y)[1],3,1))),
                                                                        ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Response Variable for LV.all 1")
              shape = plsCoVar(placeholder, i = 2, sdx = 2, sdy = 2)
              LV.all[[i]]$shapes$lollipop.LV.all.2.Response <- lolshape(shape.list = list(Y_neg = array(shape$y[,,1], dim = c(dim(shape$y)[1],3,1)),
                                                                                          Y_pos = array(shape$y[,,2], dim = c(dim(shape$y)[1],3,1))),
                                                                        ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Response Variable for LV.all 2")

            }
            if(length(dim(Predictor)) == 3) {
              shape = plsCoVar(placeholder, i = 1, sdx = 2, sdy = 2)
              LV.all[[i]]$shapes$lollipop.LV.all.1.Predictor <- lolshape(shape.list = list(X_neg = array(shape$x[,,1], dim = c(dim(shape$x)[1],3,1)),
                                                                                           X_pos = array(shape$x[,,2], dim = c(dim(shape$x)[1],3,1))),
                                                                         ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Predictor Variable for LV.all 1")
              shape = plsCoVar(placeholder, i = 2, sdx = 2, sdy = 2)
              LV.all[[i]]$shapes$lollipop.LV.all.2.Predictor <- lolshape(shape.list = list(X_neg = array(shape$x[,,1], dim = c(dim(shape$x)[1],3,1)),
                                                                                           X_pos = array(shape$x[,,2], dim = c(dim(shape$x)[1],3,1))),
                                                                         ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Predictor Variable for LV.all 2")

            }
          }
          #______________________________________________________________________________________________________

            LV.VIP.all[[i]]$vars <- colnames(data.frame(plsr.results.all[[i]]$vip))
            Predictor = data.frame(plsr.results.all[[i]]$vip)
            Predictor = if(ncol(Predictor) == 0) NULL else Predictor
          #if important variables were identified______________________________________________________________________________________________________
          if(!is.null(Predictor)){

            if(length(dim(Response)) == 3 & length(dim(Predictor)) == 3) {
              config = TRUE
            } else {
              config = FALSE
            }

            vars <- LV.VIP.all[[i]]$vars
            placeholder = pls2B(y = Response, x = Predictor, same.config = config, cv = TRUE, rounds = 100)
            LV.VIP.all[[i]]$summary <- data.frame(VIP_df = ncol(Predictor), VIP_RV = placeholder$rv, VIP_p_value = placeholder$p.value.RV)

            if(isTRUE(config)) {

              shape = plsCoVarCommonShape(placeholder, i = 1)

              LV.VIP.all[[i]]$shapes$lollipop.LV.1.Pred_and_Res <- lolshape(shape.list = list(X_neg = array(shape[,,1], dim = c(dim(shape)[1], 3, 1)),
                                                                                          X_pos = array(shape[,,2], dim = c(dim(shape)[1], 3, 1))),
                                                                        ID = c("Pred_&_Res_Neg", "Pred_&_Res_Pos"), title = "Covariance between Predictor and Response Variables for LV 1")
              saveWidget(LV.VIP.all[[i]]$shapes$lollipop.LV.1.Pred_and_Res, paste("Covariance.LV.1.Pred_and_Res", "_","Model", m, "_", "Cov", "_", i, ".html", sep = ""))

              shape = plsCoVarCommonShape(placeholder, i = 2)
              LV.VIP.all[[i]]$shapes$lollipop.LV.2.Pred_and_Res <- lolshape(shape.list = list(X_neg = array(shape[,,1], dim = c(dim(shape)[1], 3, 1)),
                                                                                          X_pos = array(shape[,,2], dim = c(dim(shape)[1], 3, 1))),
                                                                        ID = c("Predictor_Neg", "Response_Neg"), title = "Covariance between Predictor and Response Variables for LV 2")
              saveWidget(LV.VIP.all[[i]]$shapes$lollipop.LV.2.Pred_and_Res, paste("Covariance.LV.2.Pred_and_Res", "_","Model", m, "_", "Cov", "_", i, ".html", sep = ""))

            } else {

              if(length(dim(Response)) == 3) {

                shape = plsCoVar(placeholder, 1, sdx = 2, sdy = 2)
                LV.VIP.all[[i]]$shapes$lollipop.LV.1.Response <- lolshape(shape.list = list(Y_neg = array(shape$y[,,1], dim = c(dim(shape$y)[1],3,1)),
                                                                                        Y_pos = array(shape$y[,,2], dim = c(dim(shape$y)[1],3,1))),
                                                                      ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Response Variable for LV 1")
                shape = plsCoVar(placeholder, 2, sdx = 2, sdy = 2)
                LV.VIP.all[[i]]$shapes$lollipop.LV.2.Response <- lolshape(shape.list = list(Y_neg = array(shape$y[,,1], dim = c(dim(shape$y)[1],3,1)),
                                                                                        Y_pos = array(shape$y[,,2], dim = c(dim(shape$y)[1],3,1))),
                                                                      ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Response Variable for LV 2")

              }
              if(length(dim(Predictor)) == 3) {
                shape = plsCoVar(placeholder, 1, sdx = 2, sdy = 2)
                LV.VIP.all[[i]]$shapes$lollipop.LV.1.Predictor <- lolshape(shape.list = list(X_neg = array(shape$x[,,1], dim = c(dim(shape$x)[1],3,1)),
                                                                                         X_pos = array(shape$x[,,2], dim = c(dim(shape$x)[1],3,1))),
                                                                       ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Predictor Variable for LV 1")
                shape = plsCoVar(placeholder, 2, sdx = 2, sdy = 2)
                LV.VIP.all[[i]]$shapes$lollipop.LV.2.Predictor <- lolshape(shape.list = list(X_neg = array(shape$x[,,1], dim = c(dim(shape$x)[1],3,1)),
                                                                                         X_pos = array(shape$x[,,2], dim = c(dim(shape$x)[1],3,1))),
                                                                       ID = c("Neg_Cov", "Pos_Cov"), title = "Covariance among Predictor Variable for LV 2")
              }
            }
          } else {
            LV.VIP.all[[i]]$summary <- data.frame(df = "No Data", RV = "No Data", p_value = "No Data")
          }

        if(isTRUE(subset)) {
          return(list(LV.all = LV.all,
                      LV.VIP.all = LV.VIP.all,
                      LV = LV,
                      LV.VIP = LV.VIP))
        } else {
          return(list(LV.all = LV.all,
                      LV.VIP.all = LV.VIP.all))
        }
      } # end of i for parallel loop

      if(isTRUE(parallel)) {
        closeAllConnections()
      }

      if (isTRUE(print_progress)) {
        pb$tick(20/length(Models), tokens = list(what = paste("Model", m, "Finished Covariances:")))
      }


      if(length(Models[[m]]) == 1) {
        list2env(parallel.cov, envir = environment())
      } else {
        var_boot = if(isTRUE(subset)) c("LV.all", "LV.VIP.all", "LV", "LV.VIP") else (c("LV.all","LV.VIP.all"))
        booted.list = vector("list", length(var_boot))
        for(v in 1:length(var_boot)) {
          boot_temp = vector("list", nrow(parallel.cov))
          for(b in 1:nrow(parallel.cov)) {
            boot_temp[[b]] <- unlist(parallel.cov[b,v],recursive = FALSE)
            boot_temp[[b]] <- Filter(function(x) !is.null(x), boot_temp[[b]])
            boot_temp[[b]] <- unlist(boot_temp[[b]], recursive = FALSE)
            boot_temp[[b]] <- Filter(function(x) !is.null(x), boot_temp[[b]])
            boot_temp[[b]] = lapply(boot_temp[[b]], function(x) {
              x = Filter(function(y) !is.null(y), x)
              x
            })
            boot_temp[[b]] <- Filter(function(x) !purrr::is_empty(x), boot_temp[[b]])
          }
          booted.list[[v]] <- boot_temp
        }
      }
      names(booted.list) <- var_boot
      list2env(booted.list, envir = environment())
      rm(parallel.cov, booted.list, boot_temp)
      #______________________________________________________________________________________________________________________

      if (isTRUE(print_progress)) {
        pb$tick(5/length(Models), tokens = list(what = paste("Model", m, "Finalizing:")))
      }

      # Summarize Results______________________________________________________________________________________________________________________
      output[[m]]$Cov.all <- LV.all
      output[[m]]$Cov.all.VIP <- LV.VIP.all

      if(isTRUE(subset)){
        output[[m]]$Cov.subset <- LV
        output[[m]]$Cov.subset.VIP <- LV.VIP
      }

      names(output[[m]]$Cov.all) <- names(Models[[m]])
      names(output[[m]]$Cov.all.VIP) <- names(Models[[m]])
      if(isTRUE(subset)){
        for(i in seq_along(output[[m]]$Cov.subset)){
          for(j in seq_along(output[[m]]$Cov.subset[[i]])) {
            if("summary" %in% names(output[[m]]$Cov.subset[[i]])) {

            } else {
              names(output[[m]]$Cov.subset[[i]])[j] <- paste("Subset_", j, sep = "")
              names(output[[m]]$Cov.subset.VIP[[i]])[j] <- paste("Subset_", j, sep = "")
            }
          } # subsetted hypothesis test
        } # subsets
      } # hypothesis test
    } # end of if covs = TRUE statement

    output[[m]]$scores.all <- plots.vip.all
    output[[m]]$plsr.results.all <- plsr.results.all

    if(isTRUE(subset)){
      output[[m]]$scores <- plots.vip
      output[[m]]$plsr.results <- plsr.results
    }

    for(i in seq_along(output[[m]]$scores.all)){
      names(output[[m]]$scores.all)[i] <- names(Models[[m]])[i]
      names(output[[m]]$plsr.results.all)[i] <- names(Models[[m]])[i]
      if(isTRUE(subset)){
        names(output[[m]]$scores)[i] <- names(Models[[m]])[i]
        names(output[[m]]$plsr.results)[i] <- names(Models[[m]])[i]
        for(j in seq_along(output[[m]]$plsr.results[[i]])){
          if("summary" %in% names(output[[m]]$plsr.results[[i]])) {

          } else {
            names(output[[m]]$scores[[i]])[j] <- paste("Subset_", j, sep = "")
            names(output[[m]]$plsr.results[[i]])[j] <- paste("Subset_", j, sep = "")
          }
        } # subsetted hypothesis test
      } # subsets
    } # hypothesis test
  } # end of m for loop of each Model
  names(output) <- names(Models)

  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = paste("Finished! VIPs Analysis :D")))
  }

  class(output) <- "pls.bio"
  return(output)
} # end of function
