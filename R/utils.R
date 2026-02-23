#' Convert an ade4 PCA Object to Geomorph Format
#'
#' @title Convert ade4 Principal Component Analysis Object to Geomorph Object
#'
#' @description
#' This function converts a Principal Component Analysis (PCA) object from the ade4
#' package format to a format compatible with geomorph package functions. It extracts
#' and renames relevant components while maintaining the mathematical relationships
#' and statistical properties of the original PCA.
#'
#' @param dudi_pca An object of class dudi.pca from the ade4 package containing
#'   the results of a principal component analysis
#'
#' @details
#' The function performs the following conversions:
#' \itemize{
#'   \item Extracts eigenvalues and converts them to standard deviations
#'   \item Preserves the loadings matrix (variable contributions)
#'   \item Maintains the PCA scores for observations
#'   \item Retains centering and scaling information if available
#' }
#'
#' The function handles both scaled and unscaled PCA results, providing appropriate
#' defaults when centering or scaling information is not present in the original object.
#'
#' @return Returns a list of class "gm.prcomp" containing:
#' \itemize{
#'   \item sdev: Standard deviations of the principal components
#'   \item rotation: The variable loadings (eigenvectors)
#'   \item x: The rotated data (principal component scores)
#'   \item center: The centering vector used
#'   \item scale: The scaling vector used
#' }
#'
#' @importFrom ade4 dudi.pca
#' @importFrom stats prcomp
#'
#' @seealso
#' \code{\link[ade4]{dudi.pca}} for the original PCA implementation
#' \code{\link[geomorph]{gm.prcomp}} for the geomorph PCA implementation
#'
#' @author Keeling et al., 2025
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ade4)
#' # Create example data
#' data_matrix <- matrix(rnorm(100), ncol = 5)
#'
#' # Perform PCA using ade4
#' ade4_pca <- dudi.pca(data_matrix, scale = TRUE, scannf = FALSE)
#'
#' # Convert to geomorph format
#' gm_pca <- dudi_to_gm(ade4_pca)
#'
#' # Access components
#' print(head(gm_pca$x))  # PC scores
#' print(gm_pca$rotation) # Loadings
#' }
#'
#'
dudi_to_gm <- function(dudi_pca) {

  # Extract common elements from dudi.pca and rename them to match gm.prcomp structure
  output <- list(
    sdev = sqrt(dudi_pca$eig),            # Standard deviations of principal components
    rotation = as.matrix(dudi_pca$c1),   # Loadings (variable contributions to PCs)
    x = as.matrix(dudi_pca$li),          # PCA scores (coordinates of observations)
    center = if (!is.null(dudi_pca$cent)) dudi_pca$cent else rep(0, ncol(dudi_pca$tab)),  # Centering vector
    scale = if (!is.null(dudi_pca$norm)) dudi_pca$norm else rep(1, ncol(dudi_pca$tab))    # Scaling factors
  )

  # Assign a class to match gm.prcomp if needed
  class(output) <- "gm.prcomp"

  return(output)
}

#' Extract and Process Coefficients from MARS Models
#'
#' @title Extract Coefficients from Multivariate Adaptive Regression Splines
#'
#' @description
#' This function extracts and processes coefficients from Multivariate Adaptive
#' Regression Splines (MARS) models, specifically handling hinge functions and
#' aggregating coefficients for each predictor. It is particularly useful for
#' interpreting MARS models and understanding variable importance.
#'
#' @param coefficients A named vector of coefficients from a MARS model
#' @param Predictor_names A character vector containing the names of predictor
#'   variables used in the model
#'
#' @details
#' The function processes MARS model coefficients in several steps:
#' \itemize{
#'   \item Identifies and extracts base variable names from hinge functions
#'   \item Aggregates coefficients for each predictor across multiple terms
#'   \item Handles both linear and non-linear (hinge) terms
#'   \item Creates a standardized output format for easy interpretation
#' }
#'
#' The function is particularly useful for models created using the earth package
#' or similar MARS implementations, where coefficients may be spread across
#' multiple hinge functions for the same predictor.
#'
#' @return Returns a matrix with:
#' \itemize{
#'   \item Rows corresponding to predictor variables
#'   \item A single column containing the aggregated coefficients
#'   \item Row names matching the original predictor names
#' }
#'
#' @seealso
#' \code{\link[earth]{earth}}
#'
#' @importFrom stats coefficients
#'
#' @author Keeling et al., 2025
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create sample coefficients
#' coef_vector <- c(
#'   "(Intercept)" = 1.5,
#'   "h(x1-2)" = 0.5,
#'   "h(3-x1)" = -0.3,
#'   "x2" = 0.7
#' )
#'
#' # Define predictor names
#' pred_names <- c("x1", "x2")
#'
#' # Extract coefficients
#' result <- get_coef(coef_vector, pred_names)
#' print(result)
#' }
#'

get_coef <- function(coefficients, Predictor_names) {

  extract_name <- function(term) {
    if (grepl("h\\(", term)) {
      sub("h\\((.*)-.*\\)", "\\1", term)
    } else {
      term
    }
  }

  Output <- matrix(0, nrow = 1, ncol = length(Predictor_names))
  names(Output) <- Predictor_names

  for (term in names(coefficients)) {
    base_name <- extract_name(term)
    if (base_name %in% Predictor_names) {
      Output[base_name] <- Output[base_name] + coefficients[term]
    }
  }
  Output <- t(Output)
  return(Output)
}

#' Compute Confusion Matrices for Multivariate Factor Data
#'
#' @title Multivariate Confusion Matrix Analysis
#'
#' @description
#' Creates confusion matrices for multivariate categorical data, supporting both
#' matrices and data frames as input. Handles multiple categorical response
#' variables simultaneously and provides comprehensive confusion matrix statistics.
#'
#' @param predicted A matrix or data frame containing the predicted classifications.
#' @param observed A matrix or data frame containing the observed (true) classifications.
#'        Must have same dimensions as \code{predicted}.
#' @param by_variable Logical; if TRUE, computes separate confusion matrices for
#'        each variable. Default is FALSE.
#' @param var_names Character vector of variable names. If NULL, will use column
#'        names from input or generate sequential names.
#'
#' @return Returns an object of class "mv_confusion" containing:
#' \itemize{
#'   \item overall: Overall confusion matrix combining all variables
#'   \item by_variable: List of confusion matrices for each variable (if by_variable = TRUE)
#'   \item summary: List containing basic analysis information (number of variables,
#'         observations, variable names)
#' }
#'
#' @details
#' The function handles multiple categorical variables by:
#' \itemize{
#'   \item Automatically standardizing factor levels across predicted and observed values
#'   \item Creating combined states for overall confusion matrix analysis
#'   \item Supporting both factor and character inputs
#'   \item Providing comprehensive error handling and informative messages
#' }
#'
#' @seealso
#' \code{\link[caret]{confusionMatrix}} for the underlying confusion matrix computation
#'
#' \code{\link[stats]{factor}} for factor handling
#'
#' @author Keeling et al., 2025
#'
#' @keywords internal
#'
#' @references
#' Kuhn, M. (2008). Building Predictive Models in R Using the caret Package.
#' Journal of Statistical Software, 28(5), 1-26.
#'
#' Powers, D. M. W. (2011). Evaluation: From Precision, Recall and F-Factor to ROC,
#' Informedness, Markedness & Correlation. Journal of Machine Learning Technologies,
#' 2(1), 37-63.
#'
#' @importFrom caret confusionMatrix
#' @importFrom stats as.formula
#' @importFrom utils head tail
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example data
#' n <- 100
#'
#' # Generated sample data
#' observed_data <- data.frame(
#'   var1 = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
#'   var2 = factor(sample(c("0", "1"), n, replace = TRUE))
#' )
#'
#' # Create predicted data with some deliberate mistakes
#' predicted_data <- data.frame(
#'   var1 = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
#'   var2 = factor(sample(c("0", "1"), n, replace = TRUE))
#' )
#'
#' # Compute confusion matrices
#' results <- multivariate_confusion(
#'   predicted = predicted_data,
#'   observed = observed_data,
#'   by_variable = TRUE
#' )
#'
#' # Print overall results
#' print(results)
#'
#' # Access individual variable results
#' print(results$by_variable$var1)
#' print(results$by_variable$var2)
#'
#'}

multivariate_confusion <- function(predicted, observed, by_variable = FALSE, var_names = NULL) {
  # Convert inputs to data frames
  to_df <- function(x, names = NULL) {
    if (is.matrix(x)) {
      df <- as.data.frame(x)
      if (!is.null(names)) colnames(df) <- names
      return(df)
    } else if (is.data.frame(x)) {
      if (!is.null(names)) colnames(x) <- names
      return(x)
    } else {
      stop("Input must be either a matrix or data frame")
    }
  }

  # Validate dimensions and setup
  if (!identical(dim(predicted), dim(observed))) {
    stop("Predicted and observed must have same dimensions")
  }
  if (is.null(var_names)) {
    var_names <- colnames(observed) %||% colnames(predicted) %||%
      paste0("V", seq_len(ncol(predicted)))
  }

  # Convert to data frames
  pred_df <- to_df(predicted, var_names)
  obs_df <- to_df(observed, var_names)

  # Initialize results list
  results <- list()
  class_metrics <- list()

  # Process each column
  for(j in seq_along(var_names)) {
    col_name <- var_names[j]

    # Get unique classes from observed data
    obs_classes <- unique(as.character(obs_df[[j]]))

    # Special handling for binary True/False or 0/1 mapping
    if(all(obs_classes %in% c("True", "False")) && all(unique(pred_df[[j]]) %in% c(0, 1, "0", "1"))) {
      # Map 0/1 to False/True
      pred_mapped <- ifelse(pred_df[[j]] %in% c(0, "0"), "False", "True")
      obs_mapped <- obs_df[[j]]
    } else if(is.numeric(pred_df[[j]])) {
      # Map numeric to original class levels
      pred_mapped <- obs_classes[pred_df[[j]]]
      obs_mapped <- obs_df[[j]]
    } else {
      pred_mapped <- as.character(pred_df[[j]])
      obs_mapped <- as.character(obs_df[[j]])
    }

    # Calculate metrics for each class
    for(class in obs_classes) {
      pred_binary <- factor(pred_mapped == class, levels = c("FALSE", "TRUE"))
      obs_binary <- factor(obs_mapped == class, levels = c("FALSE", "TRUE"))

      conf <- confusionMatrix(pred_binary, obs_binary)
      class_metrics[[class]] <- conf$byClass
    }
  }

  # Combine all metrics into one matrix
  results <- do.call(rbind, class_metrics)

  class(results) <- c("mv_confusion", "list")
  return(results)
}

#' @title Reverse the contents of a CSGM data model
#'
#' @description This function takes a nested list object following the Models argument
#' format within the \pkg{CSGM} package and reverses the order of the response
#' and predictor variables.
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
#' @details This function takes either an AABV or PSLR object that has the argument
#' restructure.models set to FALSE (this argument automatically performs this function).
#' The function is ideal when subsetted tests of variables are performed. The function
#' will iterate over each regression model/hypothesis and each unique regression test.
#' If the regression test has subsetted variables that were each separately subject
#' to regression tests, the function will combine each subsetted test into a single
#' excel sheet within an excel file for that regression model/hypothesis.
#'
#' @returns the same nested list object but with the predictor and response
#' variables reversed.
#'
#' @seealso \code{restructure_models} \code{AABA} \code{pslr.bio}
#'
#' @importFrom pls plsr mvr
#' @importFrom plsdepot plsreg2
#'
#' @author Keeling et al., 2025
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#' # Reverse the model
#'
#' Models_revered = reverse.model(Models)
#' }
#'
#'

reverse.model <- function(Models) {
  reversed = list()
  reversed <- lapply(Models, rev)
  names(reversed) <- sapply(names(Models), function(name) {
    paste(rev(strsplit(name, "\\.")[[1]]), collapse = ".")
  })
  return(reversed)
}



#' @title Convert Principal Component Analysis Data to Original Data
#'
#' @description This function takes the eigenvalues, eigenvectors, and centered
#' principal component data to conduct a backwards reconstruction of the original
#' landmark configuration. This function is designed for reconstructing predicted
#' landmark configurations based on principal component data.
#'
#' @param Res_pca_scores a matrix or data table containing principal component scores (eigenvalues).
#' @param Res_rot a matrix or data table containing eigenvectors (rotation matrix) from a principal component analysis.
#' @param Res_pca_center a 2D matrix of a landmark configuration representing the center shape across all principal components.
#' @param Res_ncomp an optional numeric vector containing the number of principal components used in a regression.
#' @param Res_p The number of landmarks of the original landmark configuration.
#' @param Res_dim The number of dimensions of the original landmark configuration.
#'
#' @returns either 2D or 3D array of the predicted landmark configuration.
#'
#' @details This function estimates a landmark configuration based on the
#' eigenvalues, eigenvectors, and center landmark configuration from a Principal
#' Component Analysis.
#'
#' @author Keeling et al., 2025
#'
#' @importFrom geomorph arrayspecs
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#'  # Load in Corpus_Land test data
#'
#'  data(Corpus_Land)
#'
#'  # Record the number of landmarks and dimensions
#'
#'  Res_p = 600
#'  Res_dim = 3
#'
#'  # Conduct a Principal Component Analysis using the DT function
#'
#'  Response = DT(Response = Corpus_Land, Res_transform = "pca")
#'
#'  # Extract the Eigenvalues, Eigenvectors, and Center Shape
#'
#'  Res_pca_scores = Response$Res_pca_scores
#'  Res_rot = Response$Res_rot
#'  Res_pca_center = Response$Res_pca_center
#'
#'  # Conduct the reconstruction
#'
#'  output = pcland(Res_pca_scores,
#'                  Res_rot,
#'                  Res_pca_center,
#'                  Res_ncomp,
#'                  Res_p,
#'                  Res_dim)
#'
#' }
#'

pcland <- function(Res_pca_scores = NULL,
                   Res_rot = NULL,
                   Res_pca_center = NULL,
                   Res_ncomp = NULL,
                   Res_p = NULL,
                   Res_dim = NULL) {

  Res = as.data.frame(Res_pca_scores)

  if(is.null(Res_ncomp)) {

    names(Res) = paste(rep("Comp"), rep(1:ncol(Response)), sep = "")

  } else {

    names(Res) = paste(rep("Comp"), rep(Res_ncomp), sep = "")
    scores = data.frame(Res_pca_scores[,-c(1:ncol(Res))])
    pca = cbind(Res, scores)
  }

  pca = as.matrix.data.frame(Res)
  reconstruction <- pca %*% t(Res_rot)
  Prediction <- sweep(reconstruction, 2, Res_pca_center, "+")

  Prediction <- arrayspecs(Prediction, Res_p, Res_dim)
}



#' @title Trim VIP scores
#'
#' @description Trim Variable Importance in Projection (VIP) scores of predictor
#' variables by the mean, median, or scores greater than one.
#'
#' @param scores data table consisting of two variables: Response, Predictor, and Score.
#' This is an object generated by the \code{Bio.VIP} function
#'
#' @param lv_method Character string indicating the method to trim the predictor
#' variables by its VIP score. Options are "trim" (all predictor variables with
#' a score greater than 1), "median" (all predictor variables which are greater
#' than the median VIP score), and "mean" (all predictor variables which are
#' greater than the mean VIP score).
#'
#' @returns a vector containing the names of the VIP selected predictor variables.
#'
#' @details This function estimates important variable predictors using several
#' standard summary statistics such as the mean and median score as well as the
#' rule of thumb as being greater than one (the latent variable block variable
#' contribution exceeds the mean).
#'
#' @author Keeling et al., 2025
#'
#' @importFrom dplyr filter group_by summarise %>%
#'
#' @keywords internal
#'
#' @export

trim_scores <- function(scores,
                        lv_method) {

  if(lv_method == "trim") {
    scores <- scores %>% filter(Score >= 1)
    scores <- scores$Predictor
    if(is.null(scores) | is_empty(scores)) {
      lv_method = "median"
      attempted = TRUE
      print("No signficant Predictors were identified. Median method will be attempted.")
    }
  } else if(lv_method == "median") {
    scores <- scores %>% group_by(Predictor, Response) %>%
      summarise(median = median(Score), .groups = "drop") %>% filter(median >= 1)
    scores <- scores$Predictor
    if(is.null(scores | is_empty(scores))) {
      lv_method = "mean"
      print("No signficant Predictors were identified. Mean method will be attempted.")
    }
  } else if(lv_method == "mean") {
    scores <- scores %>% group_by(Predictor) %>%
      summarise(mean = mean(Score), .groups = "drop") %>% filter(mean >= 1)
    scores <- scores$Predictor
    if(is.null(scores | is_empty(scores))) {
      if(isTRUE(attempted)) {
        stop(print("No signficant Predictors were identified. Please choose alternative VIP method."))
      } else {
        print("No signficant Predictors were identified. Trim method will be attempted.")
      }
    }
  }
}

#' @title Calculate Summary Statistics from Partial Least Squares Models
#'
#' @description
#' A function designed to compute comprehensive performance statistics from partial least
#' squares (PLS) model objects through cross-validation procedures. This function serves
#' as an internal model assessment tool within Bio.VIP, calculating various measures of
#' fit and predictive ability including R2, Q2, and RMSEP. It handles different PLS model
#' types (spls, plsr2, mvr) and implements permutation testing for statistical significance.
#'
#' @param plsr.results A PLS model object from either mixOmics::spls(), plsdepot::plsreg2(),
#' or pls::mvr()
#' @param ncomp Integer specifying number of components to evaluate
#' @param CV_type Character string specifying type of cross-validation:
#' \itemize{
#'    \item{"CV"}{k-fold cross-validation with multiple repeats}
#'    \item{"LOO"}{Leave-one-out cross-validation}
#' }
#' @param repeats Integer specifying number of cross-validation repetitions
#' @param folds Integer specifying number of cross-validation folds
#' @param tuned Optional tuning parameters for sparse PLS models
#' @param parallel.me Logical value indicating whether to use parallel processing for
#' cross-validation
#' @param permut Integer specifying number of permutations for significance testing
#'
#' @returns A data frame containing model performance statistics:
#' \itemize{
#'   \item{X_var}{Explained variance in predictor variables}
#'   \item{R2}{Model fit statistic (squared correlation coefficient)}
#'   \item{Q2}{Cross-validated R2, indicates predictive ability}
#'   \item{RMSEP}{Root mean square error of prediction}
#'   \item{p_value_RMSE}{Statistical significance based on RMSEP permutation}
#'   \item{p_value_Q2}{Statistical significance based on Q2 permutation}
#' }
#'
#' @details
#' The plsr_stats function provides a unified approach to assessing PLS model performance
#' across different PLS implementations. For each model type, it computes both model fit
#' and predictive ability statistics through cross-validation procedures.
#'
#' The function first determines the PLS model type (spls, plsr2, or mvr) and adapts its
#' calculations accordingly. For sparse PLS (spls), it uses the mixOmics framework for
#' cross-validation. For plsr2 and mvr models, it implements custom cross-validation
#' procedures.
#'
#' Cross-validation is performed using either k-fold or leave-one-out methods. For k-fold
#' cross-validation, multiple repeats are performed to ensure stable estimates. During
#' cross-validation, the function:
#' 1. Splits data into training and test sets
#' 2. Fits PLS model on training data
#' 3. Predicts test set responses
#' 4. Calculates performance metrics
#'
#' The function calculates several key statistics:
#' - R2: Measures model fit as proportion of variance explained
#' - Q2: Cross-validated R2, indicates predictive ability
#' - RMSEP: Measures prediction error in original response units
#'
#' Statistical significance is assessed through permutation testing. The function:
#' 1. Randomly permutes response values
#' 2. Recalculates performance metrics
#' 3. Compares observed statistics to null distribution
#' 4. Calculates p-values
#'
#' The function implements adaptive convergence checking for permutation tests to optimize
#' computational efficiency while maintaining statistical rigor.
#'
#' @seealso \code{\link{Bio.VIP}} \code{\link{NComp}}
#'
#' @author Keeling et al., 2025
#'
#' @importFrom pls mvr plsr coefplot scores loadings explvar
#' @importFrom plsdepot plsreg2
#' @importFrom stats predict coef cor sd quantile pnorm p.adjust
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% registerDoSEQ
#'
#' @references
#' Rohart F, et al. (2017). "mixOmics: An R package for 'omics feature selection and multiple
#' data integration." PLoS Computational Biology, 13(11), e1005752.
#'
#' Wold, S., et al. (2001). PLS-regression: a basic tool of chemometrics. Chemometrics
#' and Intelligent Laboratory Systems, 58(2), 109-130.
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Calculate performance statistics for PLS analysis of symphyseal shape
#' # and biomechanical properties
#'
#' # Prepare data
#' Response <- Corpus_Land[241:360,,]  # Symphysis landmarks
#' Predictor <- Biomechanics %>%
#'   filter(Region == "LM1M2") %>%
#'   select(5:93)  # Biomechanical variables
#'
#' # Fit sparse PLS model
#' spls_model <- spls(Predictor, Response, ncomp = 5)
#'
#' # Calculate performance statistics
#' stats <- plsr_stats(
#'   plsr.results = spls_model,
#'   ncomp = 5,
#'   CV_type = "CV",
#'   repeats = 20,
#'   folds = 5,
#'   parallel.me = TRUE,
#'   permut = 500
#' )
#'
#' # View results
#' print(stats)
#' }

plsr_stats <- function(plsr.results, ncomp, vip, vip_method = "spls", cv_type = "CV", repeats = 25, folds = 5, tuned = NULL, parallel.me = TRUE, permut = 500) {

  if(class(plsr.results) == "mixo_spls") {
    explained_X_variance = max(plsr.results$prop_expl_var$X)
    y_score = data.frame(plsr.results$model[[1]])
    x_score = data.frame(plsr.results$X)

    nzv.X = (apply(x_score, 2, var) > .Machine$double.eps)
    y_pred2 = predict(plsr.results, newdata = x_score[,nzv.X])$predict

    R2 <- c_TR2(y_score,y_pred2,ncomp)

  } else {

    y_score <- data.frame(plsr.results$model[[1]])
    x_score <- data.frame(plsr.results$model[[2]])

    if(class(plsr.results) == "plsreg2") { # Multivariate PLS
      explained_X_variance <- plsr.results$expvar[ncomp, 2] * 100

      y_pred2 <- plsr.results$y.pred

      R2 <- c_TR2(y_score, y_pred2)

    } else if (class(plsr.results) == "mvr") {
      explained_X_variance <- sum(pls::explvar(plsr.results))
      y_pred2 <- plsr.results$fitted.values

      R2 <- c_TR2(y_score, y_pred2, ncomp)
      vip_method = plsr.results$method
    }
  } # if not spls which autogenerates values

  if(isTRUE(parallel.me)) {
    doParallel::registerDoParallel(detectCores() - 2)
  } else {
    foreach::registerDoSEQ()
  }

  CV_it <- foreach(cv = 1:repeats, .export = c("c_TR2", "c_R2","c_TQ2","c_RMSEP"), .packages = c("pls", "plsdepot", "caret")) %dopar% {
    rmsep <- vector("numeric", folds)
    Q2 <- vector("numeric", folds)
    Prediction <- array(NA, dim = c(nrow(y_score), ncol(y_score), ncomp))
    cv_folds <- caret::createFolds(y = 1:nrow(y_score), k = folds, list = TRUE, returnTrain = TRUE)
    for (fold_idx in 1:folds) {

      train_indices <- cv_folds[[fold_idx]]
      test_indices <- lubridate::setdiff(1:nrow(y_score), train_indices)

      X_test <- x_score[test_indices, ]
      Y_test <- y_score[test_indices, ]
      X_train <- x_score[train_indices, ]
      Y_train <- y_score[train_indices, ]

      if (class(plsr.results) == "plsreg2") {
        pls_model <- plsdepot::plsreg2(predictors = X_train, responses = Y_train, comps = ncomp, crosval = FALSE)
        pls_model$model <- list(Y = Y_train, X = X_train)
        pls_model <- pls_it(pls_model)
        Y_pred <- pls_model$Y_pred
        rmsep[fold_idx] <- mean(c_RMSEP(Y_train, Y_pred, ncomp = ncomp))
        test_Q2 <- max(c_TQ2(pls_model, Y_pred, ncomp = ncomp, Y_train))
        Q2[fold_idx] <- if(test_Q2 < 0 | test_Q2 > 1) max(c_TR2(Y_train, Y_pred, ncomp = ncomp)) else test_Q2
        for (comp in 1:ncomp) {
          Prediction[train_indices, ,comp] <- as.matrix(Y_pred[,,comp])
        }
      }
      if (class(plsr.results) == "mixo_spls") {
        nzv.X = (apply(X_train, 2, var) > .Machine$double.eps)
        pls_model <- spls(X_train, Y_train, ncomp = ncomp, mode = "regression", max.iter = 500, near.zero.var = TRUE)
        Y_pred <- predict(pls_model, newdata = X_test[,nzv.X])$predict
        pls_model$Y <- Y_train
        pls_model$type <- "mixo_spls"
        rmsep[fold_idx] <- mean(c_RMSEP(Y_test, Y_pred, ncomp))
        test_Q2 <- max(c_TQ2(pls_model, Y_pred, ncomp = ncomp, Y_test))
        Q2[fold_idx] <- if(test_Q2 < 0 | test_Q2 > 1) max(c_TR2(Y_test, Y_pred, ncomp = ncomp)) else test_Q2
        for (comp in 1:ncomp) {
          Prediction[test_indices,,comp] <- as.matrix(Y_pred[,,comp])
        }
      }
      if(class(plsr.results) == "mvr") {
        pls_model <- pls::mvr(as.matrix.data.frame(Y_train) ~ as.matrix.data.frame(X_train), ncomp = ncomp, method = vip_method, validation = "none", maxit = 500)
        pls_model$model <- list(Y = data.frame(Y_train), X = data.frame(X_train))
        Y_pred = predict(pls_model, newdata = as.matrix.data.frame(X_test), ncomp = 1:ncomp)
        pls_model <- pls_it(pls_model, X_test)
        rmsep[fold_idx] <- mean(c_RMSEP(Y_test, Y_pred, ncomp))
        test_Q2 <- max(c_TQ2(pls_model, Y_pred, ncomp = ncomp, Y_test))
        Q2[fold_idx] <- if(test_Q2 < 0 | test_Q2 > 1) max(c_TR2(Y_test, Y_pred, ncomp = ncomp)) else test_Q2
        for (comp in 1:ncomp) {
          Prediction[test_indices,,comp] <- as.matrix(Y_pred[,,comp])
        }
      }
    }
    Prediction <- apply(Prediction, c(1,2), mean)
    Prediction <- as.matrix(Prediction)
    colnames(Prediction) <- colnames(y_score)
    Q2 = max(Q2, na.rm = TRUE)
    rmsep = min(rmsep, na.rm = TRUE)


    return(list(Q2 = Q2, rmsep = rmsep, Prediction = Prediction))
  } # end of CV_it

  var_boot <- c("Q2", "rmsep", "Prediction")
    for (stat_name in var_boot) {
      assign(stat_name, lapply(CV_it, `[[`, stat_name))
    }


  if(isTRUE(parallel.me)) {
    closeAllConnections()
  }
  rm(CV_it)

  Prediction <- matrix(colMeans(do.call(rbind, Prediction), na.rm = TRUE), nrow = nrow(y_score), ncol = ncol(y_score))

  # Permutation tests
  rmsep_perm <- rep(NA, permut)
  Q2_perm <- rep(NA, permut)

  # Reduced model

  if(class(plsr.results) == "mixo_spls") {

    x_score_reduced = data.frame(x_score[,!colnames(x_score) %in% colnames(vip)])
    colnames(x_score_reduced) <- colnames(x_score)[!colnames(x_score) %in% colnames(vip)]

    if(ncol(x_score_reduced) == 0){
      keep = round(ncol(x_score) * .20,1)
      if(keep < 3){
        keep = 3
      }
      x_score_reduced = x_score[,c(1:keep)]

    } else if(ncol(x_score_reduced) <= 2){
      keep = 3 - ncol(x_score_reduced)
      keep = sample(colnames(x_score[,!colnames(x_score) %in% colnames(x_score_reduced)]),size = keep)
      keep = c(keep, colnames(x_score_reduced))
      x_score_reduced = x_score[,keep]
    }
    if(ncomp > ncol(x_score_reduced)) {
      ncomp2 = ncol(x_score_reduced)
    } else {
      ncomp2 = ncomp
    }

    x_score_reduced = as.matrix.data.frame(x_score_reduced)

    nzv.X = (apply(x_score_reduced, 2, var) > .Machine$double.eps)
    y_reduced = spls(x_score_reduced, y_score, ncomp = ncomp2, mode = "regression")
    y_reduced = predict(y_reduced, newdata = x_score_reduced[,nzv.X])$predict
    y_reduced = apply(y_reduced,c(1,2),mean)

  } else if (class(plsr.results) == "mvr"){

    x_score_reduced = data.frame(x_score[,!colnames(x_score) %in% colnames(vip)])
    colnames(x_score_reduced) <- colnames(x_score)[!colnames(x_score) %in% colnames(vip)]

    if(ncol(x_score_reduced) == 0){
      keep = round(ncol(x_score) * .20,1)
      if(keep < 3){
        keep = 3
      }
      x_score_reduced = x_score[,c(1:keep)]

    } else if(ncol(x_score_reduced) <= 2){
      keep = 3 - ncol(x_score_reduced)
      keep = sample(colnames(x_score[,!colnames(x_score) %in% colnames(x_score_reduced)]),size = keep)
      keep = c(keep, colnames(x_score_reduced))
      x_score_reduced = x_score[,keep]
    }
    if(ncomp > ncol(x_score_reduced)) {
      ncomp2 = ncol(x_score_reduced)
    } else {
      ncomp2 = ncomp
    }

    x_score_reduced = as.matrix.data.frame(x_score_reduced)

    y_reduced = pls::mvr(as.matrix.data.frame(y_score) ~ as.matrix.data.frame(x_score_reduced), ncomp = ncomp2, method = plsr.results$method, validation = "none", scale = FALSE)
    y_reduced = predict(y_reduced, newdata = as.matrix.data.frame(x_score_reduced), ncomp = 1:ncomp2)
    y_reduced = apply(y_reduced,c(1,2),mean)

  } else if (class(plsr.results) == "plsreg2") {
    x_score_reduced = data.frame(x_score[,!colnames(x_score) %in% colnames(vip)])
    colnames(x_score_reduced) <- colnames(x_score)[!colnames(x_score) %in% colnames(vip)]

    if(ncol(x_score_reduced) == 0){
      keep = round(ncol(x_score) * .20,1)
      if(keep < 3){
        keep = 3
      }
      x_score_reduced = x_score[,c(1:keep)]

    } else if(ncol(x_score_reduced) <= 2){
      keep = 3 - ncol(x_score_reduced)
      keep = sample(colnames(x_score[,!colnames(x_score) %in% colnames(x_score_reduced)]),size = keep)
      keep = c(keep, colnames(x_score_reduced))
      x_score_reduced = x_score[,keep]
    }
    if(ncomp > ncol(x_score_reduced)) {
      ncomp2 = ncol(x_score_reduced)
    } else {
      ncomp2 = ncomp
    }

    x_score_reduced = as.matrix.data.frame(x_score_reduced)

    y_reduced = plsdepot::plsreg2(predictors = x_score_reduced, responses = y_score, comps = ncomp2, crosval = FALSE)
    y_reduced = y_reduced$y.pred
  }

  convergence <- FALSE
  for (p in 1:permut) {

    if(isFALSE(convergence)) {

      samp_test <- sample(1:nrow(y_reduced), replace = FALSE)
      residuals_test <- (y_score - y_reduced)[samp_test,]
      rownames(residuals_test) <- 1:nrow(residuals_test)
      perm_y_score <- y_reduced + residuals_test
      perm_y_score <- as.matrix.data.frame(perm_y_score)

      if(class(plsr.results) == "mixo_spls") {
        nzv.X = (apply(plsr.results$X, 2, var) > .Machine$double.eps)
        pls_model <- spls(x_score, perm_y_score, ncomp = ncomp, mode = "regression", scale = FALSE)
        Y_pred <- predict(pls_model, newdata = plsr.results$X[,nzv.X])$predict
        pls_model$type <- "mixo_spls"
      }
      if(class(plsr.results) == "mvr") {
        pls_model <- pls::mvr(as.matrix.data.frame(perm_y_score) ~ as.matrix.data.frame(x_score), ncomp = ncomp, method = vip_method, validation = "none", scale = FALSE)
        pls_model$model <- list(Y = perm_y_score, X = x_score)
        Y_pred <- predict(pls_model, newdata = as.matrix.data.frame(x_score), ncomp = 1:ncomp)
        pls_model <- pls_it(pls_model, X_test, Y_test)
      }

      if(class(plsr.results) == "plsreg2") {
        pls_model <- plsdepot::plsreg2(predictors = x_score, responses = perm_y_score, comps = ncomp, crosval = FALSE)
        pls_model$model <- list(Y = perm_y_score, X = x_score)
        pls_model <- pls_it(pls_model)
        Y_pred <- pls_model$Y_pred

      }
      pls_model$Y <- perm_y_score
      test_Q2 <- max(c_TQ2(pls_model, Y_pred, ncomp = ncomp))
      Q2_perm[p] <- if(test_Q2 < 0 | test_Q2 > 1) max(c_TR2(perm_y_score, Y_pred, ncomp = ncomp)) else test_Q2
      rmsep_perm[p] <- mean(c_RMSEP(perm_y_score, Y_pred, ncomp))

      if(p >= 50) {
        se_rmsep <- sd(rmsep_perm[1:p]) / sqrt(p)
        # Check convergence
        if(se_rmsep <= 0.001) {
          convergence <- TRUE
          break
        }
      }
      if(isTRUE(convergence)) {
        break
      }
    }
    if(isTRUE(convergence)) {
      break
    }
  }

  rmsep_perm <- na.omit(rmsep_perm)
  Q2_perm <- na.omit(Q2_perm)

  Q2 <- mean(unlist(Q2), na.rm = TRUE)
  rmsep <- mean(unlist(rmsep), na.rm = TRUE)

  p_value_Q2 <- mean(p.adjust(Q2_perm >= Q2, method = "BH"), na.rm = TRUE)
  p_value_rmsep <- mean(p.adjust(rmsep_perm <= rmsep, method = "BH"), na.rm = TRUE)


  results <- data.frame(
    X_var = explained_X_variance,
    R2 = R2,
    Q2 = Q2,
    RMSEP = rmsep,
    p_RMSEP = p_value_rmsep,
    p_Q2 = p_value_Q2
  )
  return(results)
}

#' @title Transform partial least squares objects between package formats
#'
#' @description A standardization function that converts between different partial least squares (PLS)
#' object formats (pls, plsdepot, mixomics) to ensure consistent handling within Bio.VIP and AABA analyses.
#' This internal function enables seamless integration of different PLS implementations.
#'
#' @param plsr.results A PLS model object from either pls::plsr(), plsdepot::plsreg2(), or spls()
#' @param test_x Optional matrix/data.frame of test set predictor variables
#' @param test_y Optional matrix/data.frame of test set response variables
#'
#' @returns A standardized PLS object of class "mixo_spls" containing:
#' \itemize{
#'   \item{loadings:    List of X and Y variable loadings matrices}
#'   \item{variates:    List of X and Y score matrices}
#'   \item{ncomp:   Number of components}
#'   \item{Y_pred:    Array of predicted Y values}
#'   \item{mode:    Set to 'regression'}
#' }
#'
#' @details
#' This function standardizes PLS objects from different R packages into the mixOmics format
#' to enable consistent downstream processing in Bio.VIP analyses. It handles three types of input:
#' \itemize{
#'   \item{mvr objects from pls package}
#'   \item{plsr2 objects from plsdepot package}
#'   \item{mixo_spls objects from mixOmics package}
#' }
#' The function preserves all relevant model components while reformatting them to match
#' the mixOmics structure. This standardization is critical for Bio.VIP's ability to work
#' with multiple PLS implementations.
#'
#' @keywords internal
#'
#' @seealso \code{\link{Bio.VIP}} \code{\link{AABA}}
#'
#' @author Keeling et al., 2025
#'
#' @importFrom pls mvr plsr coefplot scores loadings
#' @importFrom plsdepot plsreg2
#' @importFrom stats predict coef cor sd quantile pnorm p.adjust
#'
#' @references
#' Rohart F., Gautier B., Singh A., Le Cao K.A. (2017). mixOmics: An R package for 'omics feature selection and multiple data integration. PLoS Comput Biol 13(11): e1005752
#'
#' Mevik, B.H., Wehrens, R. and Liland, K.H. (2019). pls: Partial Least Squares and Principal Component Regression. R package version 2.7-3.
#'
#' Sanchez, G. (2012). plsdepot: Partial Least Squares (PLS) Data Analysis Methods. R package version 0.1.17.
#'
#' @export

pls_it <- function(plsr.results, test_x = NULL, test_y = NULL) {

  object = list()

  if(is(plsr.results, "mvr")) {
    object$loadings <- list("X" = data.frame(as.matrix.data.frame(na.omit(plsr.results$loadings))), "Y" = data.frame(as.matrix.data.frame(na.omit(plsr.results$Yloadings))))
    object$variates <- list("X" = data.frame(as.matrix.data.frame(na.omit(plsr.results$scores))), "Y" = data.frame(as.matrix.data.frame(na.omit(plsr.results$Yscores))))
    object$keepX = ncol(plsr.results$model$X)
    object$keepY = ncol(plsr.results$model$Y)
    object$type <- "mvr"
    object$ncomp <- plsr.results$ncomp
  }
  if(is(plsr.results, "plsreg2")) {
    object$loadings <- list("X" = as.matrix(plsr.results$x.loads), "Y" = as.matrix(plsr.results$y.loads))
    object$variates <- list("X" = plsr.results$x.scores, "Y" = plsr.results$y.scores)
    object$keepX = ncol(plsr.results$model$X)
    object$keepY = ncol(plsr.results$model$Y)
    object$type <- "plsreg2"
    object$ncomp <- ncol(plsr.results$x.loads)
  }
  if(is(plsr.results, "mixo_spls")) {
    object <- plsr.results
    object$type <- "mixo_spls"
  }

  object$X <- data.frame(plsr.results$model[[2]])
  object$Y <- data.frame(plsr.results$model[[1]])

  if(object$type == "mvr") {
    #object$nzv.X = (apply(X, 2, var) > .Machine$double.eps)
    #object$nzv.Y = (apply(Y, 2, var) > .Machine$double.eps)
    newdata = if(is.null(test_x)) {as.matrix.data.frame(object$X)}else{as.matrix.data.frame(test_x)}
    object$Y_pred <- predict(plsr.results, newdata)
  }
  if(object$type == "plsreg2"){
    Y_pred = array(NA, dim = c(nrow(object$Y), ncol(object$Y), ncol(plsr.results$x.loads)))
    #object$nzv.X = (apply(X, 2, var) > .Machine$double.eps)
    #object$nzv.Y = (apply(Y, 2, var) > .Machine$double.eps)
    X = if(is.null(test_x)) {data.frame(object$X)}else{data.frame(test_x)}
    Y = if(is.null(test_y)) {data.frame(object$Y)}else{data.frame(test_y)}
    for(h in 1:ncol(plsr.results$x.loads)) {
      Y_pred[,,h] <- plsreg2(X, Y, comps = h)$y.pred
    }
    object$Y_pred <- Y_pred
  }

  object$mode <- 'regression'
  object$tol <- 1e-06
  object$max.iter = 1000
  object$scale <- FALSE

  class(object) <- "mixo_spls"

  return(object)
}


#' @title Calculate Residual Sum of Squares from PLS results
#'
#' @description An internal function that computes the Residual Sum of Squares (RSS)
#' from partial least squares regression results across multiple components. Used within
#' Bio.VIP to assess model fit quality.
#'
#' @param plsr.results A PLS model object containing component scores and loadings
#' @param ncomp Integer specifying number of components to evaluate
#'
#' @returns Matrix of RSS values with rows corresponding to components and columns to response variables
#'
#' @details
#' Calculates RSS by:
#' \itemize{
#'   \item{Extracting scores and loadings for each component}
#'   \item{Computing matrix deflation across dimensions}
#'   \item{Calculating sum of squared residuals}
#' }
#' The function handles both standard PLS and sparse PLS objects, adjusting calculations
#' accordingly. RSS values are used in downstream cross-validation statistics.
#'
#' @keywords internal
#'
#' @seealso \code{\link{Bio.VIP}} \code{\link{AABA}}
#'
#' @author Keeling et al., 2025
#'
#' @importFrom stats cor
#'
#' @references
#' Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic tool of chemometrics.
#' Chemometrics and intelligent laboratory systems, 58(2), 109-130.
#'
#' Lê Cao K-A, Rossouw D, Robert-Granié C, Besse P (2008). "A sparse PLS for
#' variable selection when integrating Omics data." Statistical Applications
#' in Genetics and Molecular Biology, 7(1), 35.
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A (2017). "mixOmics: An R package for
#' 'omics feature selection and multiple data integration." PLoS Computational
#' Biology, 13(11), e1005752.
#'
#' @export

RSS <- function(plsr.results, ncomp) {
  Output <- matrix(nrow = ncomp, ncol = ncol(plsr.results$Y))
  for (h in 1:ncomp) {
    if(class(plsr.results) == "mixo_spls") {
      tt = plsr.results$variates$X[, h]
      u = plsr.results$variates$Y[, h]
      b = plsr.results$loadings$Y[, h]
      #nx = p - keepX[h]
      #ny = q - keepY[h]
    } else if (plsr.results == "mvr"){
      tt = plsr.results$variates$X[, h]
      u = plsr.results$variates$Y[, h]
      b = plsr.results$loadings$Y[, h]
    } else if (plsr.results == "plsr2") {
      tt = plsr.results$variates$X[, h]
      u = plsr.results$variates$Y[, h]
      b = plsr.results$loadings$Y[, h]
    } else {
      Output[h,] <- sum((perm_test_y - colMeans(pred_test))^2)
    }
      # only used for matrices deflation across dimensions
      c = crossprod(as.matrix(plsr.results$X), tt)/drop(crossprod(tt))  #plsr.results$mat.c[, h]
      d = crossprod(as.matrix(plsr.results$Y), tt)/drop(crossprod(tt))  #plsr.results$mat.d[, h]
      e = crossprod(as.matrix(plsr.results$Y), u)/drop(crossprod(u))

      # deflate matrices
      Y = as.matrix(plsr.results$Y) - tt %*% t(d)

      Output[h, ] <- colSums((Y)^2)

  } # end h for loop
  return(Output)
}

#' @title Calculate Cross-validated Q2 Statistic from PLS results
#'
#' @description An internal function that computes the Q2 statistic (cross-validated R2)
#' for partial least squares models, providing a measure of predictive ability. Used within
#' Bio.VIP for model validation.
#'
#' @param plsr.results A PLS model object containing model components
#' @param Y_pred Matrix/array of predicted Y values
#' @param ncomp Integer specifying number of components to evaluate
#' @param Y_test optional n x p (obs, vars) matrix if using for test predictions from training model
#'
#' @returns Numeric vector of Q2 values for each component
#'
#' @details
#' Q2 statistic represents the cross-validated R2 value, calculated as:
#' 1 - (PRESS/RSS) where PRESS is the predicted residual sum of squares.
#' The function uses RSS calculations and predicted values to assess how well
#' the model predicts new observations. Q2 values above 0.5 typically indicate
#' good predictive ability.
#'
#' @keywords internal
#'
#' @seealso \code{\link{Bio.VIP}} \code{\link{AABA}}
#'
#' @importFrom stats cor
#'
#' @author Keeling et al., 2025
#'
#' @references
#' Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic tool of chemometrics.
#' Chemometrics and intelligent laboratory systems, 58(2), 109-130.
#'
#' Lê Cao K-A, Rossouw D, Robert-Granié C, Besse P (2008). "A sparse PLS for
#' variable selection when integrating Omics data." Statistical Applications
#' in Genetics and Molecular Biology, 7(1), 35.
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A (2017). "mixOmics: An R package for
#' 'omics feature selection and multiple data integration." PLoS Computational
#' Biology, 13(11), e1005752.
#'
#' @export

c_Q2 <- function(plsr.results, Y_pred, ncomp, Y_test = NULL) {
  RSS_info <- RSS(plsr.results, ncomp)
  Y_actual = if(!is.null(Y_test)) Y_test else plsr.results$Y
  PRESS.inside = matrix(nrow = ncomp, ncol = ncol(Y_actual))
  Q2 = numeric(ncomp)
  if(is.null(ncomp)) {
    Resid = as.matrix(Y_actual) - Y_pred
    PRESS.inside = sum(Resid^2)
    Q2 <- 1 - max(PRESS.inside / RSS_info)
  } else {
    for(h in 1:ncomp) {
      PRESS.inside[h, ] = sum((Y_actual - Y_pred[,,h])^2)
      Q2[h] <- 1 - max(PRESS.inside[h, ] / RSS_info[h, ])
    }
  }
  return(Q2)
}

#' @title Calculate Total Q2 Statistic from PLS results
#'
#' @description An internal function that computes the total Q2 statistic by selecting
#' either the maximum or mean Q2 value across components based on response variable type.
#' Used within Bio.VIP for overall model validation.
#'
#' @param plsr.results A PLS model object containing model components
#' @param Y_pred Matrix/array of predicted Y values
#' @param ncomp Integer specifying number of components to evaluate
#' @param Y_test optional n x p (obs, vars) matrix if using for test predictions from training model
#'
#' @returns Single numeric value representing total Q2 statistic
#'
#' @details
#' Extends Q2 calculation by providing a single summary statistic:
#' \itemize{
#'   \item{For PCA-type responses (identified by "Comp" in column names), uses maximum Q2}
#'   \item{For other responses, uses mean Q2 across components}
#' }
#' This adaptation allows appropriate handling of different response variable types
#' in multivariate analyses.
#'
#' @keywords internal
#'
#' @seealso \code{\link{Bio.VIP}} \code{\link{AABA}}
#'
#' @importFrom stats cor
#'
#' @author Keeling et al., 2025
#'
#' @references
#' Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic tool of chemometrics.
#' Chemometrics and intelligent laboratory systems, 58(2), 109-130.
#'
#' Lê Cao K-A, Rossouw D, Robert-Granié C, Besse P (2008). "A sparse PLS for
#' variable selection when integrating Omics data." Statistical Applications
#' in Genetics and Molecular Biology, 7(1), 35.
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A (2017). "mixOmics: An R package for
#' 'omics feature selection and multiple data integration." PLoS Computational
#' Biology, 13(11), e1005752.
#'
#' @export

c_TQ2 <- function(plsr.results, Y_pred, ncomp, Y_test = NULL) {
  RSS_info <- RSS(plsr.results, ncomp)
  Y_actual = if(!is.null(Y_test)) Y_test else plsr.results$Y
  PRESS.inside = matrix(nrow = ncomp, ncol = ncol(Y_actual))
  Q2 = numeric(ncomp)
  if(is.null(ncomp)) {
    Resid = as.matrix(Y_actual) - Y_pred
    PRESS.inside = sum(Resid^2)
    Q2 <- 1 - max(PRESS.inside / RSS_info)
  } else {
    for(h in 1:ncomp) {
      PRESS.inside[h, ] = sum((Y_actual - Y_pred[,,h])^2)
      Q2[h] <- 1 - max(PRESS.inside[h, ] / RSS_info[h, ])
    }
  }
  # Decide to use max or mean Q2 based on the response variable type
  if (any(grepl("Comp", colnames(plsr.results$Y)))) {
    TQ2 <- max(Q2, na.rm = TRUE)
  } else {
    TQ2 <- mean(Q2, na.rm = TRUE)
  }
  return(TQ2)
}


#' @title Calculate R-square Statistic from PLS results
#'
#' @description An internal function that computes R2 values for each component in a
#' partial least squares model. Used within Bio.VIP for assessing model fit quality.
#'
#' @param Y_actual Matrix of observed response values
#' @param Y_pred Matrix/array of predicted response values
#' @param ncomp Integer specifying number of components to evaluate
#'
#' @returns Numeric vector of R2 values for each component
#'
#' @details
#' Calculates R2 as squared correlation between actual and predicted values.
#' For multiple response variables, computes diagonal correlation matrix.
#' When no ncomp specified, returns single R2 value; otherwise returns R2
#' for each component up to ncomp.
#'
#' @seealso \code{\link{Bio.VIP}} \code{\link{AABA}}
#'
#' @importFrom stats cor
#'
#' @author Keeling et al., 2025
#'
#' @references
#' Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic tool of chemometrics.
#' Chemometrics and intelligent laboratory systems, 58(2), 109-130.
#'
#' Lê Cao K-A, Rossouw D, Robert-Granié C, Besse P (2008). "A sparse PLS for
#' variable selection when integrating Omics data." Statistical Applications
#' in Genetics and Molecular Biology, 7(1), 35.
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A (2017). "mixOmics: An R package for
#' 'omics feature selection and multiple data integration." PLoS Computational
#' Biology, 13(11), e1005752.
#'
#' @export

c_R2 <- function(Y_actual, Y_pred, ncomp = NULL) {
  Y_actual <- as.matrix(Y_actual)
  R2 = as.numeric(ncomp)
  if(is.null(ncomp)) {
    R2 = diag(cor(Y_actual, Y_pred))^2
  } else {
    for (h in 1:ncomp) {
      if(any(grepl("Comp", colnames(Y_actual)))) {
        R2[h] = max(diag(cor(Y_actual, Y_pred[, , h]))^2)
      } else {
        R2[h] = mean(diag(cor(Y_actual, Y_pred[, , h]))^2)
      }
    }
  }
  return(as.numeric(R2))
}

#' @title Multivariate Coefficient of Determination (R^2) Calculation
#'
#' @description Calculates appropriate R2 based on data characteristics:
#' 1. Correlation-based R2 for PCA/ranked data
#' 2. Regularized trace-based R2 with Ledoit-Wolf shrinkage
#' 3. Pillai's Trace-based R2 for multivariate data
#' 4. Determinant based approach
#' 5. Pseudo (Eigen) determinant R2 for multivariate data
#' 6. Averaged R2
#'
#' @param observed Matrix of observed values
#' @param predicted Matrix of predicted values
#'
#' @return a numeric value containing the coefficient of determination
#'
#' @details A function to calculate the coefficient of determination for multivariate
#' data. This function conducts all six methodologies for different multivariate
#' data contexts, normalizes negative R2 (poor fits for the multivariate data) and
#' finds the maximum R2 that is between 0 and 1. If no result is found (poor model
#' performance), 0 is selected.
#'
#' @keywords internal
#'
#' @export

m_R2 <- function(observed, predicted) {
  observed <- as.matrix(observed)
  predicted <- as.matrix(predicted)
  residuals <- observed - predicted
  # Get dimensions
  n <- nrow(observed)
  p <- ncol(observed)
  lambda = 1e-6

  # Method 1: Correlation-based R² for ranked/PCA data (Squared Pearson correlation Coefficient)
    r2_test1 <- max(diag(cor(observed, predicted))^2)

  # Method 2: "Regularized trace R²" or "Shrinkage-based multivariate R2"
    shrinkage <- min(0.5, p/n)  # Ledoit-Wolf type shrinkage parameter

    diag_mean_res <- mean(diag(var(residuals)))
    diag_mean_obs <- mean(diag(var(observed)))

    cov_residuals <- (1-shrinkage) * var(residuals) + shrinkage * diag(diag_mean_res, p)
    cov_observed <- (1-shrinkage) * var(observed) + shrinkage * diag(diag_mean_obs, p)

    r2_test2 <- 1 - (sum(diag(cov_residuals)) / sum(diag(cov_observed)))

  # Method 3: "Pillai's trace criterion" or "Multivariate trace R²"
    cov_residuals <- cov(residuals) + diag(lambda, p)
    cov_observed <- cov(observed) + diag(lambda, p)
    r2_test3 <- 1 - (sum(diag(cov_residuals)) / sum(diag(cov_observed)))

  # Method 4: Determinant R2, Multivariate determinant R2
    cov_residuals <- cov(residuals)
    cov_observed <- cov(observed)
    cov_residuals <- cov_residuals + diag(lambda, p)
    cov_observed <- cov_observed + diag(lambda, p)

    # Calculate R² using determinants
    r2_test4 <- 1 - (det(cov_residuals) / det(cov_observed))

  # Method 5: Pseudo determinant R2 with eigenvalues
    cov_residuals <- cov(residuals)
    cov_observed <- cov(observed)

    # Calculate eigenvalues
    eig_residuals <- eigen(cov_residuals, only.values = TRUE)$values
    eig_observed <- eigen(cov_observed, only.values = TRUE)$values

    # Filter out near-zero eigenvalues (adjust tolerance as needed)
    tol <- 1e-10
    eig_residuals <- eig_residuals[eig_residuals > tol]
    eig_observed <- eig_observed[eig_observed > tol]

    # Use product of non-zero eigenvalues as pseudodeterminant
    pseudo_det_residuals <- prod(eig_residuals)
    pseudo_det_observed <- prod(eig_observed)

    # Calculate R² using pseudodeterminants
    r2_test5 <- 1 - (pseudo_det_residuals / pseudo_det_observed)

  # Method 6: Averaged squared correlation coefficient
    r2_values <- numeric(ncol(observed))

    for(i in 1:ncol(observed)) {
      SST <- sum((observed[,i] - mean(observed[,i]))^2)
      SSE <- sum((observed[,i] - predicted[,i])^2)
      r2_values[i] <- 1 - (SSE/SST)
    }
    r2_test6 = mean(r2_values, na.rm = TRUE)

  r2 = numeric(6)
  r2 = na.omit(c(r2_test1,r2_test2,r2_test3,r2_test4,r2_test5,r2_test6))

  r2 = pmax(0, pmin(1, r2))

  r2 = data.frame(R2 = r2) %>% filter(R2 < 1 & R2 > 0)
  r2 = if(nrow(r2) == 0) 0 else max(r2)

  return(r2)
}

#' @title Calculate Total R-square Statistic from PLS results
#'
#' @description An internal function that computes a single summary R2 statistic from
#' partial least squares results, using either maximum or mean R2 based on response type.
#' Used within Bio.VIP for overall model assessment.
#'
#' @param Y_actual Matrix of observed response values
#' @param Y_pred Matrix/array of predicted response values
#' @param ncomp Optional integer specifying number of components to evaluate
#'
#' @returns Single numeric value representing total R2 statistic
#'
#' @details
#' Similar to total Q2 calculation, provides single summary R2:
#' \itemize{
#'   \item{For PCA-type responses (identified by "Comp" in column names), uses maximum R2}
#'   \item{For other responses, uses mean R2 across components}
#' }
#' This adaptation ensures appropriate summarization for different response variable types.
#'
#' @keywords internal
#'
#' @seealso \code{\link{Bio.VIP}} \code{\link{AABA}}
#'
#' @importFrom stats cor
#'
#' @author Keeling et al., 2025
#'
#' @references Mevik, B. H., & Wehrens, R. (2007). The pls package: principal component and partial
#' least squares regression in R. Journal of Statistical Software, 18(2), 1-24.
#'
#' @export

c_TR2 <- function(Y_actual, Y_pred, ncomp = NULL) {
  Y_actual <- as.matrix(Y_actual)
  R2 = as.numeric(ncomp)
  if(is.null(ncomp)) {
    R2 = diag(cor(Y_actual, Y_pred))^2
  } else {
    for (h in 1:ncomp) {
      if(any(grepl("Comp", colnames(Y_actual)))) {
        R2[h] = max(diag(cor(Y_actual, Y_pred[, , h]))^2)
      } else {
        R2[h] = mean(diag(cor(Y_actual, Y_pred[, , h]))^2)
      }
    }
  }
  if(any(grepl("Comp", colnames(Y_actual)))) {
    R2 = max(as.numeric(R2))
  } else {
    R2 = mean(as.numeric(R2))
  }
  return(R2)
}


#' @title Calculate RMSEP Statistic from PLS results
#'
#' @description An internal function that computes the Root Mean Square Error of
#' Prediction (RMSEP) for partial least squares models. Used within Bio.VIP for
#' assessing prediction accuracy.
#'
#' @param Y_actual Matrix of observed response values
#' @param Y_pred Matrix/array of predicted response values
#' @param ncomp Optional integer specifying number of components to evaluate
#'
#' @returns Numeric vector of RMSEP values, either single value or one per component
#'
#' @details
#' Calculates RMSEP as square root of mean squared prediction error.
#' When ncomp is specified, computes RMSEP for each component up to ncomp.
#' RMSEP provides measure of prediction accuracy in original response variable units.
#'
#' @keywords internal
#'
#' @seealso \code{\link{Bio.VIP}} \code{\link{AABA}}
#'
#' @importFrom stats cor
#'
#' @author Keeling et al., 2025
#'
#' @references
#' Mevik, B. H., & Wehrens, R. (2007). The pls package: principal component and partial
#' least squares regression in R. Journal of Statistical Software, 18(2), 1-24.
#'
#' @export

c_RMSEP <- function(Y_actual, Y_pred, ncomp = NULL) {
  if(is.null(ncomp)) {
    Y_actual <- as.matrix(Y_actual)

    RMSEP <- numeric(ncomp)
    RMSEP <- sqrt(mean((Y_actual - Y_pred)^2))
  } else {

    Y_actual <- as.matrix(Y_actual)

    RMSEP <- vector("numeric",ncomp)

    for (h in 1:ncomp) {
      RMSEP[h] <- sqrt(mean((Y_actual - Y_pred[,,h])^2))
    }
  }
  return(as.numeric(RMSEP))
}

#' @title NComp: Calculate best number of latent variables for a Partial Least Squares Analysis
#'
#' @description
#' A function designed to identify the optimal number of latent variables (components)
#' in Partial Least Squares analysis by evaluating model performance through cross-validation.
#' This function supports multiple PLS methods and implements both one-sigma and
#' randomization approaches for component selection. It optimizes component selection
#' by maximizing predictive power (Q2) while minimizing prediction error (RMSEP),
#' making it particularly useful for high-dimensional biological data analysis.
#'
#' @param model List object containing response and predictor data in a model format:
#' model$model((1)) = Response data
#' model$model((2)) = Predictor data
#' @param method Character string specifying component selection method:
#' \itemize{
#'   \item{"onesigma":    Selects simplest model within one standard error of minimum RMSEP}
#'   \item{"randomization":   Uses permutation testing to assess component significance}
#' }
#' @param vip_method Character string specifying PLS method:
#' \itemize{
#'   \item{"plsr2":   Two-block PLS via plsdepot}
#'   \item{"pls":   Standard PLS via pls package}
#'   \item{"spls":    Sparse PLS via mixOmics}
#' }
#' @param cv_type Character string specifying cross-validation type:
#' \itemize{
#'   \item{"CV":    k-fold cross-validation}
#'   \item{"LOO":   Leave-one-out cross-validation}
#' }
#' @param scale Logical value indicating whether to scale variables
#' @param center Logical value indicating whether to center variables
#' @param folds Integer specifying number of cross-validation folds
#' @param max_comps Integer specifying maximum number of components to evaluate
#' @param repeats Integer specifying number of cross-validation repetitions
#' @param nperm Integer specifying number of permutations for randomization test
#' @param alpha Numeric value specifying significance level for randomization test
#' @param plot Logical value indicating whether to generate diagnostic plots
#' @param save.plot Logical value indicating whether to save diagnostic plots
#'
#' @returns A list containing:
#' \itemize{
#'   \item{ncomp:   Optimal number of components}
#'   \item{plot:    If plot=TRUE, diagnostic plot showing RMSEP vs components}
#' }
#'
#' @details
#' NComp implements a systematic approach to determining the optimal number of
#' components in PLS analysis. For the one-sigma method (method="onesigma"), it
#' identifies the simplest model whose Root Mean Square Error of Prediction (RMSEP)
#' is within one standard error of the minimum RMSEP. This approach balances model
#' complexity against predictive power.
#'
#' For the randomization method (method="randomization"), the function performs
#' permutation testing to assess the significance of each additional component.
#' Components are retained until their contribution becomes non-significant at the
#' specified alpha level.
#'
#' Cross-validation is implemented through either k-fold (cv_type="CV") or
#' leave-one-out (cv_type="LOO") approaches. For k-fold cross-validation, the
#' function performs multiple repetitions to ensure stable estimates. Model
#' performance is assessed through both Q2 (predictive power) and RMSEP
#' (prediction error) statistics.
#'
#' The function adapts its approach based on the specified vip_method. For plsr2,
#' it handles two-block PLS specifically. For spls, it incorporates variable
#' selection considerations. The function also supports standard PLS implementations
#' through the pls package.
#'
#' When plot=TRUE, the function generates diagnostic plots showing RMSEP across
#' different numbers of components, with the selected optimal number clearly marked.
#'
#' @seealso \code{\link{Bio.VIP}} \code{\link{AABA}}
#'
#' @author Keeling et al., 2025
#'
#' @importFrom caret createFolds trainControl
#' @importFrom stats runif
#' @importFrom pls mvr
#' @importFrom plsdepot plsreg2
#'
#' @references
#' Wold, S., Johansson, E., & Cocchi, M. (1993). PLS: Partial Least Squares Projections
#' to Latent Structures. 3D QSAR in Drug Design, 1, 523-550.
#'
#' Mevik, B. H., & Wehrens, R. (2007). The pls Package: Principal Component and
#' Partial Least Squares Regression in R. Journal of Statistical Software, 18(2), 1-24.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Determine optimal number of components for analyzing relationship
#' # between symphyseal shape and biomechanical properties
#'
#' # Prepare data
#' Response <- Corpus_Land[241:360,,]  # Symphysis landmarks
#' Predictor <- Biomechanics %>%
#'   filter(Region == "LM1M2") %>%
#'   select(5:93)  # Biomechanical variables
#'
#' # Create model structure
#' model <- list(
#'   model = list(
#'     Response,
#'     Predictor
#'   )
#' )
#'
#' # Find optimal number of components using one-sigma method
#' ncomp_onesig <- NComp(
#'   model = model,
#'   method = "onesigma",
#'   vip_method = "plsr2",
#'   cv_type = "CV",
#'   folds = 5,
#'   repeats = 20,
#'   plot = TRUE
#' )
#'
#' # Find optimal number of components using randomization
#' ncomp_rand <- NComp(
#'   model = model,
#'   method = "randomization",
#'   vip_method = "plsr2",
#'   cv_type = "CV",
#'   folds = 5,
#'   repeats = 20,
#'   nperm = 500,
#'   plot = TRUE
#' )
#'
#' # Compare results
#' print(paste("One-sigma components:", ncomp_onesig$ncomp))
#' print(paste("Randomization components:", ncomp_rand$ncomp))
#' }

NComp <- function(model,
                  method = "onesigma",
                  vip_method = "spls",
                  cv_type = "CV",
                  scale = FALSE,
                  center = FALSE,
                  folds = 5,
                  max_comps = 10,
                  repeats = 25,
                  nperm = 1000,
                  alpha = 0.01,
                  plot = FALSE,
                  save.plot = FALSE,
                  parallel.me = FALSE) {

  randomiz.test <- function(residualsNew, residualsReference, nperm) {
    d <- residualsNew^2 - residualsReference^2
    md <- mean(d)
    N <- length(d)

    signs <- round(matrix(runif(N * nperm), N, nperm)) * 2 - 1
    dsigns <- d * signs
    mdsigns <- colMeans(dsigns)

    count <- sum(mdsigns >= md)
    (count + .5) / (nperm + 1)
  }

  X <- as.matrix.data.frame(model[[2]])
  Y <- as.matrix.data.frame(model[[1]])

  max_comps = if(ncol(Y) < max_comps || ncol(X) < max_comps) min(c(ncol(Y),ncol(X))) else max_comps

  if(isTRUE(parallel.me)) {
    doParallel::registerDoParallel(detectCores() - 2)
  } else {
    foreach::registerDoSEQ()
  }

  CV_it <- foreach(cv = 1:repeats, .export = c("c_TR2", "c_R2","c_TQ2","c_RMSEP"), .packages = c("pls", "plsdepot", "caret")) %dopar% {

     rmsep <- matrix(NA, nrow = folds, ncol = max_comps)
     Q2 <- matrix(NA, nrow = folds, ncol = max_comps)
     if (cv_type == "LOO") {
       n <- nrow(X)
       folds <- n
       indices <- 1:n
     } else if (cv_type == "CV") {
       cv_folds <- createFolds(y = 1:nrow(Y), k = folds, list = TRUE, returnTrain = TRUE)
     }

     for (fold_idx in 1:folds) {
       if(cv_type == "CV"){
         train_indices <- cv_folds[[fold_idx]]
         test_indices <- setdiff(1:nrow(X), train_indices)
       } else if (cv_type == "LOO"){
         train_indices <- indices[-i]
         test_indices <- indices[i]
       }

       X_train <- X[train_indices, ]
       Y_train <- Y[train_indices, ]
       X_test <- X[test_indices, ]
       Y_test <- Y[test_indices, ]

       if(vip_method == "plsr2") {
         pls_model <- plsdepot::plsreg2(predictors = X_train, responses = Y_train, comps = max_comps, crosval = FALSE)
         pls_model$model <- list(Y = Y_train, X = X_train)
         pls_model <- pls_it(pls_model)
         Y_pred <- pls_model$Y_pred
         rmsep[fold_idx,] <- c_RMSEP(Y_train, Y_pred, ncomp = max_comps)
         test_Q2 <- c_Q2(pls_model, Y_pred, ncomp = max_comps)
         Q2[fold_idx,] <- if(median(test_Q2) < 0 | median(test_Q2) > 1) c_R2(Y_train, Y_pred, ncomp = max_comps) else test_Q2
       } else if(vip_method == "spls") {
         nzv.X = (apply(X_test, 2, var) > .Machine$double.eps)
         pls_model <- spls(X_train, Y_train, ncomp = max_comps, mode = "regression", max.iter = 500, near.zero.var = TRUE)
         Y_pred <- predict(pls_model, newdata = X_test[,nzv.X])$predict
         pls_model$Y <- Y_train
         pls_model$type <- "mixo_spls"
         rmsep[fold_idx,] <- c_RMSEP(Y_test, Y_pred, ncomp = max_comps)
         test_Q2 <- c_Q2(pls_model, Y_pred, ncomp = max_comps, Y_test)
         Q2[fold_idx,] <- if(median(test_Q2) < 0 | median(test_Q2) > 1) c_R2(Y_test, Y_pred, ncomp = max_comps) else test_Q2
       } else {
         pls_model <- pls::mvr(as.matrix.data.frame(Y_train) ~ as.matrix.data.frame(X_train), ncomp = max_comps, method = vip_method, validation = "none", maxit = 500)
         pls_model$model <- list(Y = data.frame(Y_train), X = data.frame(X_train))
         Y_pred = predict(pls_model, newdata = as.matrix.data.frame(X_test), ncomp = 1:max_comps)
         pls_model <- pls_it(pls_model, X_test)
         rmsep[fold_idx,] <- c_RMSEP(Y_test, Y_pred, ncomp = max_comps)
         test_Q2 <- c_Q2(pls_model, Y_pred, ncomp = max_comps, Y_test)
         Q2[fold_idx,] <- if(median(test_Q2) < 0 | median(test_Q2) > 1) c_R2(Y_test, Y_pred, ncomp = max_comps) else test_Q2
       }
     }
     return(list(Q2 = Q2, rmsep = rmsep))
  } # end of CV_it

  var_boot <- c("Q2", "rmsep")
  for (stat_name in var_boot) {
    assign(stat_name, lapply(CV_it, `[[`, stat_name))
  }

  if(isTRUE(parallel.me)) {
    closeAllConnections()
  }
  rm(CV_it)

  Q2 <- do.call(rbind, lapply(Q2,function(x) colMeans(x)))
  rmsep <- do.call(rbind, lapply(rmsep,function(x) colMeans(x)))

  mean_rmsep <- colMeans(rmsep, na.rm = TRUE)
  se_rmsep <- apply(rmsep, 2, sd, na.rm = TRUE) / sqrt(repeats)

  mean_q2 <- colMeans(Q2, na.rm = TRUE)
  se_q2 <- apply(Q2, 2, sd, na.rm = TRUE) / sqrt(repeats)

  # Optimize based on highest Q² and lowest RMSEP
  if (method == "onesigma") {
    # One sigma method: choose the simplest model within one standard error of the minimum RMSEP
    min_rmsep <- min(mean_rmsep)
    one_sigma_threshold <- min_rmsep + se_rmsep[which.min(mean_rmsep)]
    optimal_ncomps_rmsep <- which(mean_rmsep <= one_sigma_threshold)[1]

    # Choose the model with the highest cumulative Q²
    optimal_ncomps_q2 <- which.max(mean_q2)

    # Combine criteria: prioritize Q², but ensure RMSEP is within acceptable range
    optimal_ncomps <- ifelse(mean_rmsep[optimal_ncomps_q2] <= one_sigma_threshold, optimal_ncomps_q2, optimal_ncomps_rmsep)
  } else {
    # Randomization test method
    pvals <- sapply(seq_len(max_comps), function(ncomp) {
      randomiz.test(rmsep_values[, ncomp], rmsep_values[, which.min(mean_rmsep)], nperm = nperm)
    })
    idx <- which(pvals > alpha)
    optimal_ncomps_rmsep <- as.numeric(min(c(idx, which.min(mean_rmsep))))

    # Choose the model with the highest cumulative Q²
    optimal_ncomps_q2 <- which.max(mean_q2)

    # Combine criteria: prioritize Q², but ensure RMSEP is within acceptable range
    optimal_ncomps <- ifelse(mean_rmsep[optimal_ncomps_q2] <= mean_rmsep[optimal_ncomps_rmsep], optimal_ncomps_q2, optimal_ncomps_rmsep)
  }

  # Plot RMSEP values
  if (isTRUE(plot)) {
    rmsep_df <- data.frame(
      Components = 1:max_comps,
      Mean_RMSEP = mean_rmsep,
      Lower_CI = mean_rmsep - se_rmsep,
      Upper_CI = mean_rmsep + se_rmsep
    )
    Comp_plot = ggplot(data = rmsep_df, aes(x = Components, y = Mean_RMSEP)) +
      geom_line(linewidth = 3) +
      geom_point(size = 6) +
      geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "blue", alpha = 0.3) +
      geom_vline(xintercept = optimal_ncomps, linetype = "dashed", color = "darkorange", linewidth = 2) +
      labs(title = "Cross-Validated Number of Components", x = "Number of Components", y = "Mean RMSEP") +
      theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5)) +
      ggpubr::theme_pubclean()
    print(Comp_plot)

    if(isTRUE(save.plot)) {
      optimal_ncomps = list(ncomps = optimal_ncomps, plot = Comp_plot)
    }
  }
  return(optimal_ncomps)
}

#' @title VIP: Estimate Variable Importance in Projections
#'
#' @description Use a VIP threshold to calculate variable importance in projections.
#'
#' @param plsr.results a PLS object from the following PLS packages: mixOmics, pls, plsdepot
#' @param ncomp the number of important latent components to consider
#' @param Scores a logical value indicating whether or not if the VIP scores are wanted (TRUE)
#' or if this function should output an already trimmed predictor dataset which
#' surpassed the VIPs threshold
#' @param vip_threshold a threshold calculation method for which VIPs are found to be significant.
#' The general rule of thumb is for the VIPs to be greater than one which is the
#' a priori assumption if they were to contribute to the latent space at all. The function
#' calculates the vips for each meaningful component, then takes the maximum value from
#' all of the latent spaces. Then, vip trimming is available through vip_threshold.
#' Options available are the mean vip ("mean"), the median vip ("median"), and the vips > 1 ("trim").
#' If vip_threshold = "NULL", no variable trimming will be performed, just the scores.
#'
#' @details This function is applicable to Partial Least Squares objects from the
#' following packages: mixOmics, pls, plsdepot.
#'
#' @keywords internal
#'
#' @export

VIP <- function(plsr.results, ncomp, Scores = FALSE, vip_threshold = "NULL") {
  Predictor = plsr.results$model[[2]]
  if("VIP" %in% names(plsr.results)) {
    W = plsr.results$raw.wgs
    H = ncomp
    q = ncol(plsr.results$model[[1]])
    p = ncol(plsr.results$model[[2]])
    vip = matrix(0, nrow = p, ncol = H)
    cor2 = cor(as.matrix.data.frame(plsr.results$model[[1]]), plsr.results$x.scores, use = "pairwise")^2
    cor2 = as.matrix(cor2, nrow = q)
    vip[, 1] = W[, 1]^2
    if (H > 1) {
      for (h in 2:H) {
        if (q == 1) {
          Rd = cor2[, 1:h]
          vip[, h] = Rd %*% t(W[, 1:h]^2)/sum(Rd)
        }
        else {
          Rd = apply(cor2[, 1:h], 2, sum)
          vip[, h] = Rd %*% t(W[, 1:h]^2)/sum(Rd)
        }
      }
    }
    vip = sqrt(p * vip)
    vip = as.matrix.data.frame(vip)
    rownames(vip) = colnames(plsr.results$model$X)
    colnames(vip) = paste0("t", 1:H)

    vip = data.frame(Predictor = rownames(vip), vip)

    if(vip_threshold == "trim") {
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score")) %>%
        filter(Score >= 1) %>%
        group_by(Predictor) %>%
        summarise(Score = mean(Score))
    } else if (vip_threshold == "median"){
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score")) %>%
        filter(Score >= median(Score)) %>%
        group_by(Predictor) %>%
        summarise(Score = mean(Score))
    } else if (vip_threshold == "mean"){
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score")) %>%
        filter(Score >= mean(Score)) %>%
        group_by(Predictor) %>%
        summarise(Score = mean(Score))
    } else if (vip_threshold == "NULL") {
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score"))
    }
    if(isTRUE(Scores)) {
      X_selected = vip
    } else {
      vip = unique(vip$Predictor)
      X_selected <- Predictor[, vip]
    }
  }
  if("loading.weights" %in% names(plsr.results)) {
    W = plsr.results$loading.weights
    H = plsr.results$ncomp
    q = length(plsr.results$Ymeans)
    p = length(plsr.results$Xmeans)
    vip = matrix(0, nrow = p, ncol = H)
    cor2 = cor(as.matrix.data.frame(plsr.results$model[[1]]), plsr.results$scores, use = "pairwise")^2
    cor2 = as.matrix(cor2, nrow = q)
    vip[, 1] = W[, 1]^2
    if (H > 1) {
      for (h in 2:H) {
        if (q == 1) {
          Rd = cor2[, 1:h]
          vip[, h] = Rd %*% t(W[, 1:h]^2)/sum(Rd)
        }
        else {
          Rd = apply(cor2[, 1:h], 2, sum)
          vip[, h] = Rd %*% t(W[, 1:h]^2)/sum(Rd)
        }
      }
    }
    vip = sqrt(p * vip)
    vip = as.matrix.data.frame(vip)
    rownames(vip) = colnames(plsr.results$model$X)
    colnames(vip) = paste0("t", 1:H)

    vip = data.frame(Predictor = rownames(vip), vip)

    if(vip_threshold == "trim") {
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score")) %>%
        filter(Score >= 1) %>%
        group_by(Predictor) %>%
        summarise(Score = mean(Score))
    } else if (vip_threshold == "median"){
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score")) %>%
        filter(Score >= median(Score)) %>%
        group_by(Predictor) %>%
        summarise(Score = mean(Score))
    } else if (vip_threshold == "mean"){
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score")) %>%
        filter(Score >= mean(Score)) %>%
        group_by(Predictor) %>%
        summarise(Score = mean(Score))
    } else if (vip_threshold == "NULL") {
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score"))
    }
    if(isTRUE(Scores)) {
      X_selected = vip
    } else {
      vip = unique(vip$Predictor)
      X_selected <- Predictor[, vip]
    }
  }
  if("loadings.star" %in% names(plsr.results)) {
    W = plsr.results$loadings$X
    H = plsr.results$ncomp
    q = ncol(plsr.results$Y)
    p = ncol(plsr.results$X)
    vip = matrix(0, nrow = p, ncol = H)
    cor2 = cor(plsr.results$Y, plsr.results$variates$X, use = "pairwise")^2
    cor2 = as.matrix(cor2, nrow = q)
    vip[, 1] = W[, 1]^2
    if (H > 1) {
      for (h in 2:H) {
        if (q == 1) {
          Rd = cor2[, 1:h]
          vip[, h] = Rd %*% t(W[, 1:h]^2)/sum(Rd)
        }
        else {
          Rd = apply(cor2[, 1:h], 2, sum)
          vip[, h] = Rd %*% t(W[, 1:h]^2)/sum(Rd)
        }
      }
    }
    vip = sqrt(p * vip)
    vip = as.matrix.data.frame(vip)
    rownames(vip) = colnames(plsr.results$model$X)
    colnames(vip) = paste0("t", 1:H)

    vip = data.frame(Predictor = rownames(vip), vip)

    if(vip_threshold == "trim") {
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score")) %>%
        filter(Score >= 1) %>%
        group_by(Predictor) %>%
        summarise(Score = mean(Score))
    } else if (vip_threshold == "median"){
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score")) %>%
        filter(Score >= median(Score)) %>%
        group_by(Predictor) %>%
        summarise(Score = mean(Score))
    } else if (vip_threshold == "mean"){
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score")) %>%
        filter(Score >= mean(Score)) %>%
        group_by(Predictor) %>%
        summarise(Score = mean(Score))
    } else if (vip_threshold == "NULL") {
      vip = data.frame(pivot_longer(vip, cols = starts_with("t"),
                                    names_to = "Response", values_to = "Score"))
    }
    if(isTRUE(Scores)) {
      X_selected = vip
    } else {
      vip = unique(vip$Predictor)
      X_selected <- Predictor[, vip]
    }
  }
  return(X_selected)
} # end of function


#' @title Create a k-fold segmentation
#'
#' @description Creates a list of row indices for a dataset by which a k-fold
#' segmentation is performed.
#'
#' @param data a data table or data frame for k-fold segmentation
#' @param k a numeric indice indicating how many folds should be performed.
#'
#' @keywords internal
#'
#' @export


create_segments <- function(data, k = 5) {

  shuffled_indices <- sample(seq_len(nrow(data)))

  fold_size <- floor(nrow(data) / k)
  segments <- vector("list", k)

  for (i in seq_len(k)) {
    if (i != k) {
      segments[[i]] <- shuffled_indices[((i - 1) * fold_size + 1):(i * fold_size)]
    } else {
      segments[[i]] <- shuffled_indices[((i - 1) * fold_size + 1):nrow(data)]
    }
  }
  return(segments)
}

#' Create and Save a Correlation Network Plot
#'
#' This function generates a network visualization of correlations between two groups
#' of variables (rows and columns of the input matrix). It automatically filters
#' weak correlations, adjusts visual parameters based on network density, and saves
#' a high-resolution PNG image.
#'
#' @param matrix A numeric matrix where rows represent variables in Group 1 and
#'   columns represent variables in Group 2.
#' @param g1name A character string. The label for the first group (rows). Defaults to "Group 1".
#' @param g2name A character string. The label for the second group (columns). Defaults to "Group 2".
#' @param title A character string. The main title of the plot. Defaults to
#'   'Correlation Network Plot Comparison'.
#' @param threshold A numeric value. The absolute correlation coefficient cutoff.
#'   Edges with weights below this value are excluded. If the threshold is higher
#'   than the maximum observed correlation, it is automatically lowered to the maximum.
#'   Defaults to 0.7.
#' @param colors A character vector of length 2. Hex codes or names for the node
#'   colors (Group 1, Group 2). Defaults to c("#FF8C00", "#1E90FF").
#' @param res Numeric. The resolution (DPI) for the saved PNG. Defaults to 600.
#' @param width Numeric. The width of the saved image in inches. Defaults to 10.
#' @param height Numeric. The height of the saved image in inches. Defaults to 8.
#'
#' @return A \code{ggplot} object representing the network graph.
#'
#' @section Side Effects:
#' Saves a PNG file to the working directory with the naming convention:
#' \code{Cor_Network_[Group1]_vs_[Group2].png}.
#'
#' @importFrom igraph graph_from_data_frame delete_vertices degree V vcount
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text scale_edge_width theme_graph
#' @importFrom ggplot2 aes scale_color_manual labs theme element_text ggsave
#' @importFrom ragg agg_png
#' @importFrom rlang %||%
#'
#' @export

cn <- function(matrix,
               g1name = "Group 1",
               g2name = "Group 2",
               title = 'Correlation Network Plot Comparison',
               threshold = 0.6,
               colors = c("orange", "royalblue"),
               res = 600,
               width = 10,
               height = 8) {


  matrix[is.na(matrix)] <- 0
  max_corr <- max(abs(matrix), na.rm = TRUE)

  # If the user's threshold is too high, lower it automatically
  if (max_corr < threshold) {
    message(sprintf("Notice: No correlations met the threshold of %.2f.", threshold))
    message(sprintf("Auto-lowering threshold to %.2f (maximum observed correlation).", max_corr))
    threshold <- max_corr
  }

  n_vars_g1 <- nrow(matrix)
  n_vars_g2 <- ncol(matrix)
  total_vars <- n_vars_g1 + n_vars_g2

  # Build edge list including direction
  edges_df <- expand.grid(from_idx = 1:n_vars_g1, to_idx = 1:n_vars_g2)
  edges_df$weight <- as.vector(matrix)
  edges_df <- edges_df[!is.na(edges_df$weight) & abs(edges_df$weight) >= threshold, ]

  if(nrow(edges_df) == 0) {
    message("No correlations met the threshold.")
    return(NULL)
  }

  # Create Direction column for the legend
  edges_df$dir <- ifelse(edges_df$weight > 0, "Positive", "Negative")
  edges_df$dir <- factor(edges_df$dir, levels = c("Positive", "Negative"))
  edges_df$to_idx <- edges_df$to_idx + n_vars_g1
  edges_df$weight <- abs(edges_df$weight)

  # Node data
  nodes_df <- data.frame(
    id = 1:total_vars,
    name = c(rownames(matrix) %||% paste0("Row", 1:n_vars_g1),
             colnames(matrix) %||% paste0("Col", 1:n_vars_g2)),
    group = c(rep(g1name, n_vars_g1), rep(g2name, n_vars_g2))
  )

  g <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes_df)
  g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])

  n_final <- igraph::vcount(g)

  # Auto-adjust visual parameters
  if (n_final <= 10) {
    p_size <- 12; t_size <- 8
  } else if (n_final >= 11) {
    p_size <- 11; t_size <- 7
  } else if (n_final >= 30) {
    p_size <- 10; t_size <- 6
  } else if (n_final >= 40) {
    p_size <- 9; t_size <- 5
  } else if (n_final >= 50) {
    p_size <- 8; t_size <- 4
  } else if (n_final >= 60) {
    p_size <- 7; t_size <- 3
  } else {
    p_size <- 5; t_size <- 1
  }

  # Plot with ggraph
  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight, linetype = dir),
      color = "grey",
      alpha = 0.5
    ) +
    ggraph::geom_node_point(ggplot2::aes(color = group), size = p_size) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      repel = TRUE,
      #size = t_size,
      fontface = "bold",
      family = "sans"
    ) +
    ggraph::scale_edge_linetype_manual(
      values = c("Positive" = "solid", "Negative" = "dashed"),
      drop = FALSE
    ) +
    # Aesthetic Controls
    ggplot2::scale_color_manual(values = colors) +
    ggraph::scale_edge_width(range = c(0.2, 2.0)) +
    # Legend Refinement
    ggplot2::labs(
      title = title,
      color = expression(bold("Comparisons")),
      edge_width = "Strength",
      edge_linetype = "Direction",
    ) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
    )

    # Save with ragg for multiple device compatibility
    clean_name1 <- gsub("[^A-Za-z0-9]", "", g1name)
    clean_name2 <- gsub("[^A-Za-z0-9]", "", g2name)
    filename <- paste0("CN_", clean_name1, "_vs_", clean_name2, ".png")
    ggplot2::ggsave(filename, plot = p, width = width, height = height, dpi = res, device = ragg::agg_png)

  return(p)
}

#' @title Sparse Partial Least Squares from MixOmics Package
#'
#' @description
#' A standalone implementation of Sparse Partial Least Squares (sPLS) regression
#' designed to replace the mixOmics dependency. It performs simultaneous variable
#' selection and dimensionality reduction. This function structures the output
#' to be compatible with the custom VIP and statistics functions in this package.
#'
#' @references
#' Le Cao K-A, Rossouw D, Robert-Granie C, Besse P (2008). A sparse PLS for
#' variable selection when integrating Omics data. Statistical Applications in
#' Genetics and Molecular Biology 7(1).
#'
#' @param X Numeric matrix or data frame of predictors.
#' @param Y Numeric matrix or data frame of responses.
#' @param ncomp Integer, number of latent components (default 2).
#' @param mode Character, PLS mode. Default is "regression".
#' @param keepX Numeric vector of length ncomp, number of variables to keep in X loadings.
#' @param keepY Numeric vector of length ncomp, number of variables to keep in Y loadings.
#' @param scale Logical, whether to scale the data (default TRUE).
#' @param near.zero.var Logical, whether to remove zero variance columns (default TRUE).
#' @param max.iter Integer, maximum iterations for convergence (default 500).
#' @param tol Numeric, tolerance for convergence (default 1e-06).
#'
#' @return An object of class "mixo_spls" (and "list") containing loadings,
#' variates (scores), and model data compatible with downstream VIP analysis.
#'
#' @export
spls <- function(X, Y,
                 ncomp = 2,
                 mode = "regression",
                 keepX = NULL,
                 keepY = NULL,
                 scale = TRUE,
                 near.zero.var = TRUE,
                 max.iter = 500,
                 tol = 1e-06) {

  # --- Helper: Zero Variance Removal ---
  clean_data <- function(mat) {
    if (is.null(mat)) return(NULL)
    mat <- as.matrix(mat)
    if (near.zero.var) {
      vars <- apply(mat, 2, var, na.rm = TRUE)
      keep_cols <- vars > .Machine$double.eps
      if (sum(keep_cols) < 2) warning("Variables with near-zero variance detected and removed, leaving < 2 variables.")
      return(mat[, keep_cols, drop = FALSE])
    }
    return(mat)
  }

  # --- Helper: Soft Thresholding for Sparsity ---
  select_var <- function(vec, keep) {
    if (is.null(keep) || keep >= length(vec)) return(vec)
    absv <- abs(vec)
    threshold <- sort(absv, decreasing = TRUE)[keep]
    vec[absv < threshold] <- 0
    return(vec)
  }

  # --- 1. Data Preparation ---
  X <- clean_data(X)
  Y <- clean_data(Y)

  X.names <- colnames(X)
  Y.names <- colnames(Y)

  # Store original matrices for model output
  X_orig <- X
  Y_orig <- Y

  # Center and Scale (and save means/sd for prediction/VIP)
  X_mean <- colMeans(X, na.rm = TRUE)
  X_sd <- apply(X, 2, sd, na.rm = TRUE)
  if(!scale) X_sd <- rep(1, ncol(X))

  Y_mean <- colMeans(Y, na.rm = TRUE)
  Y_sd <- apply(Y, 2, sd, na.rm = TRUE)
  if(!scale) Y_sd <- rep(1, ncol(Y))

  X <- scale(X, center = X_mean, scale = X_sd)
  Y <- scale(Y, center = Y_mean, scale = Y_sd)

  # Prepare output containers
  n_X <- ncol(X)
  n_Y <- ncol(Y)

  loadings_X <- matrix(0, nrow = n_X, ncol = ncomp)
  loadings_Y <- matrix(0, nrow = n_Y, ncol = ncomp)
  variates_X <- matrix(0, nrow = nrow(X), ncol = ncomp)
  variates_Y <- matrix(0, nrow = nrow(Y), ncol = ncomp)
  mat_c      <- matrix(0, nrow = n_Y, ncol = ncomp) # Y loadings (for prediction)
  mat_d      <- matrix(0, nrow = n_X, ncol = ncomp) # X loadings (for deflation)

  # Row/Col names
  rownames(loadings_X) <- X.names
  rownames(loadings_Y) <- Y.names
  colnames(loadings_X) <- colnames(loadings_Y) <- paste0("comp", 1:ncomp)

  # Default keepX/keepY
  if (is.null(keepX)) keepX <- rep(n_X, ncomp)
  if (is.null(keepY)) keepY <- rep(n_Y, ncomp)
  if (length(keepX) < ncomp) keepX <- rep(keepX, length.out = ncomp)
  if (length(keepY) < ncomp) keepY <- rep(keepY, length.out = ncomp)

  # Create deflated matrices for iteration
  X_curr <- X
  Y_curr <- Y

  # --- 2. Iterative NIPALS / SVD Loop ---
  for (h in 1:ncomp) {

    # 2a. SVD on Cross-Product Matrix
    M <- crossprod(X_curr, Y_curr)
    svd_M <- svd(M, nu = 1, nv = 1)
    u_h <- svd_M$u[, 1]
    v_h <- svd_M$v[, 1]

    # 2b. Apply Sparsity
    u_h <- select_var(u_h, keepX[h])
    v_h <- select_var(v_h, keepY[h])

    # Normalize
    if(sum(u_h^2) > 0) u_h <- u_h / sqrt(sum(u_h^2))
    if(sum(v_h^2) > 0) v_h <- v_h / sqrt(sum(v_h^2))

    # 2c. Calculate Scores
    xi_h <- X_curr %*% u_h
    omega_h <- Y_curr %*% v_h

    # 2d. Calculate Deflation Vectors (Loadings C and D)
    denom <- sum(xi_h^2)
    c_h <- crossprod(Y_curr, xi_h) / denom # Y loading
    d_h <- crossprod(X_curr, xi_h) / denom # X loading

    # 2e. Deflation
    X_curr <- X_curr - xi_h %*% t(d_h)
    Y_curr <- Y_curr - xi_h %*% t(c_h)

    # 2f. Store results
    loadings_X[, h] <- u_h
    loadings_Y[, h] <- v_h
    variates_X[, h] <- xi_h
    variates_Y[, h] <- omega_h
    mat_c[, h]      <- c_h
    mat_d[, h]      <- d_h
  }

  # --- 3. Format Output ---
  res <- list(
    call = match.call(),
    X = X_orig,
    Y = Y_orig,
    ncomp = ncomp,
    mode = mode,
    keepX = keepX,
    keepY = keepY,

    # MixOmics standard slots
    loadings = list(X = loadings_X, Y = loadings_Y),
    variates = list(X = data.frame(variates_X), Y = data.frame(variates_Y)),
    mat.c = mat_c, # Required for prediction
    mat.d = mat_d, # Required for prediction/deflation

    # Aliases for VIP and other utils.R functions
    loading.weights = loadings_X,
    scores = data.frame(variates_X),

    # Scaling info (Required for VIP and Predict)
    Xmeans = X_mean,
    Ymeans = Y_mean,
    Xsigma = X_sd, # Used in prediction scaling
    Ysigma = Y_sd,

    names = list(sample = rownames(X), colnames = list(X = X.names, Y = Y.names)),
    tol = tol,
    max.iter = max.iter,
    iter = rep(0, ncomp)
  )

  # Model slot required by pls.bio.R (Ordered Y then X)
  res$model <- list(Y = Y_orig, X = X_orig)

  # Explained Variance
  expl_X <- apply(variates_X, 2, var) / sum(apply(X, 2, var))
  expl_Y <- apply(variates_Y, 2, var) / sum(apply(Y, 2, var))
  res$prop_expl_var <- list(X = expl_X, Y = expl_Y)

  class(res) <- "mixo_spls"

  return(res)
}

#' @title Predict Method for sPLS through MixOmics Package
#'
#' @description
#' S3 method to enable prediction for mixo_spls objects created by the standalone spls function.
#' This is required for plsr_stats to function correctly.
#'
#' @param object A mixo_spls object.
#' @param newdata Matrix or data frame of new predictor data.
#' @param ... Additional arguments.
#'
#' @export
predict.mixo_spls <- function(object, newdata, ...) {

  newdata <- as.matrix(newdata)

  # Scale new data using training statistics
  # Handling cases where scale=FALSE (sigma=1)
  newdata <- scale(newdata, center = object$Xmeans, scale = object$Xsigma)

  # Containers
  ncomp <- object$ncomp
  n_samples <- nrow(newdata)
  n_Y <- ncol(object$Y)

  # Array to store predictions: Samples x Variables x Components
  pred_array <- array(0, dim = c(n_samples, n_Y, ncomp),
                      dimnames = list(rownames(newdata), colnames(object$Y), paste0("dim", 1:ncomp)))

  # Iterative prediction based on deflation
  X_curr <- newdata
  Y_pred_cum <- matrix(0, nrow = n_samples, ncol = n_Y)

  for (h in 1:ncomp) {
    # Calculate score for new data
    t_h <- X_curr %*% object$loadings$X[, h]

    # Predict Y part for this component (t * c')
    # c is stored in mat.c
    y_pred_h <- t_h %*% t(object$mat.c[, h])

    # Accumulate predictions
    Y_pred_cum <- Y_pred_cum + y_pred_h

    # Deflate X using d (stored in mat.d)
    X_curr <- X_curr - t_h %*% t(object$mat.d[, h])

    # Unscale the prediction to original units
    # Y_raw = Y_scaled * sigma + mean
    Y_pred_raw <- sweep(Y_pred_cum, 2, object$Ysigma, "*")
    Y_pred_raw <- sweep(Y_pred_raw, 2, object$Ymeans, "+")

    pred_array[, , h] <- Y_pred_raw
  }

  # Return list with $predict element to match plsr_stats expectation
  list(predict = pred_array)
}

#' @title tune.spls: Estimate the maximum number of variables in X to keep
#'
#' @description
#' Minimal implementation of tune.spls to satisfy dependency calls.
#' Returns the largest keepX to maximize information retention.
#'
#' @export
tune.spls <- function(X, Y, ncomp, test.keepX, ...) {
  choice.keepX <- rep(max(test.keepX), ncomp)
  names(choice.keepX) <- paste0("comp", 1:ncomp)
  list(choice.keepX = choice.keepX)
}


#' @title condense multiple data sheets into one
#'
#' @description This function takes a nested list object containing data table
#' results from multiple subsetted tests with matching column names and condenses
#' it down into one single excel file with each sheet as a separate regression model.
#' This is an internal function specifically used to reconstruct the summary
#' outputs for the AABV and PSLR functions.The function is ideal only when
#' subsetted tests of variables are performed.
#'
#' @param AABV.summary Either an AABV or PSLR function output.
#'
#' @param write.results Logical value whether to summarize the results to a
#' readable Excel files. Each model containing the associated regression results will
#' be a separate Excel file with each regression result as a separate Excel sheet.
#' Default is TRUE.
#'
#' @param key Character string. Identifier for output files and directories.
#' Default is NULL.
#'
#' @details This function takes either an AABV or PSLR object that has the argument
#' restructure.models set to FALSE (this argument automatically performs this function).
#' The function is ideal when subsetted tests of variables are performed. The function
#' will iterate over each regression model/hypothesis and each unique regression test.
#' If the regression test has subsetted variables that were each separately subject
#' to regression tests, the function will combine each subsetted test into a single
#' excel sheet within an excel file for that regression model/hypothesis.
#'
#' @returns a nested list object containing each overall regression hypothesis/model
#' as a nested list of individual, subsetted regression tests.
#'
#' @author Keeling et al., 2025
#'
#' @importFrom writexl write_xlsx
#' @importFrom dplyr filter mutate
#' @importFrom utils write.csv
#'
#' @export
#'
#' @import geomorph writexl readxl
#'
#' @examples
#'
#' #' \dontrun{
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
#' \itemize{
#'   \item{Dimensions}{Variables 1:3, 5:8, 10:15}
#'   \item{Orientation}{Variables 9, 24}
#'   \item{Cortical Thickness}{Variables 4, 16:18, 30:60}
#'   \item{Bending Resistance}{Variables 19:23}
#'   \item{Breaking Strength}{Variables 25:29}
#' }
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
#' # Run the PSLR function to conduct parallel regressions.
#'
#' output <- PSLR(Models, # the regression model created above
#'                subset = TRUE, # there are subsets in this regression test
#'                paired_subsets = FALSE, # Since response variable contains no unique subsets, this is not necessary
#'                point_set, # this is the subset model that was created above
#'                key = "Corpus", # Name for the output file folder.
#'                restructured.results = FALSE # The function restructure_models performs this
#'                ) # keep everything else as default
#'
#'
#' }
#'
#' # Restructure the results
#'
#' output <- restructure_models(output)
#'
#' # View the results
#'
#' view(output[[1]][[1]]) # view the first set of regression tests from the first model

restructure_models <- function(AABV.summary,
                               write.results = TRUE,
                               key = NULL) {
  Dir = getwd()
  if(!is.null(key)) {
    path = paste("AABV_Restructured-", key)
    path = file.path(Dir, path)
    if(!dir.exists(path)) {
      dir.create(path)
    }

  } else {
    path = file.path(Dir, "AABV_Restructured")
    if(!dir.exists("AABV_Restructured")) {
      dir.create("AABV_Restructured")
    }
  }
  setwd(path)

  # Initialize a new list to store the restructured data
  restructured_output <- list()

  # Loop through each model in the output list
  for (model_name in names(AABV.summary)) {
    # Estimate model structure
    num_comparisons <- length(AABV.summary[[model_name]])
    num_subsets <- nrow(AABV.summary[[model_name]][[1]])
    # Initialize the model's list in the restructured output
    restructured_output[[model_name]] <- vector("list", num_subsets)

    # Initialize empty dataframes for each subset
    for (subset in 1:num_subsets) {
      restructured_output[[model_name]][[subset]] <- data.frame(matrix(nrow = num_comparisons, ncol = ncol(AABV.summary[[model_name]][[1]])))
      colnames(restructured_output[[model_name]][[subset]]) <- colnames(AABV.summary[[model_name]][[1]])
    }
    # Loop through each comparison
    for (comparison in 1:num_comparisons) {
      # Get the dataframe for the current comparison
      comparison_df <- AABV.summary[[model_name]][[comparison]]
      if (!is.data.frame(comparison_df)) {
        stop("The data structure is not as expected. Each comparison should be a dataframe.")
      }
      # Loop through each subset
      for (subset in 1:num_subsets) {
        # Place the corresponding row of the comparison dataframe into the correct position
        restructured_output[[model_name]][[subset]][comparison, ] <- comparison_df[subset, ]
      }
    }
    # Name the list elements for each model
    names(restructured_output[[model_name]]) <- paste0("Subset_", 1:num_subsets)
  }
  if(isTRUE(write.results)) {
    for(m in seq_along(restructured_output)) {
      write_xlsx(restructured_output[[m]], path = paste("Restructured AABV_Summary of Model_",m,".xlsx", sep = ""))
    }
  }
  return(restructured_output)
}

