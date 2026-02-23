#' @title reg.bio - A Robust Multiple, Multivariate Regression Algorithm
#'
#' @description
#' The reg.bio function is a comprehensive regression tool designed for conducting complex
#' multivariate regression analyses in cross-sectioanl geometric morphometric
#' applications. It incorporates traditional regression approaches and modern machine learning methods
#' by implementing both classification and prediction capabilities within a unified
#' framework. The function particularly excels at handling high-dimensional data common
#' in morphometric studies, offering robust solutions for both linear and non-linear
#' relationships between variables.
#'
#' At its core, the function utilizes support vector machines (SVM) through either the
#' caret or e1071 packages, enhanced with permutation testing capabilities to provide
#' statistical rigor in morphometric analyses. It automatically handles data
#' preprocessing, model validation, and statistical testing, making it particularly
#' suitable for analyzing shape variation and morphological data. The function can
#' process both single and multiple response variables, making it versatile for various
#' research scenarios in geometric morphometrics.
#'
#' The function conducts model validation techniques through
#' cross-validation and bootstrapping, along with permutation testing for statistical
#' significance. These features ensure robust results even with the complex,
#' multivariate data structures common in morphometric analyses. The function also
#' provides comprehensive model diagnostics and performance metrics, crucial for
#' understanding the reliability and accuracy of morphometric predictions.
#'
#' @param Response A matrix or data frame containing the response variables to be predicted
#' @param Predictor A matrix or data frame containing the predictor variables
#' @param type Character; type of analysis to perform: "pred" for prediction or "class" for classification
#' @param package Character; specifies the package to use: "svm" (e1071) or "caret"
#' @param p_calc Character; method for p-value calculation ("RMSE", "R2", "F", "gcv", or "coef")
#' @param n_boot Integer; number of bootstrap iterations (default: 0)
#' @param perm Integer; number of permutations for significance testing (default: 500)
#' @param CV Logical; whether to perform cross-validation (default: TRUE)
#' @param which_classify Vector; optional indices for specific observations to classify
#' @param full_model Logical value indicating whether to use full model without variable trimming. Default is FALSE to perform test using trimmed model rather than using colinear and/or uninformative variables.
#' @param VIPS Vector; optional variable importance scores
#' @param parallel Logical value whether to use parallel processing (default: TRUE)
#' @param core_choice Character value specifying number of cores to use.
#' The following options are available:
#' \itemize{
#'    \item{"detect"} {:  total available cores - 2}
#'    \item{"high"} {:  28 cores}
#'    \item{"medium-high"} {:   20 cores}
#'    \item{"medium"} {:  14 cores}
#'    \item{"medium-low"} {:  10 cores}
#'    \item{"low"} {:   8 cores}
#' } ("detect" , "high", "medium-high", "medium", "medium-low", "low", or specific number)
#' @param print_progress Logical value whether to display progress bars (default: TRUE)
#' @param ... Additional arguments passed for \code{DT}, \pkg{caret}, and \pkg{e1071} parameters
#'
#' @details
#' The reg.bio function operates through a sophisticated multi-step process that ensures
#' robust analysis of morphometric data. Initially, it performs extensive data validation
#' and preprocessing, including automatic handling of factor levels and data type
#' conversions. The function supports both prediction and classification tasks,
#' automatically adjusting its internal algorithms based on the nature of the response
#' variables.
#'
#' For model fitting, the function employs a flexible approach that can utilize either
#' the caret or e1071 package's implementation of support vector machines. When using
#' caret, it takes advantage of the package's extensive model tuning capabilities,
#' automatically optimizing model parameters through cross-validation. With e1071, it
#' implements direct SVM fitting with customizable kernel functions and parameter settings.
#'
#' The validation process integrates several critical components. When cross-validation
#' is enabled, the function implements a k-fold validation scheme, systematically
#' evaluating model performance across different data partitions. The bootstrapping
#' option provides an additional layer of validation, generating confidence intervals
#' for model parameters and predictions. These validation procedures are particularly
#' important in morphometric analyses where understanding the stability and reliability
#' of shape predictions is crucial.
#'
#' Statistical significance is assessed through a sophisticated permutation testing
#' framework. The function generates null distributions by randomly permuting the
#' response variables while maintaining the correlation structure of the predictors.
#' This approach provides robust p-values that account for the multivariate nature of
#' morphometric data. The permutation process includes an adaptive convergence check,
#' potentially reducing computational time while maintaining statistical rigor.
#'
#' The function implements parallel processing capabilities to handle computationally
#' intensive tasks efficiently. It automatically manages parallel processes based on
#' the specified core choice, with built-in options for different levels of parallel
#' processing intensity. This feature is particularly valuable when working with large
#' morphometric datasets or when performing extensive permutation testing.
#'
#' Variable importance analysis is integrated into the function, allowing researchers
#' to identify which morphometric variables contribute most significantly to the
#' observed patterns. This analysis can be performed through various methods, including
#' coefficient analysis and dedicated variable importance calculations, depending on
#' the chosen modeling approach.
#'
#' The function also handles PCA transformations for response
#' variables, which is particularly useful in geometric morphometric analyses where
#' dimensional reduction is often necessary. It automatically manages the transformation
#' and back-transformation of data, ensuring that results are presented in the original
#' coordinate space while taking advantage of the statistical benefits of PCA.
#'
#' Error handling and input validation are comprehensive throughout the function, with
#' informative messages and warnings to guide users in correctly specifying analyses.
#' The function also implements memory management strategies to handle large datasets
#' efficiently, including cleanup of intermediate results and optimal storage of model
#' components.
#'
#' @return A list of class "pred" containing:
#' \itemize{
#'   \item Model:         Formula representation of the fitted model
#'   \item Coefficients:  Matrix of model coefficients
#'   \item Prediction:    Matrix or array of predicted values
#'   \item summary (when type == "pred") - Data frame containing model statistics including:
#'     \itemize{
#'       \item n_perm:      Number of permutations used
#'       \item n_boot:      Number of bootstrap iterations
#'       \item df:          Degrees of freedom
#'       \item F_stat:      F-statistic
#'       \item RMSE:        Root mean square error
#'       \item R2:          R-squared value
#'       \item BIC_test:    Bayesian Information Criterion
#'       \item p_value:     Statistical significance of the model.
#'     }
#'   \item summary (when type == "class") - Data frame containing model statistics including:
#'     \itemize{
#'       \item n_perm:                Number of permutations used.
#'       \item n_boot:                Number of bootstrap iterations.
#'       \item df:                    Degrees of freedom.
#'       \item Sensitivity:           Proportion of actual positives correctly identified (True Positive Rate).
#'       \item Specificity:           Proportion of actual negatives correctly identified (True Negative Rate).
#'       \item Pos Pred Value:        Proportion of predicted positives that are actual positives.
#'       \item Neg Pred Value:        Proportion of predicted negatives that are actual negatives.
#'       \item Precision:             How precise the positive predictions are.
#'       \item Recall:                Ability of the model to identify positive cases.
#'       \item F1:                    Harmonic mean of Precision and Recall
#'       \item Prevalence:            Proportion of positive cases in the dataset.
#'       \item Detection Rate:        Proportion of true positives out of the total observations.
#'       \item Detection Prevalence:  Proportion of predicted positive cases out of the total observations.
#'       \item Balanced Accuracy:     Mean of Sensitivity and Specificity, useful for imbalanced datasets.
#'       \item p_value: Indicates     Statistical significance of the model.
#'     }
#' }
#'
#' @importFrom caret trainControl train predict.train createDataPartition createFolds varImp
#' @importFrom e1071 svm tune.svm
#' @importFrom stats mahalanobis cov sd median p.adjust pnorm predict coef quantile na.omit cor as.formula
#' @importFrom glmnet cv.glmnet coef.glmnet
#' @importFrom randomForest randomForest
#' @importFrom dplyr filter mutate arrange %>% across where intersect
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom progress progress_bar
#' @importFrom kernlab ksvm
#' @importFrom kerndwd kerndwd
#' @importFrom tidyr pivot_longer
#' @importFrom geomorph two.d.array gm.prcomp arrayspecs
#' @importFrom purrr is_empty
#' @importFrom utils head tail
#' @importFrom methods is
#' @importFrom earth earth
#'
#' @seealso
#' \code{\link[caret]{train}} for the underlying caret implementation
#' \code{\link[e1071]{svm}} for the underlying SVM implementation
#'
#' @author Keeling et al., 2025
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage for prediction
#' result <- pred(Response = y_data,
#'               Predictor = x_data,
#'               type = "pred")
#'
#' # Classification with cross-validation
#' class_result <- pred(Response = class_data,
#'                     Predictor = predictors,
#'                     type = "class",
#'                     CV = TRUE,
#'                     n_boot = 100)
#'
#' # Parallel processing with custom cores
#' parallel_result <- pred(Response = y_data,
#'                        Predictor = x_data,
#'                        parallel = TRUE,
#'                        core_choice = "medium-high")
#' }
#'

reg.bio <- function(Response,
                 Predictor,
                 type = "pred",
                 p_calc = "gcv",
                 n_boot = 25,
                 perm = 500,
                 CV = TRUE,
                 which_classify = NULL,
                 full_model = TRUE,
                 VIPS = NULL,
                 parallel = TRUE,
                 core_choice = "detect",
                 print_progress = TRUE,
                 ...) {

  closeAllConnections() # stop all background processes
  args = list(...) # create list for extra arguments
  output = list() # create function return object
  #__________________________________________________________________________________

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

  # Print Progress__________________________________________________________________________________
  if (isTRUE(print_progress)) {
    pb_start <- progress::progress_bar$new(
      format = "Processing :what [:bar] :percent ETA: :eta",
      clear = FALSE, total = 100, width = 80)
  }

  if (isTRUE(print_progress)) {
    pb_start$tick(10, tokens = list(what = "Regression Data Transformation:"))
  }
  # Parallel processing cores________________________________________________________________

  if(core_choice == "detect") {
    core_choice = parallel::detectCores()-3
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

  clust <- parallel::makeCluster(core_choice)

  # Data Transformation prior to regression______________________________________________________________

  if (isTRUE(print_progress)) {
    pb_start$tick(10, tokens = list(what = "Data Transformation:"))
  }

  dt_params = names(formals(DT))

  dt_parameters <- list(Response = Response, Predictor = Predictor, dt_method = "procdist",
                        Res_transform = "NULL", Pred_transform = "NULL",
                        Res_ncomp = NULL, Pred_ncomp = NULL,
                        Res_point_set = NULL, Pred_point_set = NULL)

  dt_parameters <- modifyList(dt_parameters, args) %>% .[intersect(names(.), dt_params)]

  Output <- do.call(DT, dt_parameters)

  if(length(dim(Response)) == 3) {
    res_3d = TRUE
  } else {
    res_3d = FALSE
  }

  Predictor = data.frame(Output$Predictor)
  Predictor <- Predictor %>% mutate(across(where(is.logical), as.integer))
  Response = data.frame(Output$Response)

  if(type == "class") {
    Response = data.frame(Response)
    Response <- as.data.frame(lapply(Response, function(x) {
      factor(x, levels = unique(x))
    }))
  } else {
    Response = data.frame(Response)
    Response <- Response %>% mutate(across(where(is.logical), as.integer))
  }
  Response <- data.frame(Response)

  if(startsWith("pca", dt_parameters$Res_transform) || startsWith("bgpca", dt_parameters$Res_transform)) {
    Res_dim = Output$Res_dim
    Res_p = Output$Res_p
    Res_pca = Output$Res_pca
    Res_rot = as.matrix.data.frame(Output$Res_rot)
    Res_pca_scores = Output$Res_pca_scores
    Res_pca_center = Output$Res_pca_center
  }

  if(grepl("pca", dt_parameters$Res_transform) & grepl("pca",dt_parameters$Pred_transform)) {
    rank = TRUE
  } else {
    rank = FALSE
  }
  # Caret defaults: modifying parameters if needed___________________________________________________________

  if (isTRUE(print_progress)) {
    pb_start$tick(10, tokens = list(what = "Setting Regression Parameters:"))
  }

  if(type == "pred") {

  caret_params = c("x","y","method","preProcess","weights","metric","maximize",
                   "trControl","tuneGrid","tuneLength","mtry","replace","sampsize",
                   "nodesize","maxnodes","importance","localImp","nPerm","proximity",
                   "oob.prox","norm.votes","keep.forest","keep.inbag")

  caret_parameters = list(x = Predictor, y = apply(Response, 2, as.factor),
                          method = "gcvEarth", tuneLength = 10,
                          metric = "RMSE", tuneGrid = NULL, maximize = FALSE,
                          trControl = caret::trainControl(method = "repeatedcv", repeats = 5, number = 5))

  caret_parameters = modifyList(caret_parameters, args) %>% .[intersect(names(.), caret_params)]

  family = "gaussian"

  } else {

    caret_params = c("x","y","method","preProcess","weights","metric","maximize",
                     "trControl","tuneGrid","tuneLength","mtry","replace","sampsize",
                     "nodesize","maxnodes","importance","localImp","nPerm","proximity",
                     "oob.prox","norm.votes","keep.forest","keep.inbag")

    caret_parameters = list(x = Predictor, y = apply(Response,2,as.factor),
                            method = 'rf', tuneLength = 10,
                            metric = "Accuracy",
                            tuneGrid = NULL, maximize = TRUE,
                            preProcess = NULL,
                            trControl = caret::trainControl(method = "none"))

    caret_parameters = modifyList(caret_parameters, args) %>% .[intersect(names(.), caret_params)]

    if(is.logical(Response)) {
      family = "binomial"
    } else {
      family = "multinomial"
    }
  }

  if(!is.null(which_classify)) { # in case some observations shouldn't be considered in regression
    res_class = data.frame(Response[which_classify])
    pred_class = Predictor[which_classify,]
    Response = Response[-c(which_classify)]
    Predictor = Predictor[-c(which_classify),]
    pred_special = matrix(NA, ncol = ncol(res_class), nrow = nrow(res_class))
  }

  # Start Regression____________________________________________________________________________________________________________________

  if (isTRUE(print_progress)) {
    pb_start$tick(10, tokens = list(what = "Tuning Regression Parameters:"))
  }

  # Full Model Regression______________________________________________________________________________________________________________

  if(!"trControl" %in% names(caret_parameters)){
    trControl = suppressWarnings(caret::trainControl(method = "repeatedcv", number = 5, repeats = 5))
    caret_parameters$trControl <- trControl
  }

  if(isTRUE(full_model)) {

    pred_train_o = matrix(NA, ncol = ncol(Response), nrow = nrow(Response))
    coefs = matrix(0, nrow = ncol(Response), ncol = ncol(Predictor))
    colnames(coefs) <- colnames(Predictor)
    place_coef = matrix(NA, nrow = 1, ncol = ncol(Predictor))
    gcv = matrix(NA, nrow = 1, ncol = ncol(Response))
    if(p_calc == "gcv") {
      #p_calc = "coef"
    }
    lasso.coefs = list()
    coef_id = FALSE

    models = list()
    reg.parameter = list()

    doParallel::registerDoParallel(cores = core_choice)

    for (j in 1:ncol(Response)) { # loop over each response variable individually

      caret_parameters = modifyList(caret_parameters, list(y = Response[,j],
                                                           x = as.matrix.data.frame(Predictor)))

      pred_perm = suppressWarnings(do.call(caret::train, caret_parameters))

      pred_train_o[,j] = predict.train(pred_perm, newdata = Predictor)

      if(!is.null(which_classify)) {
        pred_special[,j] = predict.train(pred_perm, newdata = pred_class)
      }

      #models[[j]] <- pred_perm # takes a lot of storage

      if(p_calc == "gcv") {

        gcv[1,j] <- if(startsWith(pred_perm$method, "svm")) {
          tryCatch({pred_perm$finalModel@error},
                   error = function(e) {
                     n = nrow(Response[,j])
                     edf = pred_perm$finalModel@nSV
                     rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                     rss / (n * (1 - edf/n)^2)
                   })
        } else if("fit" %in% names(pred_perm$finalModel)) {
          mean(unlist(lapply(pred_perm$finalModel$fit, function(x) x$gcv)))
        } else if("gcv" %in% names(pred_perm$finalModel)){
          pred_perm$finalModel$gcv
        } else if(length(coef(pred_perm$finalModel)) > 0) {
          n = nrow(Response[,j])
          edf = length(colMeans(abs(coef(pred_perm$finalModel))) >= median(abs(coef(pred_perm$finalModel)), na.rm = TRUE))
          rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
          rss / (n * (1 - edf/n)^2)
        } else if(length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall))) > 0){
          edf = length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall)))
          n = nrow(Response[,j])
          rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
          rss / (n * (1 - edf/n)^2)
        } else if("mse" %in% names(pred_perm$finalModel)) {
          mean(pred_perm$finalModel$mse)
        } else {
          n = nrow(Response[,j])
          edf = if(!is.null(VIPS)) {length(VIPS)} else {
            lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
            coef_lasso <- coef(lasso, s = lasso$lambda.min)
            length(which(coef_lasso[rownames(coef_lasso) != "(Intercept)"] != 0))
          }
          rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
          rss / (n * (1 - edf/n)^2)
        }
      }
      if(p_calc == "coef") {
        if("svm" %in% pred_perm$method) {
          place_coef = tryCatch({pred_perm$finalModel@coef},
                                error = function(e){ stats::coef(pred_perm$finalModel)})

          if(is.list(place_coef)){
            names(place_coef) = unique(Response[,j])
            for(i in seq_along(place_coef)){
              place_coef[[i]] <- sum(-1 >= place_coef[[i]] | 1 <= place_coef[[i]])
              coefs[j,] <- as.matrix.data.frame(place_coef) / pred_perm$finalModel@nSV
            }
          }
          coef_id = TRUE

        } else if("K" %in% colnames(pred_perm$results)) {
          named <- colnames(Predictor)

          place_coef <- stats::coef(pred_perm$finalModel)

          coefs[j,] <- place_coef

          coef_id = TRUE
        } else if ("lambda" %in% colnames(pred_perm$bestTune) & !isTRUE(coef_id)) {
          place_coef <- stats::coef(pred_perm$finalModel, s = pred_perm$bestTune$lambda)[-1]
          place_coef[is.na(place_coef) | place_coef == "."] <- 0
          place_coef <- get_coef(place_coef, colnames(Predictor))
          coef_id = TRUE

        } else if("fit" %in% names(pred_perm$finalModel)){
          named <- names(stats::coef(pred_perm$finalModel))[-1]
          coefficients_list = list()
          combined_array <- array(NA, dim = c(length(pred_perm$finalModel$fit), length(colnames(Predictor))))
          colnames(combined_array) <- colnames(Predictor)
          for (i in 1:length(pred_perm$finalModel$fit)) {
            resample_name <- paste0("Resample", sprintf("%02d", i))
            coeff <- t(pred_perm$finalModel$fit[[resample_name]]$coefficients)
            if(length(coeff) > 2) {
              coeff = coeff[1,2:ncol(coeff)]
            } else {
              coeff = coeff[1,]
            }
            coefficients_list[[i]] <- get_coef(coeff, colnames(Predictor))
            combined_array[i, ] <- coefficients_list[[i]]
          }
          place_coef <- apply(combined_array, 2, mean, na.rm = TRUE)
          coefs[j,] <- place_coef
          colnames(coefs) <- colnames(Predictor)
          coef_id = TRUE
        } else if("coefficients" %in% names(pred_perm$finalModel) | "coef" %in% names(pred_perm$finalModel)) {
          place_coef = t(pred_perm$finalModel$coefficients)
          if(length(place_coef) > 2) {
            place_coef = place_coef[1,2:ncol(place_coef)]
          } else {
            place_coef = place_coef[1,]
          }
          place_coef[is.na(place_coef) | place_coef == "."] <- 0
          place_coef <- get_coef(place_coef, colnames(Predictor))
          coefs[j,] <- place_coef
          colnames(coefs) <- colnames(Predictor)
          coef_id = TRUE
        } else if(length(stats::coef(pred_perm$finalModel)) > 0) {
          place_coef <- stats::coef(pred_perm$finalModel)[-1]
          place_coef[is.na(place_coef) | place_coef == "."] <- 0
          place_coef <- get_coef(place_coef, colnames(Predictor))
          coefs[j,] <- place_coef
          colnames(coefs) <- colnames(Predictor)
          coef_id = TRUE
        } else if(length(caret::varImp(pred_perm$finalModel)) > 0){
          place_coef = t(caret::varImp(pred_perm$finalModel))
          colnames(coefs) <- colnames(Predictor)
          coefs[j,colnames(place_coef)] <- place_coef
          coef_id = TRUE
        } else {
          lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
          coef_lasso <- coef(lasso, s = lasso$lambda.min)
          coefs[j,] <- place_coef
          colnames(coefs) <- colnames(Predictor)
          coef_id = TRUE
        }
      }
    } # end of j for loop

    if(type == "class") {
      probs = multivariate_confusion(pred_train_o, Response)[,11]
      probs_real = rep(0.999, times = length(probs))
    }

    if(isTRUE(parallel)) {
      closeAllConnections()
    }
  } # end of full model
  #_________________________________________________________________________________________________________________
  if (isTRUE(print_progress)) {
    pb_start$tick(10, tokens = list(what = "Estimating Important Variables"))
  }
  #Add in Variable Importance in Projections or Variable Regulations_______________________________________________

  if(is.null(VIPS)) { # using Elastic Net Regularization
    VIPS = numeric(0)
    lasso.coefs = list()
    for(j in 1:ncol(Response)) {
      lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
      coef_lasso <- coef(lasso, s = lasso$lambda.min)
      lasso.coefs[[j]] = which(coef_lasso[rownames(coef_lasso) != "(Intercept)"] != 0)
    } # end of j for loop
    lasso.coefs = unique(unlist(lasso.coefs))
    VIPS <- lasso.coefs
    rm(lasso.coefs)
  }
  VIPS = unique(VIPS)
  VIPS = colnames(Predictor[,VIPS])
  pred_reduced = Predictor[,!colnames(Predictor) %in% VIPS]
  pred_reduced = data.frame(pred_reduced)
  colnames(pred_reduced) <- colnames(Predictor)[!colnames(Predictor) %in% VIPS]

  if(ncol(pred_reduced) == 0){
    keep = round(ncol(Predictor) * .20,1)
    if(keep < 3){
      keep = 3
      keep = sample(VIPS,size = keep)
    }
    pred_reduced = Predictor[,c(keep)]
  } else if (ncol(pred_reduced) <= 1){
    keep = 2 - ncol(pred_reduced)
    keep = sample(VIPS,size = keep)
    pred_reduced = Predictor[,c(keep,colnames(pred_reduced))]
  }

  if(type == "class") {
    reduced_fit <- apply(data.frame(matrix(NA, nrow = nrow(Response), ncol = ncol(Response))),2,as.factor)
  } else {
    reduced_fit = matrix(NA, nrow = nrow(Response), ncol = ncol(Response))
  }

  # Conduct a trimmed regression____________________________________________________________________________________________

  if(isFALSE(full_model)) {

    lasso.coefs = list()
    coef_id = FALSE

    models = list()
    reg.parameter = list()

    Predictor = data.frame(Predictor[,VIPS])

    pred_train_o = matrix(NA, ncol = ncol(Response), nrow = nrow(Response))
    coefs = matrix(0, nrow = ncol(Response), ncol = ncol(Predictor))
    place_coef = matrix(NA, nrow = 1, ncol = ncol(Predictor))
    gcv = matrix(NA, nrow = 1, ncol = ncol(Response))

    doParallel::registerDoParallel(cores = core_choice)

    for (j in 1:ncol(Response)) { # loop over each response variable individually

        caret_parameters = modifyList(caret_parameters, list(y = Response[,j], x = as.matrix.data.frame(Predictor)))

        pred_perm = suppressWarnings(do.call(caret::train, caret_parameters))

        pred_train_o[,j] = predict.train(pred_perm, newdata = Predictor)

        if(!is.null(which_classify)) {
          pred_special[,j] = predict.train(pred_perm, newdata = pred_class)
        }

        if(p_calc == "gcv") {

          gcv[1,j] <- if(startsWith(pred_perm$method, "svm")) {
            tryCatch({pred_perm$finalModel@error},
                     error = function(e) {
                       n = nrow(Response[,j])
                       edf = pred_perm$finalModel@nSV
                       rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                       rss / (n * (1 - edf/n)^2)
                     })
          } else if("fit" %in% names(pred_perm$finalModel)) {
            mean(unlist(lapply(pred_perm$finalModel$fit, function(x) x$gcv)))
          } else if("gcv" %in% names(pred_perm$finalModel)){
            pred_perm$finalModel$gcv
          } else if(length(coef(pred_perm$finalModel)) > 0) {
            n = nrow(Response[,j])
            edf = length(colMeans(abs(coef(pred_perm$finalModel))) >= median(abs(coef(pred_perm$finalModel)), na.rm = TRUE))
            rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
            rss / (n * (1 - edf/n)^2)
          } else if(length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall))) > 0){
            edf = length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall)))
            n = nrow(Response[,j])
            rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
            rss / (n * (1 - edf/n)^2)
          } else if("mse" %in% names(pred_perm$finalModel)) {
            mean(pred_perm$finalModel$mse)
          } else {
            n = nrow(Response[,j])
            edf = if(!is.null(VIPS)) {length(VIPS)} else {
              lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
              coef_lasso <- coef(lasso, s = lasso$lambda.min)
              length(which(coef_lasso[rownames(coef_lasso) != "(Intercept)"] != 0))
            }
            rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
            rss / (n * (1 - edf/n)^2)
          }
        }
        if(p_calc == "coef") {
          if("svm" %in% pred_perm$method) {
            place_coef = tryCatch({pred_perm$finalModel@coef},
                                  error = function(e){ stats::coef(pred_perm$finalModel)})

            if(is.list(place_coef)){
              names(place_coef) = unique(Response[,j])
              for(i in seq_along(place_coef)){
                place_coef[[i]] <- sum(-1 >= place_coef[[i]] | 1 <= place_coef[[i]])
                coefs[j,] <- as.matrix.data.frame(place_coef) / pred_perm$finalModel@nSV
              }
            }
            coef_id = TRUE

          } else if("K" %in% colnames(pred_perm$results)) {
            named <- colnames(Predictor)

            place_coef <- stats::coef(pred_perm$finalModel)

            coefs[j,] <- place_coef

            coef_id = TRUE
          } else if ("lambda" %in% colnames(pred_perm$bestTune) & !isTRUE(coef_id)) {
            place_coef <- stats::coef(pred_perm$finalModel, s = pred_perm$bestTune$lambda)[-1]
            place_coef[is.na(place_coef) | place_coef == "."] <- 0
            place_coef <- get_coef(place_coef, colnames(Predictor))
            coef_id = TRUE

          } else if("fit" %in% names(pred_perm$finalModel)){
            named <- names(stats::coef(pred_perm$finalModel))[-1]
            coefficients_list = list()
            combined_array <- array(NA, dim = c(length(pred_perm$finalModel$fit), length(colnames(Predictor))))
            colnames(combined_array) <- colnames(Predictor)
            for (i in 1:length(pred_perm$finalModel$fit)) {
              resample_name <- paste0("Resample", sprintf("%02d", i))
              coeff <- t(pred_perm$finalModel$fit[[resample_name]]$coefficients)
              if(length(coeff) > 2) {
                coeff = coeff[1,2:ncol(coeff)]
              } else {
                coeff = coeff[1,]
              }
              coefficients_list[[i]] <- get_coef(coeff, colnames(Predictor))
              combined_array[i, ] <- coefficients_list[[i]]
            }
            place_coef <- apply(combined_array, 2, mean, na.rm = TRUE)
            coefs[j,] <- place_coef
            colnames(coefs) <- colnames(Predictor)
            coef_id = TRUE
          } else if("coefficients" %in% names(pred_perm$finalModel) | "coef" %in% names(pred_perm$finalModel)) {
            place_coef = t(pred_perm$finalModel$coefficients)
            if(length(place_coef) > 2) {
              place_coef = place_coef[1,2:ncol(place_coef)]
            } else {
              place_coef = place_coef[1,]
            }
            place_coef[is.na(place_coef) | place_coef == "."] <- 0
            place_coef <- get_coef(place_coef, colnames(Predictor))
            coefs[j,] <- place_coef
            colnames(coefs) <- colnames(Predictor)
            coef_id = TRUE
          } else if(length(stats::coef(pred_perm$finalModel)) > 0) {
            place_coef <- stats::coef(pred_perm$finalModel)[-1]
            place_coef[is.na(place_coef) | place_coef == "."] <- 0
            place_coef <- get_coef(place_coef, colnames(Predictor))
            coefs[j,] <- place_coef
            colnames(coefs) <- colnames(Predictor)
            coef_id = TRUE
          } else if(length(caret::varImp(pred_perm$finalModel)) > 0){
            place_coef = t(caret::varImp(pred_perm$finalModel))
            colnames(coefs) <- colnames(Predictor)
            coefs[j,colnames(place_coef)] <- place_coef
            coef_id = TRUE
          } else {
            lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
            coef_lasso <- coef(lasso, s = lasso$lambda.min)
            coefs[j,] <- place_coef
            colnames(coefs) <- colnames(Predictor)
            coef_id = TRUE
          }
        }
    } # j loop
    pred_perm = NULL

    if(type == "class") {
      probs = multivariate_confusion(pred_train_o, Response)[,11]
    }

    if(isTRUE(parallel)) {
      closeAllConnections()
    }
  }

  if(ncol(Predictor) >= nrow(Predictor) || ncol(Predictor) >= nrow(Predictor) * 0.9) {
    lasso.coefs = list()
    for(j in 1:ncol(Response)) {
      lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
      coef_lasso <- coef(lasso, s = lasso$lambda.min)
      lasso.coefs[[j]] = which(coef_lasso[rownames(coef_lasso) != "(Intercept)"] != 0)
      reg.parameter[[j]] <- length(lasso.coefs[[j]])
    }
    reg.parameter = round(mean(unlist(reg.parameter)))

    p_test = reg.parameter
    n_test = nrow(Response)
    print("Too many Predictor variables detected. Regularization will be applied.")
  } else {
    p_test = ncol(Predictor)
    n_test = nrow(Response)
  }


  # Permutation Preliminary Convergence Test________________________________________________________________________

  if (isTRUE(print_progress)) {
    pb_start$tick(10, tokens = list(what = "Preliminary Convergence Test"))
  }

  doParallel::registerDoParallel(cores = core_choice)

  converge_perm = 0.1 * perm

  for(j in 1:ncol(Response)) {

      reduced_caret = modifyList(caret_parameters, list(y = Response[,j], x = pred_reduced))

      reduced_model <- tryCatch({suppressWarnings(do.call(caret::train, reduced_caret))}, error = function (e) {
        print("Invalid parameters for Caret package. Please consider changing the train control parameters.")
      })

      reduced_fit[,j] <- predict.train(reduced_model)

  }

  # Preliminary test for permutation convergence ____________________________________________________________________________________

  prelim_perm <- foreach(1:converge_perm, .packages = c("e1071", "glmnet", "caret", "doParallel", "parallel", "kernlab", "kerndwd", "tidyverse", "randomForest", "CSGM")) %dopar% {

    samp_test <- sample(1:nrow(reduced_fit), replace=FALSE)

    if(type == "class") {
      perm_test_y = data.frame(reduced_fit[samp_test,])

      probs_random = multivariate_confusion(perm_test_y, Response)[,11]

    } else {
      residuals_test <- data.frame(Response - reduced_fit)[samp_test,]
      if(is.null(ncol(residuals_test))) {
        residuals_test = data.frame(Response = residuals_test)
      }
      rownames(residuals_test) <- 1:nrow(residuals_test)
      perm_test_y <- reduced_fit + residuals_test
      perm_test_y = as.matrix.data.frame(perm_test_y)
    }
    pred_test = matrix(NA, nrow = nrow(perm_test_y), ncol = ncol(perm_test_y))
    Output = numeric(0)

    perm_gcv = matrix(NA, nrow = 1, ncol = ncol(Response))
    place_coef = matrix(NA, nrow = 1, ncol = ncol(Predictor))
    perm_coefs = matrix(NA, nrow = ncol(Response), ncol = ncol(Predictor))
    colnames(perm_coefs) <- colnames(Predictor)

    for(j in 1:ncol(perm_test_y)) {

      pred_perm <- tryCatch({suppressWarnings(caret::train(y = perm_test_y[,j], x = data.frame(Predictor),
                                                          method = caret_parameters$method, trControl = trainControl(method = "none")))
      }, error = function(e) {
        suppressWarnings(caret::train(y = perm_test_y[,j], x = data.frame(Predictor),
                                      method = caret_parameters$method, trControl = trainControl(method = "cv", number = 5, repeats = 1)))
      })
      pred_test[,j] <- predict.train(pred_perm)

      if(p_calc == "gcv") {

        perm_gcv[1,j] <- if(startsWith(pred_perm$method, "svm")) {
          tryCatch({pred_perm$finalModel@error},
                   error = function(e) {
                     n = nrow(Response[,j])
                     edf = pred_perm$finalModel@nSV
                     rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                     rss / (n * (1 - edf/n)^2)
                   })
        } else if("fit" %in% names(pred_perm$finalModel)) {
          mean(unlist(lapply(pred_perm$finalModel$fit, function(x) x$gcv)))
        } else if("gcv" %in% names(pred_perm$finalModel)){
          pred_perm$finalModel$gcv
        } else if(length(coef(pred_perm$finalModel)) > 0) {
          n = nrow(Response[,j])
          edf = length(colMeans(abs(coef(pred_perm$finalModel))) >= median(abs(coef(pred_perm$finalModel)), na.rm = TRUE))
          rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
          rss / (n * (1 - edf/n)^2)
        } else if(length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall))) > 0){
          edf = length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall)))
          n = nrow(Response[,j])
          rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
          rss / (n * (1 - edf/n)^2)
        } else if("mse" %in% names(pred_perm$finalModel)) {
          mean(pred_perm$finalModel$mse)
        } else {
          n = nrow(Response[,j])
          edf = if(!is.null(VIPS)) {length(VIPS)} else {
            lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
            coef_lasso <- coef(lasso, s = lasso$lambda.min)
            length(which(coef_lasso[rownames(coef_lasso) != "(Intercept)"] != 0))
          }
          rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
          rss / (n * (1 - edf/n)^2)
        }
      }
      if(p_calc == "coef") {
        if("svm" %in% pred_perm$method) {
          place_coef = tryCatch({pred_perm$finalModel@coef},
                                error = function(e){ stats::coef(pred_perm$finalModel)})

          if(is.list(place_coef)){
            names(place_coef) = unique(Response[,j])
            for(i in seq_along(place_coef)){
              place_coef[[i]] <- sum(-1 >= place_coef[[i]] | 1 <= place_coef[[i]])
              perm_coefs[j,] <- as.matrix.data.frame(place_coef) / pred_perm$finalModel@nSV
            }
          }
          coef_id = TRUE

        } else if("K" %in% colnames(pred_perm$results)) {
          named <- colnames(Predictor)

          place_coef <- stats::coef(pred_perm$finalModel)

          perm_coefs[j,] <- place_coef

          coef_id = TRUE
        } else if ("lambda" %in% colnames(pred_perm$bestTune) & !isTRUE(coef_id)) {
          place_coef <- stats::coef(pred_perm$finalModel, s = pred_perm$bestTune$lambda)[-1]
          place_coef[is.na(place_coef) | place_coef == "."] <- 0
          place_coef <- get_coef(place_coef, colnames(Predictor))
          coef_id = TRUE

        } else if("fit" %in% names(pred_perm$finalModel)){
          named <- names(stats::coef(pred_perm$finalModel))[-1]
          coefficients_list = list()
          combined_array <- array(NA, dim = c(length(pred_perm$finalModel$fit), length(colnames(Predictor))))
          colnames(combined_array) <- colnames(Predictor)
          for (i in 1:length(pred_perm$finalModel$fit)) {
            resample_name <- paste0("Resample", sprintf("%02d", i))
            coeff <- t(pred_perm$finalModel$fit[[resample_name]]$coefficients)
            if(length(coeff) > 2) {
              coeff = coeff[1,2:ncol(coeff)]
            } else {
              coeff = coeff[1,]
            }
            coefficients_list[[i]] <- get_coef(coeff, colnames(Predictor))
            combined_array[i, ] <- coefficients_list[[i]]
          }
          place_coef <- apply(combined_array, 2, mean, na.rm = TRUE)
          perm_coefs[j,] <- place_coef
          colnames(perm_coefs) <- colnames(Predictor)
          coef_id = TRUE
        } else if("coefficients" %in% names(pred_perm$finalModel) | "coef" %in% names(pred_perm$finalModel)) {
          place_coef = t(pred_perm$finalModel$coefficients)
          if(length(place_coef) > 2) {
            place_coef = place_coef[1,2:ncol(place_coef)]
          } else {
            place_coef = place_coef[1,]
          }
          place_coef[is.na(place_coef) | place_coef == "."] <- 0
          place_coef <- get_coef(place_coef, colnames(Predictor))
          perm_coefs[j,] <- place_coef
          colnames(perm_coefs) <- colnames(Predictor)
          coef_id = TRUE
        } else if(length(stats::coef(pred_perm$finalModel)) > 0) {
          place_coef <- stats::coef(pred_perm$finalModel)[-1]
          place_coef[is.na(place_coef) | place_coef == "."] <- 0
          place_coef <- get_coef(place_coef, colnames(Predictor))
          perm_coefs[j,] <- place_coef
          colnames(perm_coefs) <- colnames(Predictor)
          coef_id = TRUE
        } else if(length(caret::varImp(pred_perm$finalModel)) > 0){
          place_coef = t(caret::varImp(pred_perm$finalModel))
          colnames(perm_coefs) <- colnames(Predictor)
          perm_coefs[j,colnames(place_coef)] <- place_coef
          coef_id = TRUE
        } else {
          lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
          coef_lasso <- coef(lasso, s = lasso$lambda.min)
          perm_coefs[j,] <- place_coef
          colnames(perm_coefs) <- colnames(Predictor)
          coef_id = TRUE
        }
      }
    } # end of j for loop

    if(type == "class") {
      probs_perm = multivariate_confusion(pred_test, perm_test_y)[,11]
    }

    if(p_calc == "coef") {
      Output <- perm_coefs
      rm(perm_coefs, place_coef)
    }
    if(p_calc == "gcv") {
      Output <- mean(perm_gcv, na.rm = TRUE)
    }
    if(p_calc == "F") {

      if(type == "class") {
        SSR_perm <- sum((probs_random - mean(probs_perm))^2)
        SSE_perm <- sum((probs_random - probs_perm)^2)
        MSR <- as.numeric(SSR_perm / p_test)
        MSE <- as.numeric(SSE_perm / (n_test - p_test))
        F_test <- MSR/MSE
      } else {
        SSR_perm <- sum((pred_test - colMeans(perm_test_y))^2)
        SSE_perm <- sum((perm_test_y - pred_test)^2)
        MSR <- as.numeric(SSR_perm / p_test)
        MSE <- as.numeric(SSE_perm / (n_test - p_test))
        F_test <- MSR/MSE
      }
      Output <- F_test
    }
    if(p_calc == "RMSE") {
      if(type == "class") {
        RMSE <- sqrt(mean((probs_random - probs_perm)^2))
      } else {
        RMSE <- sqrt(mean((perm_test_y - pred_test)^2))
      }
      Output <- RMSE
    }

    if(p_calc == "R2") {
      if(type == "class") {
        R2 <- m_R2(probs_random, probs_perm)
      } else {
        R2 <- m_R2(perm_test_y, pred_test)
      }
      Output <- R2
    }

    pred_test = NULL

    return(Output)
  }

  convergence = FALSE

  if(p_calc == "coef") {
    perm_coefs <- prelim_perm
    array_coefs <- simplify2array(perm_coefs)
    sd_coefs <- apply(array_coefs, c(1, 2), sd)
    rm(array_coefs)
    se <- median(sd_coefs, na.rm = TRUE)
  } else {
    Output <- unlist(prelim_perm)
    se <- mean(sd(Output, na.rm = TRUE), na.rm = TRUE)
  }

  if(p_calc == "F") {

    if(round(se,3) <= 2) {
      perm <- converge_perm
      convergence <- TRUE
      print(paste("Early Convergence of the Permutations Reached. Permutation Test set to:", perm))

    } else {
      print(paste("Early Convergence of the Permutations Not Reached. Full Permutation Test Will Be Computed."))
      perm = perm - converge_perm
    }
  } else {

    if(round(se,3) <= 0.05) {
      perm <- converge_perm
      convergence <- TRUE
      print(paste("Early Convergence of the Permutations Reached. Permutation Test set to:", perm))

    } else {
      print(paste("Early Convergence of the Permutations Not Reached. Full Permutation Test Will Be Computed."))
      perm = perm - converge_perm
    }
  }

  rm(prelim_perm)

  if(isTRUE(parallel)) {
    closeAllConnections()
  }






  if (isTRUE(print_progress)) {
    pb_start$tick(10, tokens = list(what = paste("Regression Test:")))
  }


  # ______________________________________________________________________________________________________________

  # Conduct bootstrapping and/or cross-validation_________________________________________________________________________
  if (n_boot >= 1) {

    if(isTRUE(parallel)) {
      doParallel::registerDoParallel(cores = core_choice)
    } else {
      foreach::registerDoSEQ()
    }

    i = 1

    boot_iter = foreach::foreach(i = 1:n_boot, .packages =
                                   c("e1071", "glmnet", "caret", "plyr", "doParallel", "kernlab", "kerndwd", "tidyverse", "randomForest", "CSGM")) %dopar% {

       Prediction = matrix(NA, nrow = nrow(Response), ncol = ncol(Response))
       gcv_test2 = list()
       coefs_test2 = list()

       if(isTRUE(CV)) {

         if(type == "class") {
           cv_folds <- caret::createMultiFolds(as.factor(Response[[1]]), k = 5, times = 1)
         } else {
           cv_folds <- caret::createFolds(y = 1:nrow(Response), k = 5, list = TRUE, returnTrain = TRUE)
         }

         RMSE_test2 = vector("numeric", 5)
         SSR_test2 = vector("numeric", 5)
         SSE_test2 = vector("numeric", 5)
         SST_test2 = vector("numeric", 5)
         R2_test2 = vector("numeric", 5)
         MSR = vector("numeric", 5)
         MSE = vector("numeric", 5)
         F_test2 = vector("numeric", 5)
         BIC_test2 = vector("numeric", 5)
         if(!is.null(which_classify)) {
           pred_special_cv = list(Fold_1 = matrix(nrow = nrow(res_class), ncol = ncol(res_class)),
                                  Fold_2 = matrix(nrow = nrow(res_class), ncol = ncol(res_class)),
                                  Fold_3 = matrix(nrow = nrow(res_class), ncol = ncol(res_class)),
                                  Fold_4 = matrix(nrow = nrow(res_class), ncol = ncol(res_class)),
                                  Fold_5 = matrix(nrow = nrow(res_class), ncol = ncol(res_class)))
         }
         # Conduct the folds

         for (fold in 1:5) {
           train_indices <- cv_folds[[fold]]
           test_indices <- setdiff(1:nrow(Response), train_indices)

           X_train <- as.matrix.data.frame(Predictor[train_indices, ])
           Y_train <- as.matrix.data.frame(Response[train_indices, ])
           X_test <- as.matrix.data.frame(Predictor[test_indices, ])
           Y_test <- as.matrix.data.frame(Response[test_indices, ])
           pred_test = matrix(NA, nrow(Y_test), ncol = ncol(Y_test))
           pred_train = matrix(NA, nrow(Y_train), ncol = ncol(Y_train))

           n_test = nrow(X_test)
           if(nrow(X_test) <= ncol(X_test)) {
             p_test <- length(VIPS)
             if(p_test == 0) {
               p_test <- ncol(X_test)
             }
           } else {
             p_test <- ncol(X_test)
           }

           gcv = matrix(NA, nrow = 1, ncol = ncol(Y_test))
           coefs = matrix(NA, nrow = ncol(Y_test), ncol = ncol(X_test))

           for(j in 1:ncol(Y_train)) {

               caret_parameters_boot <- modifyList(caret_parameters, list(x = X_train, y  = Y_train[,j], trControl = trainControl(method = "none")))
               pred_perm <- suppressWarnings(do.call(caret::train, caret_parameters_boot))
               pred_test[,j] <- predict.train(pred_perm, X_test)
               pred_train[,j] <- predict.train(pred_perm, X_train)

               if(p_calc == "gcv") {

                 gcv[1,j] <- if(startsWith(pred_perm$method, "svm")) {
                   tryCatch({pred_perm$finalModel@error},
                            error = function(e) {
                              n = nrow(Response[,j])
                              edf = pred_perm$finalModel@nSV
                              rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                              rss / (n * (1 - edf/n)^2)
                            })
                 } else if("fit" %in% names(pred_perm$finalModel)) {
                   mean(unlist(lapply(pred_perm$finalModel$fit, function(x) x$gcv)))
                 } else if("gcv" %in% names(pred_perm$finalModel)){
                   pred_perm$finalModel$gcv
                 } else if(length(coef(pred_perm$finalModel)) > 0) {
                   n = nrow(Response[,j])
                   edf = length(colMeans(abs(coef(pred_perm$finalModel))) >= median(abs(coef(pred_perm$finalModel)), na.rm = TRUE))
                   rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                   rss / (n * (1 - edf/n)^2)
                 } else if(length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall))) > 0){
                   edf = length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall)))
                   n = nrow(Response[,j])
                   rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                   rss / (n * (1 - edf/n)^2)
                 } else if("mse" %in% names(pred_perm$finalModel)) {
                   mean(pred_perm$finalModel$mse)
                 } else {
                   n = nrow(Response[,j])
                   edf = if(!is.null(VIPS)) {length(VIPS)} else {
                     lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
                     coef_lasso <- coef(lasso, s = lasso$lambda.min)
                     length(which(coef_lasso[rownames(coef_lasso) != "(Intercept)"] != 0))
                   }
                   rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                   rss / (n * (1 - edf/n)^2)
                 }
               }
               if(p_calc == "coef") {
                 if("svm" %in% pred_perm$method) {
                   place_coef = tryCatch({pred_perm$finalModel@coef},
                                         error = function(e){ stats::coef(pred_perm$finalModel)})

                   if(is.list(place_coef)){
                     names(place_coef) = unique(Response[,j])
                     for(i in seq_along(place_coef)){
                       place_coef[[i]] <- sum(-1 >= place_coef[[i]] | 1 <= place_coef[[i]])
                       coefs[j,] <- as.matrix.data.frame(place_coef) / pred_perm$finalModel@nSV
                     }
                   }
                   coef_id = TRUE

                 } else if("K" %in% colnames(pred_perm$results)) {
                   named <- colnames(Predictor)

                   place_coef <- stats::coef(pred_perm$finalModel)

                   coefs[j,] <- place_coef

                   coef_id = TRUE
                 } else if ("lambda" %in% colnames(pred_perm$bestTune) & !isTRUE(coef_id)) {
                   place_coef <- stats::coef(pred_perm$finalModel, s = pred_perm$bestTune$lambda)[-1]
                   place_coef[is.na(place_coef) | place_coef == "."] <- 0
                   place_coef <- get_coef(place_coef, colnames(Predictor))
                   coef_id = TRUE

                 } else if("fit" %in% names(pred_perm$finalModel)){
                   named <- names(stats::coef(pred_perm$finalModel))[-1]
                   coefficients_list = list()
                   combined_array <- array(NA, dim = c(length(pred_perm$finalModel$fit), length(colnames(Predictor))))
                   colnames(combined_array) <- colnames(Predictor)
                   for (i in 1:length(pred_perm$finalModel$fit)) {
                     resample_name <- paste0("Resample", sprintf("%02d", i))
                     coeff <- t(pred_perm$finalModel$fit[[resample_name]]$coefficients)
                     if(length(coeff) > 2) {
                       coeff = coeff[1,2:ncol(coeff)]
                     } else {
                       coeff = coeff[1,]
                     }
                     coefficients_list[[i]] <- get_coef(coeff, colnames(Predictor))
                     combined_array[i, ] <- coefficients_list[[i]]
                   }
                   place_coef <- apply(combined_array, 2, mean, na.rm = TRUE)
                   coefs[j,] <- place_coef
                   colnames(coefs) <- colnames(Predictor)
                   coef_id = TRUE
                 } else if("coefficients" %in% names(pred_perm$finalModel) | "coef" %in% names(pred_perm$finalModel)) {
                   place_coef = t(pred_perm$finalModel$coefficients)
                   if(length(place_coef) > 2) {
                     place_coef = place_coef[1,2:ncol(place_coef)]
                   } else {
                     place_coef = place_coef[1,]
                   }
                   place_coef[is.na(place_coef) | place_coef == "."] <- 0
                   place_coef <- get_coef(place_coef, colnames(Predictor))
                   coefs[j,] <- place_coef
                   colnames(coefs) <- colnames(Predictor)
                   coef_id = TRUE
                 } else if(length(stats::coef(pred_perm$finalModel)) > 0) {
                   place_coef <- stats::coef(pred_perm$finalModel)[-1]
                   place_coef[is.na(place_coef) | place_coef == "."] <- 0
                   place_coef <- get_coef(place_coef, colnames(Predictor))
                   coefs[j,] <- place_coef
                   colnames(coefs) <- colnames(Predictor)
                   coef_id = TRUE
                 } else if(length(caret::varImp(pred_perm$finalModel)) > 0){
                   place_coef = t(caret::varImp(pred_perm$finalModel))
                   colnames(coefs) <- colnames(Predictor)
                   coefs[j,colnames(place_coef)] <- place_coef
                   coef_id = TRUE
                 } else {
                   lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
                   coef_lasso <- coef(lasso, s = lasso$lambda.min)
                   coefs[j,] <- place_coef
                   colnames(coefs) <- colnames(Predictor)
                   coef_id = TRUE
                 }
               }
             if(type == "class") {
               if(!is.null(which_classify)) {
                 pred_special_cv[[fold]][,j] = predict(pred_perm, newdata = pred_class)
               } # end which_classify loop
             } # end class loop
           } # end of j for loop

           if(type == "class") {

             Prediction[test_indices,] <- pred_test

             probs_test <- list()
             probs_train <- list()

             probs_test <- multivariate_confusion(pred_test, Y_test)[,11]
             probs_train <- multivariate_confusion(pred_train, Y_train)[,11]

             R2_test2[fold] <- m_R2(probs_train, probs_test)
             RMSE_test2[fold] <- sqrt(mean((probs_train - probs_test)^2))
             SSR_test2[fold] <- sum((probs_test - mean(probs_train))^2)
             SSE_test2[fold] <- sum((probs_train - probs_test)^2)
           } else {
             RMSE_test2[fold] <- sqrt(mean((as.matrix(Y_test) - as.matrix(pred_test))^2))
             R2_test2[fold] <- m_R2(Y_test, pred_test)
             SSR_test2[fold] <- sum((pred_test - colMeans(Y_test))^2)
             SSE_test2[fold] <- sum((Y_test - pred_test)^2)
           }


           MSR[fold] = as.numeric(SSR_test2[fold] / p_test)
           MSE[fold] = as.numeric(SSE_test2[fold] / abs(n_test - p_test))
           F_test2[fold] <- MSR[fold] / MSE[fold]
           BIC_test2[fold] <- abs((n_test * log(SSE_test2[fold]/n_test)) + (log(n_test) * p_test))
           if(p_calc == "gcv") {
             gcv_test2[[fold]] <- colMeans(gcv, na.rm = TRUE)
           }
           if(p_calc == "coef"){
             coefs_test2[[fold]] <- colMeans(coefs, na.rm = TRUE)
           }
         } # end of cross-validation

         if(type == "class") {
           if(!is.null(which_classify)) {
             pred_special_cv <- data.frame(do.call(rbind, pred_special_cv) %>%
                                             sapply(., function(x) {
                                               pred_matrix <- matrix(x, ncol=ncol(res_class), nrow = nrow(res_class))
                                               apply(pred_matrix, 2, function(y) names(which.max(table(y))))
                                             }))
           }
         } else {
           pred_special_cv = data.frame(Response = c(0,0,0,0,0))
         }

         RMSE_test2 <- min(RMSE_test2, na.rm = TRUE)
         R2_test2  <- max(R2_test2, na.rm = TRUE)
         F_test2 <- mean(F_test2, na.rm = TRUE)
         BIC_test2 <- min(abs(BIC_test2), na.rm = TRUE)
         if(p_calc == "gcv"){
           gcv_test2 <- colMeans(cbind(unlist(gcv_test2)), na.rm = TRUE)
         } else {
           gcv_test2 = 0
         }
         if(p_calc == "coef"){
           coefs_test2 <- colMeans(cbind(unlist(coefs_test2)), na.rm = TRUE)
         } else {
           coefs_test2 = 0
         }


      #__________________________________________________________________________
      #__________________________________________________________________________
       } else { # for just bootstrapping


         index <- caret::createDataPartition(p = .70, y = Response[,1], list = FALSE)
         test_x <- as.matrix.data.frame(Predictor[-index,])
         test_y <- as.matrix.data.frame(Response[-index,])
         train_x <- as.matrix.data.frame(Predictor[index,])
         train_y <- as.matrix.data.frame(Response[index,])

         n_test = nrow(test_x)
         p_test = ncol(test_x)

         if(ncol(test_x) >= nrow(test_x) || ncol(test_x)/nrow(test_x) >= 0.80) {
           p_test = length(VIPS)
         }

         if(!is.null(which_classify)) {
           pred_special_cv = matrix(nrow = nrow(res_class), ncol = ncol(res_class))
         }

         pred_train <- matrix(nrow = nrow(train_y), ncol = ncol(train_y))
         pred_test <- matrix(nrow = nrow(test_y), ncol = ncol(test_y))


         for(j in 1:ncol(train_y)) {

             caret_parameters_boot <- modifyList(caret_parameters, list(x = train_x, y  = train_y[,1], trControl = trainControl(method = "none")))

             pred_perm <- suppressWarnings(do.call(caret::train, caret_parameters_boot))

             pred_test[,j] <- predict.train(pred_perm, test_x)

           if(type == "class") {
             if(!is.null(which_classify)) {
               pred_special_cv[,j] = predict(pred_perm, newdata = pred_class)
             } # end which_classify loop
           } # end class loop

             if(p_calc == "gcv") {

               gcv[1,j] <- if(startsWith(pred_perm$method, "svm")) {
                 tryCatch({pred_perm$finalModel@error},
                          error = function(e) {
                            n = nrow(Response[,j])
                            edf = pred_perm$finalModel@nSV
                            rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                            rss / (n * (1 - edf/n)^2)
                          })
               } else if("fit" %in% names(pred_perm$finalModel)) {
                 mean(unlist(lapply(pred_perm$finalModel$fit, function(x) x$gcv)))
               } else if("gcv" %in% names(pred_perm$finalModel)){
                 pred_perm$finalModel$gcv
               } else if(length(coef(pred_perm$finalModel)) > 0) {
                 n = nrow(Response[,j])
                 edf = length(colMeans(abs(coef(pred_perm$finalModel))) >= median(abs(coef(pred_perm$finalModel)), na.rm = TRUE))
                 rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                 rss / (n * (1 - edf/n)^2)
               } else if(length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall))) > 0){
                 edf = length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall)))
                 n = nrow(Response[,j])
                 rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                 rss / (n * (1 - edf/n)^2)
               } else if("mse" %in% names(pred_perm$finalModel)) {
                 mean(pred_perm$finalModel$mse)
               } else {
                 n = nrow(Response[,j])
                 edf = if(!is.null(VIPS)) {length(VIPS)} else {
                   lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
                   coef_lasso <- coef(lasso, s = lasso$lambda.min)
                   length(which(coef_lasso[rownames(coef_lasso) != "(Intercept)"] != 0))
                 }
                 rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                 rss / (n * (1 - edf/n)^2)
               }
             }
             if(p_calc == "coef") {
               if("svm" %in% pred_perm$method) {
                 place_coef = tryCatch({pred_perm$finalModel@coef},
                                       error = function(e){ stats::coef(pred_perm$finalModel)})

                 if(is.list(place_coef)){
                   names(place_coef) = unique(Response[,j])
                   for(i in seq_along(place_coef)){
                     place_coef[[i]] <- sum(-1 >= place_coef[[i]] | 1 <= place_coef[[i]])
                     coefs[j,] <- as.matrix.data.frame(place_coef) / pred_perm$finalModel@nSV
                   }
                 }
                 coef_id = TRUE

               } else if("K" %in% colnames(pred_perm$results)) {
                 named <- colnames(Predictor)

                 place_coef <- stats::coef(pred_perm$finalModel)

                 coefs[j,] <- place_coef

                 coef_id = TRUE
               } else if ("lambda" %in% colnames(pred_perm$bestTune) & !isTRUE(coef_id)) {
                 place_coef <- stats::coef(pred_perm$finalModel, s = pred_perm$bestTune$lambda)[-1]
                 place_coef[is.na(place_coef) | place_coef == "."] <- 0
                 place_coef <- get_coef(place_coef, colnames(Predictor))
                 coef_id = TRUE

               } else if("fit" %in% names(pred_perm$finalModel)){
                 named <- names(stats::coef(pred_perm$finalModel))[-1]
                 coefficients_list = list()
                 combined_array <- array(NA, dim = c(length(pred_perm$finalModel$fit), length(colnames(Predictor))))
                 colnames(combined_array) <- colnames(Predictor)
                 for (i in 1:length(pred_perm$finalModel$fit)) {
                   resample_name <- paste0("Resample", sprintf("%02d", i))
                   coeff <- t(pred_perm$finalModel$fit[[resample_name]]$coefficients)
                   if(length(coeff) > 2) {
                     coeff = coeff[1,2:ncol(coeff)]
                   } else {
                     coeff = coeff[1,]
                   }
                   coefficients_list[[i]] <- get_coef(coeff, colnames(Predictor))
                   combined_array[i, ] <- coefficients_list[[i]]
                 }
                 place_coef <- apply(combined_array, 2, mean, na.rm = TRUE)
                 coefs[j,] <- place_coef
                 colnames(coefs) <- colnames(Predictor)
                 coef_id = TRUE
               } else if("coefficients" %in% names(pred_perm$finalModel) | "coef" %in% names(pred_perm$finalModel)) {
                 place_coef = t(pred_perm$finalModel$coefficients)
                 if(length(place_coef) > 2) {
                   place_coef = place_coef[1,2:ncol(place_coef)]
                 } else {
                   place_coef = place_coef[1,]
                 }
                 place_coef[is.na(place_coef) | place_coef == "."] <- 0
                 place_coef <- get_coef(place_coef, colnames(Predictor))
                 coefs[j,] <- place_coef
                 colnames(coefs) <- colnames(Predictor)
                 coef_id = TRUE
               } else if(length(stats::coef(pred_perm$finalModel)) > 0) {
                 place_coef <- stats::coef(pred_perm$finalModel)[-1]
                 place_coef[is.na(place_coef) | place_coef == "."] <- 0
                 place_coef <- get_coef(place_coef, colnames(Predictor))
                 coefs[j,] <- place_coef
                 colnames(coefs) <- colnames(Predictor)
                 coef_id = TRUE
               } else if(length(caret::varImp(pred_perm$finalModel)) > 0){
                 place_coef = t(caret::varImp(pred_perm$finalModel))
                 colnames(coefs) <- colnames(Predictor)
                 coefs[j,colnames(place_coef)] <- place_coef
                 coef_id = TRUE
               } else {
                 lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
                 coef_lasso <- coef(lasso, s = lasso$lambda.min)
                 coefs[j,] <- place_coef
                 colnames(coefs) <- colnames(Predictor)
                 coef_id = TRUE
               }
             }
        } # end of j for loop

         Prediction[index,] <- pred_train
         Prediction[-index,] <- pred_test

         if(type == "class") {

           probs_train <- list()
           probs_test <-list()
           probs_test <- multivariate_confusion(pred_test, test_y)[,11]
           probs_train <- multivariate_confusion(pred_train, train_y)[,11]

           RMSE_test2 <- sqrt(mean((probs_train - probs_test)^2))
           SSR_test2 <- sum((probs_test - mean(probs_train))^2)
           SSE_test2 <- sum((probs_train - probs_test)^2)
           SST_test2 <- sum((probs_train - mean(probs_train))^2)

         } else {

           RMSE_test2 <- sqrt(mean((as.matrix(test_y) - as.matrix(pred_test))^2))
           SSR_test2 <- sum((pred_test - colMeans(test_y))^2)
           SSE_test2 <- sum((test_y - pred_test)^2)
           SST_test2 <- sum((test_y - colMeans(test_y))^2)
         }

         R2_test2 <- m_R2(test_y, pred_test)
         MSR = as.numeric(SSR_test2 / p_test)
         MSE = as.numeric(SSE_test2 / (n_test - p_test))
         F_test2 <- MSR / MSE
         BIC_test2 <- abs((n_test * log(SSE_test2/n_test)) + (log(n_test) * p_test))
         if(p_calc == "gcv") {
           gcv_test2 <- min(gcv_test2, na.rm = TRUE)
         } else {
           gcv_test2 = 0
         }
         if(p_calc == "coef"){
           coefs_test2 = colMeans(abs(coefs), na.rm = TRUE)
         } else {
           coefs_test2 = 0
         }
       } # end of CV or boot
       return(list(R2_test2 = R2_test2,
                   RMSE_test2 = RMSE_test2,
                   BIC_test2 = BIC_test2,
                   F_test2 = F_test2,
                   Prediction = Prediction,
                   pred_special_cv = pred_special_cv,
                   gcv_test = gcv_test2,
                   coefs_test2 = coefs_test2))

     } # end of boot_iter for loop

     if(isTRUE(parallel)) {
        closeAllConnections()
     }

    if(type == "class") {

      var_boot = c("R2_test2","RMSE_test2",
                   "BIC_test2", "F_test2",
                   "Prediction", "pred_special_cv",
                   "gcv_test2", "coefs_test2")

      for (stat_name in var_boot) {
        assign(stat_name, lapply(boot_iter, `[[`, stat_name))
      }

      Prediction <- data.frame(do.call(rbind, Prediction) %>%
                                 sapply(., function(x) {
                                   pred_matrix <- matrix(x, ncol=ncol(Response), nrow = nrow(Response))
                                   apply(pred_matrix, 2, function(y) names(which.max(table(y))))
                                 }))
      if(p_calc == "coef"){
        coefs <- data.frame(do.call(rbind, coefs_test2) %>%
                              sapply(., function(x) {
                                pred_matrix <- matrix(x, ncol=ncol(Predictor), nrow = ncol(Response))
                                apply(pred_matrix, 2, function(y) names(which.max(table(y))))
                              }))
      }


      summary = multivariate_confusion(Prediction, Response)


      if(!is.null(which_classify)) {

        pred_special <- data.frame(do.call(rbind, pred_special_cv) %>%
                                     sapply(., function(x) {
                                       pred_matrix <- matrix(x, ncol=ncol(res_class), nrow = nrow(res_class))
                                       apply(pred_matrix, 2, function(y) names(which.max(table(y))))
                                     }))
        colnames(Prediction) <- colnames(Response)
        colnames(pred_special) <- colnames(Response)
        Prediction <- rbind(pred_special, Prediction)
      }

    } else {

      var_boot = c("R2_test2","RMSE_test2",
                   "BIC_test2", "F_test2",
                   "Prediction", "pred_special_cv",
                   "gcv_test", "coefs_test2")

      for (stat_name in var_boot) {
        assign(stat_name, lapply(boot_iter, `[[`, stat_name))
      }

      Prediction <- matrix(colMeans(do.call(rbind, Prediction), na.rm = TRUE), nrow = nrow(Response), ncol = ncol(Response))
      if(p_calc == "coef"){
        coefs <- colMeans(rbind(unlist(coefs_test2)), na.rm = TRUE)
      }
      if(p_calc == "gcv"){
        gcv <- mean(unlist(coefs_test2), na.rm = TRUE)
      }
      RMSE_test <- mean(unlist(RMSE_test2), na.rm = TRUE)
      R2_test <- mean(unlist(R2_test2), na.rm = TRUE)
      BIC_test <- mean(unlist(BIC_test2), na.rm = TRUE)
      F_test <- mean(unlist(F_test2), na.rm = TRUE)

      #CIU_test_RMSE <-  mean(quantile(unlist(RMSE_test), probs = 0.975, na.rm = TRUE), na.rm = TRUE)
      #CIL_test_RMSE <-  mean(quantile(unlist(RMSE_test), probs = 0.025, na.rm = TRUE), na.rm = TRUE)
    }

  } else { #full model

    n_test <- nrow(Predictor)

    Prediction = data.frame(pred_train_o)

    if(type == "class") {

      summary = multivariate_confusion(pred_train_o, Response)

      if(!is.null(which_classify)) {
        colnames(pred_special) <- colnames(Response)
        colnames(pred_train_o) <- colnames(Response)
        Prediction <- rbind(pred_special, pred_train_o)
      }

    } else {

      RMSE_test <- sqrt(mean((as.matrix(Response) - as.matrix(pred_train_o))^2))

      R2_test  <- m_R2(Response, pred_train_o)

      SSR_test <- sum((pred_train_o - colMeans(Response))^2)
      SSE_test <- sum((Response - pred_train_o)^2)

      MSR = as.numeric(SSR_test / p_test)
      MSE = as.numeric(SSE_test / (n_test - p_test))
      F_test <- MSR / MSE
      BIC_test <- abs((n_test * log(SSE_test/n_test)) + (p_test * log(n_test)))

      #CIU_test_RMSE <- quantile(RMSE_test, probs = 0.975)
      #CIL_test_RMSE <- quantile(RMSE_test, probs = 0.025)
      if(p_calc == "gcv"){
        gcv_test <- min(gcv, na.rm = TRUE)
      }
    }

  } # end of if full model

  #Permutation testing_____________________________________________________________________________________________

  if (isTRUE(print_progress)) {
    pb_start$tick(15, tokens = list(what = paste("RRPP Permutations:")))
  }

  if(isFALSE(convergence)) {

    if(isTRUE(parallel)) {
      doParallel::registerDoParallel(cores = core_choice)
    } else {
      foreach::registerDoSEQ()
    }

    perm_iter <- foreach(1:perm, .packages = c("e1071", "glmnet", "caret", "doParallel", "parallel", "kernlab", "kerndwd", "tidyverse", "randomForest", "CSGM")) %dopar% {

      samp_test <- sample(1:nrow(reduced_fit), replace=FALSE)

      if(type == "class") {
        perm_test_y = data.frame(Response = reduced_fit[samp_test,])
        rownames(perm_test_y) <- 1:nrow(perm_test_y)

        probs_random = multivariate_confusion(perm_test_y, Response)[,11]

      } else {
        residuals_test <- data.frame(Response - reduced_fit)[samp_test,]
        if(is.null(ncol(residuals_test))) {
          residuals_test = data.frame(Response = residuals_test)
        }
        rownames(residuals_test) <- 1:nrow(residuals_test)
        perm_test_y <- reduced_fit + residuals_test
        perm_test_y = as.matrix.data.frame(perm_test_y)
      }
      pred_test = matrix(NA, nrow = nrow(perm_test_y), ncol = ncol(perm_test_y))
      Output = numeric(0)

      perm_gcv = matrix(NA, nrow = 1, ncol = ncol(Response))
      place_coef = matrix(NA, nrow = 1, ncol = ncol(Predictor))
      perm_coefs = matrix(NA, nrow = ncol(Response), ncol = ncol(Predictor))
      colnames(perm_coefs) <- colnames(Predictor)

      for(j in 1:ncol(perm_test_y)) {

          pred_perm <- tryCatch({suppressWarnings(caret::train(y = perm_test_y[,j], x = data.frame(Predictor),
                                                              method = caret_parameters$method, trControl = trainControl(method = "none")))
          }, error = function(e) {
            suppressWarnings(caret::train(y = perm_test_y[,j], x = data.frame(Predictor),
                                          method = caret_parameters$method, trControl = trainControl(method = "repeatedcv", number = 5, repeats = 5)))
          }
          )
          pred_test[,j] <- predict.train(pred_perm)

          if(p_calc == "gcv") {

            perm_gcv[1,j] <- if(startsWith(pred_perm$method, "svm")) {
              tryCatch({pred_perm$finalModel@error},
                       error = function(e) {
                         n = nrow(Response[,j])
                         edf = pred_perm$finalModel@nSV
                         rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
                         rss / (n * (1 - edf/n)^2)
                       })
            } else if("fit" %in% names(pred_perm$finalModel)) {
              mean(unlist(lapply(pred_perm$finalModel$fit, function(x) x$gcv)))
            } else if("gcv" %in% names(pred_perm$finalModel)){
              pred_perm$finalModel$gcv
            } else if(length(coef(pred_perm$finalModel)) > 0) {
              n = nrow(Response[,j])
              edf = length(colMeans(abs(coef(pred_perm$finalModel))) >= median(abs(coef(pred_perm$finalModel)), na.rm = TRUE))
              rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
              rss / (n * (1 - edf/n)^2)
            } else if(length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall))) > 0){
              edf = length(caret::varImp(pred_perm$finalModel) %>% filter(Overall >= mean(Overall)))
              n = nrow(Response[,j])
              rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
              rss / (n * (1 - edf/n)^2)
            } else if("mse" %in% names(pred_perm$finalModel)) {
              mean(pred_perm$finalModel$mse)
            } else {
              n = nrow(Response[,j])
              edf = if(!is.null(VIPS)) {length(VIPS)} else {
                lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
                coef_lasso <- coef(lasso, s = lasso$lambda.min)
                length(which(coef_lasso[rownames(coef_lasso) != "(Intercept)"] != 0))
              }
              rss = if(type == "class") (caret::confusionMatrix(predict(pred_perm),sapply(Response[,j],as.factor))$overall[1])^2 else sum(Response[,j] - predict(pred_perm))^2
              rss / (n * (1 - edf/n)^2)
            }
          }
          if(p_calc == "coef") {
            if("svm" %in% pred_perm$method) {
              place_coef = tryCatch({pred_perm$finalModel@coef},
                                    error = function(e){ stats::coef(pred_perm$finalModel)})

              if(is.list(place_coef)){
                names(place_coef) = unique(Response[,j])
                for(i in seq_along(place_coef)){
                  place_coef[[i]] <- sum(-1 >= place_coef[[i]] | 1 <= place_coef[[i]])
                  perm_coefs[j,] <- as.matrix.data.frame(place_coef) / pred_perm$finalModel@nSV
                }
              }
              coef_id = TRUE

            } else if("K" %in% colnames(pred_perm$results)) {
              named <- colnames(Predictor)

              place_coef <- stats::coef(pred_perm$finalModel)

              perm_coefs[j,] <- place_coef

              coef_id = TRUE
            } else if ("lambda" %in% colnames(pred_perm$bestTune) & !isTRUE(coef_id)) {
              place_coef <- stats::coef(pred_perm$finalModel, s = pred_perm$bestTune$lambda)[-1]
              place_coef[is.na(place_coef) | place_coef == "."] <- 0
              place_coef <- get_coef(place_coef, colnames(Predictor))
              coef_id = TRUE

            } else if("fit" %in% names(pred_perm$finalModel)){
              named <- names(stats::coef(pred_perm$finalModel))[-1]
              coefficients_list = list()
              combined_array <- array(NA, dim = c(length(pred_perm$finalModel$fit), length(colnames(Predictor))))
              colnames(combined_array) <- colnames(Predictor)
              for (i in 1:length(pred_perm$finalModel$fit)) {
                resample_name <- paste0("Resample", sprintf("%02d", i))
                coeff <- t(pred_perm$finalModel$fit[[resample_name]]$coefficients)
                if(length(coeff) > 2) {
                  coeff = coeff[1,2:ncol(coeff)]
                } else {
                  coeff = coeff[1,]
                }
                coefficients_list[[i]] <- get_coef(coeff, colnames(Predictor))
                combined_array[i, ] <- coefficients_list[[i]]
              }
              place_coef <- apply(combined_array, 2, mean, na.rm = TRUE)
              perm_coefs[j,] <- place_coef
              colnames(perm_coefs) <- colnames(Predictor)
              coef_id = TRUE
            } else if("coefficients" %in% names(pred_perm$finalModel) | "coef" %in% names(pred_perm$finalModel)) {
              place_coef = t(pred_perm$finalModel$coefficients)
              if(length(place_coef) > 2) {
                place_coef = place_coef[1,2:ncol(place_coef)]
              } else {
                place_coef = place_coef[1,]
              }
              place_coef[is.na(place_coef) | place_coef == "."] <- 0
              place_coef <- get_coef(place_coef, colnames(Predictor))
              perm_coefs[j,] <- place_coef
              colnames(perm_coefs) <- colnames(Predictor)
              coef_id = TRUE
            } else if(length(stats::coef(pred_perm$finalModel)) > 0) {
              place_coef <- stats::coef(pred_perm$finalModel)[-1]
              place_coef[is.na(place_coef) | place_coef == "."] <- 0
              place_coef <- get_coef(place_coef, colnames(Predictor))
              perm_coefs[j,] <- place_coef
              colnames(perm_coefs) <- colnames(Predictor)
              coef_id = TRUE
            } else if(length(caret::varImp(pred_perm$finalModel)) > 0){
              place_coef = t(caret::varImp(pred_perm$finalModel))
              colnames(perm_coefs) <- colnames(Predictor)
              perm_coefs[j,colnames(place_coef)] <- place_coef
              coef_id = TRUE
            } else {
              lasso <- glmnet::cv.glmnet(x = as.matrix.data.frame(Predictor), y = Response[,j], alpha = 0.5, family = family)
              coef_lasso <- coef(lasso, s = lasso$lambda.min)
              perm_coefs[j,] <- place_coef
              colnames(perm_coefs) <- colnames(Predictor)
              coef_id = TRUE
            }
          }
      } # end of j for loop

      if(type == "class") {
        probs_perm = multivariate_confusion(pred_test, perm_test_y)[,11]
      }

      if(p_calc == "coef") {
        Output <- colMeans(abs(coefs))
        rm(perm_coefs, place_coef)
      }
      if(p_calc == "gcv") {
        Output <- mean(perm_gcv, na.rm = TRUE)
      }
      if(p_calc == "F") {

        if(type == "class") {

          SSR_perm <- sum((probs_perm - mean(probs_random))^2)
          SSE_perm <- sum((probs_random - probs_perm)^2)
          MSR <- as.numeric(SSR_perm / p_test)
          MSE <- as.numeric(SSE_perm / (n_test - p_test))
          F_test <- MSR/MSE
        } else {
          SSR_perm <- sum((pred_test - colMeans(perm_test_y))^2)
          SSE_perm <- sum((perm_test_y - pred_test)^2)
          MSR <- as.numeric(SSR_perm / p_test)
          MSE <- as.numeric(SSE_perm / (n_test - p_test))
          F_test <- MSR/MSE
        }
        Output <- F_test
      }
      if(p_calc == "RMSE") {
        if(type == "class") {
          RMSE <- sqrt(mean((probs_random - probs_perm)^2))
        } else {
          RMSE <- sqrt(mean((perm_test_y - pred_test)^2))
        }
        Output <- RMSE
      }

      if(p_calc == "R2") {
        if(type == "class") {
          R2 <- m_R2(as.matrix(probs_random), as.matrix(probs_perm))
        } else {
          R2 <- m_R2(perm_test_y, pred_test)
        }
        Output <- R2
      }

      pred_test = NULL

      return(Output)
    }

    if(isTRUE(parallel)) {
      closeAllConnections()
    }

    if(p_calc == "coef") {
      perm_coefs <- cbind(perm_coefs, unlist(perm_iter))
    } else {
      Output <- c(Output, perm_iter)
    }
    rm(perm_iter)
  } # if convergence was not reached

  if(p_calc == "coef") {
    result_matrix <- matrix(0, nrow = nrow(coefs), ncol = ncol(coefs))
    list_p_values <- vector("list", nrow(coefs))
    p_value = numeric(0)

    array_coefs <- simplify2array(perm_coefs)
    mean_coefs <- apply(array_coefs, c(1, 2), mean)
    sd_coefs <- apply(array_coefs, c(1, 2), sd)
    relevant_predictors = unique(which(abs(coefs) >= median(abs(coefs), na.rm = TRUE), arr.ind = TRUE)[,2])


    Z <- (abs(coefs[,relevant_predictors]) - abs(mean_coefs[,relevant_predictors])) / abs(sd_coefs[,relevant_predictors])

    p_value <- suppressWarnings(mean(p.adjust(p = 1.65 > abs(Z), method = "BH"), na.rm = TRUE))

  }

  if(p_calc == "gcv") {
    Original = min(gcv, na.rm = TRUE)

    Z = (Original - mean(Output)) / sd(Output)

    p_value <- 1 - pnorm(abs(Z))

  }
  if(p_calc == "RMSE") {
    Original = RMSE_test

    Z = (Original - mean(Output)) / sd(Output)

    p_value <- 1 - pnorm(abs(Z))

    #alt p_value <- suppressWarnings(mean(p.adjust(p = Output <= Original, method = "BH"), na.rm = TRUE))
  }
  if(p_calc == "R2") {
    Original = R2_test
    Z = (Original - mean(Output)) / sd(Output)

    p_value <- 1 - pnorm(abs(Z))

    #alt p_value <- suppressWarnings(mean(p.adjust(p = Output >= Original, method = "BH"), na.rm = TRUE))
  }
  if(p_calc == "F") {
    Original = F_test
    Z = (Original - mean(Output)) / sd(Output)

    p_value <- 1 - pnorm(abs(Z))

    # alt p_value <- suppressWarnings(mean(p.adjust(p = Output >= Original, method = "BH"), na.rm = TRUE))
  }

  #____________________________________________________________________________________________________________________
  if (isTRUE(print_progress)) {
    pb_start$tick(10, tokens = list(what = paste("Summarizing Regression Results...")))
  }

    if(type == "class") {

      summary = data.frame(
        n_perm = perm,
        n_boot = n_boot,
        df = ncol(Predictor) - 1,
        summary,
        p_value = p_value
      )

    } else {

      summary = data.frame(
        n_perm = perm,
        n_boot = n_boot,
        df = ncol(Predictor) - 1,
        F_stat = F_test,
        RMSE = RMSE_test,
        R2 = R2_test,
        BIC_test = BIC_test,
        p_value = p_value
      )
    }


    if(!startsWith("pca", dt_parameters$Res_transform) || !startsWith("bgpca", dt_parameters$Res_transform)) {
      Prediction = pred_train_o
    } else {

      Res = as.data.frame(Res_pca_scores)

      if(is.null(Res_ncomp)) {

        names(Res) = paste(rep("Comp"), rep(1:ncol(Response)), sep = "")

      } else {

        names(Res) = paste(rep("Comp"), rep(1:Res_ncomp), sep = "")
        scores = data.frame(Res_pca_scores[,-c(1:ncol(Res))])
        pca = cbind(Res, scores)
      }

      pca = as.matrix.data.frame(Res)
      reconstruction <- pca %*% t(Res_rot)
      Prediction <- sweep(reconstruction, 2, Res_pca_center, "+")

      if(isTRUE(res_3d)) {
        Prediction <- arrayspecs(Prediction, Res_p, Res_dim)
      }
    }

    output <- list(Model = paste(paste0(rep(colnames(data.frame(Response))), collapse = "+"), "~",
                                             paste0(rep(colnames(data.frame(Predictor))), collapse = "+")),
                   Prediction = Prediction,
                   Classification = if(type == "class" && !is.null(which_classify)) pred_special else NULL,
                   summary = summary
                   )
    rownames(output$summary) <- "Regression Statistics"

  class(output) <- "reg.bio"
  #____________________________________________________________________________________________________________________
  if (isTRUE(print_progress)) {
    pb_start$tick(10, tokens = list(what = "Finished!! :D"))
    Sys.sleep(1)
  }

  class(output) <- "reg.bio"
  return(output)
}
