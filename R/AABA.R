#' @title AABA - Analysis of Associations Between Attributes
#'
#' @description
#' AABA is an umbrella function used to analyze relationships between multivariate
#' and multidimensional data. This function takes a model-based hypothesis
#' approach to conduct multiple tests to evaluate one or more hypothesis models
#' related to integration, correlation, covariation, classification, and/or prediction.
#' Data can be subsetted and transformed using the \code{DT} function (see details).
#' AABA has the ability to conduct three subsequent series of analyses: (1)
#' correlation via \code{cor.bio}, (2) covariation/integration via \code{pls.bio},
#' and (3) multivariate prediction or classification via \code{pred}.
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
#' of response and predictor subsets
#' @param cors Logical value indicating whether to conduct correlation analysis
#' @param vips Logical value indicating whether to conduct PLS and variable importance analysis
#' @param covs Logical value indicating whether to conduct a shape specific covariation analysis
#' @param regs Logical value indicating whether to conduct prediction/classification analysis
#' @param PSLR Logical value indicating whether to conduct Procrustes shape linear regression.
#' @param VIP_trim Logical value indicating whether to trim regression models by VIP scores calculated in vips. Only applicable when vips = TRUE as well as either regs or PSLR is set to TRUE.
#' @param key Optional character string for naming output directories
#' @param parallel Logical value indicating whether to use parallel processing
#' @param core_choice Character value specifying number of cores to use
#' @param print_progress Logical value indicating whether to display progress bars
#' @param ... Additional arguments passed to underlying functions. Can include:
#'   \itemize{
#'     \item{cor.bio parameters: method}
#'     \item{pls.bio parameters: Res_transform, Pred_transform, Res_ncomp, Pred_ncomp}
#'     \item{reg.bio parameters: package, p_calc, n_boot, perm, CV, which_classify}
#'     \item{pslr.bio parameters: write.results, restructured.results}
#'   }
#'
#' @returns A list containing results from all enabled analyses:
#' \itemize{
#'   \item{Cor.Results}{If cors=TRUE, correlation analysis results:
#'     \itemize{
#'       \item{correlation results}{Correlation matrices}
#'       \item{Mantel/Canonical Correlation statistics}{Matrix correlation tests}
#'       \item{Subset analyses}{If point_set provided}
#'       \item{Correlograms}{In folder directory "AABA Results"}
#'     }
#'   }
#'   \item{VIP.Results}{If vips=TRUE, PLS and variable importance results:
#'     \itemize{
#'       \item{VIP scores}{Variable importance in projection scores}
#'       \item{PLS model statistics}{degrees of freedom, R2, Q2, RMSEP, p-value}
#'       \item{Visualization plots}{VIP score plots}
#'     }
#'   }
#'   \item{Reg.Results}{If regs=TRUE, prediction/classification results:
#'     \itemize{
#'       \item{Model performance}{Accuracy/error metrics}
#'       \item{Predictions}{Predicted values/classes}
#'       \item{Statistical significance}{p-values, regression statistics}
#'     }
#'   }
#'   \item{PSLR.Results}{If PSLR=TRUE and Response contains landmarks:
#'     \itemize{
#'       \item{Shape regression results}{Procrustes regression statistics}
#'       \item{Shape changes}{Predicted shape deformations via \code{lolshape}}
#'     }
#'   }
#' }
#'
#' @details
#'
#' AABA integrates multiple statistical approaches to provide a comprehensive analysis
#' of relationships between multivariate datasets. The function processes each hypothesis
#' model sequentially through up to four analytical stages, depending on the enabled options.
#' Multiple tests of subsetted data for a single test within a hypothesis model
#' is also supported through the point_set parameter, enabling investigation of
#' specific variable combinations or anatomical regions. When paired_subsets=TRUE,
#' all possible combinations of response and predictor subsets are analyzed,
#' providing a comprehensive view of relationships between data sets.
#'
#' When cors = TRUE, several correlation and integration statistics will be estimated
#' to evaluate the associations between multivariate data sets via via \code{cor.bio}.
#' If the multivariate data sets structures are identical, a permuted Mantel test
#' is conducted, otherwise canonical correlation analysis is performed. Correlograms
#' of correlation matrices are generated as an output to visualize variable associations.
#' See \code{\link{cor.bio}} for more details.
#'
#' When vips = TRUE, a user specified partial least squares (PLS) and variable importance
#' in projection (VIP) analysis is performed for each test within each hypothesis model,
#' This part of the function also generates summary statistics from the PLS analyses
#' as well as VIP plots to assess variable importance when vip_plots = TRUE.
#' See \code{\link{pls.bio}} for more details.
#'
#' When covs = TRUE, a shape specific two-block PLS (2B-PLS) covariation analysis
#' is conducted and subsequently will generate interactive .html files in a user
#' specified directory to view the major shape changes associated with each latent
#' variable using the \code{lolshape} function. See \code{\link{pls.bio}} for more details.
#'
#' When regs = TRUE, the function will conduct a permuted, cross-validated, and/or
#' bootstrapped prediction or classification regression via the \code{reg.bio} function.
#' This function, like the others, permits for both data subsetting as well offers
#' different methods to test for statistical significance using a flexible
#' Randomized Residual Permutation Procedure (RRPP) which has been adapted for
#' multivariate data. See \code{\link{reg.bio}} for more details.
#'
#' When PSLR = TRUE, a Procrustes Shape Linear Regression (PSLR) is performed via
#' the \code{pslr.bio} function. This function specifically handles landmark data as
#' response (dependent) data, conducting shape regression analysis while
#' preserving geometric relationships. See \code{\link{pslr.bio}} for more details.
#'
#' Data transformation can be parameterized through dt_parameters within the Models
#' structure as the third object within a test within a hypothesis model or
#' through arguments passed via the argument (...). These parameters modify how
#' data is processed before each analytical stage, allowing flexibility in
#' handling different data types.
#'
#' @seealso \code{\link{cor.bio}} \code{\link{pls.bio}} \code{\link{reg.bio}} \code{\link{pslr.bio}}
#'
#' @author Keeling et al., 2025
#'
#' @references
#' Placeholder
#'
#' @importFrom stats prcomp cor sd median p.adjust pnorm
#' @importFrom dplyr filter mutate arrange group_by summarise intersect
#' @importFrom tidyr pivot_longer
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom progress progress_bar
#' @importFrom utils modifyList
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Analysis with geometric morphometric data
#' # Test relationships between mandibular shape and biomechanical properties
#'
#' # Create model structure
#' Models <- list(
#'   'Mandible_Analysis' = list(
#'     'Symphysis_Biomechanics' = list(
#'       'Response' = Corpus_Land[241:360,,],  # Symphysis landmarks
#'       'Predictor' = Biomechanics %>%
#'         filter(Region == "LM1M2") %>%
#'         select(5:93)  # Biomechanical variables
#'     )
#'   )
#' )
#'
#' # Define variable subsets
#' point_set <- list(
#'   'Mandible_Analysis' = list(
#'     'Symphysis_Biomechanics' = list(
#'       'Response' = list('Anterior' = 1:60, 'Posterior' = 61:120),
#'       'Predictor' = list(
#'         'Dimensions' = c(5:8, 10:15),
#'         'Cortical_Thickness' = c(4, 16:18, 30:60)
#'       )
#'     )
#'   )
#' )
#'
#' # Run complete analysis including shape regression
#' results_shape <- AABA(
#'   Models = Models,
#'   point_set = point_set,
#'   cors = TRUE,
#'   vips = TRUE,
#'   regs = TRUE,
#'   PSLR = TRUE,  # Include shape regression
#'   vip_method = "spls",
#'   parallel = TRUE
#' )
#'
#' # Example 2: Analysis with standard multivariate data
#' # Test relationships between biomechanical properties at different cross-sections
#'
#' # Create model structure for non-shape analysis
#' Models_standard <- list(
#'   'Biomechanics_Analysis' = list(
#'     'CrossSection_Comparison' = list(
#'       'Response' = Biomechanics %>%
#'         filter(Region == "LM1M2") %>%
#'         select(5:93),
#'       'Predictor' = Biomechanics %>%
#'         filter(Region == "LP3P4") %>%
#'         select(5:93)
#'     )
#'   )
#' )
#'
#' # Run analysis without shape regression
#' results_standard <- AABA(
#'   Models = Models_standard,
#'   cors = TRUE,
#'   vips = TRUE,
#'   regs = TRUE,
#'   PSLR = FALSE,  # Skip shape regression
#'   vip_method = "spls",
#'   parallel = TRUE
#' )
#' }

AABA <- function(Models,
                 point_set = NULL,
                 paired_subsets = FALSE,
                 cors = TRUE,
                 vips = TRUE,
                 covs = TRUE,
                 regs = TRUE,
                 PSLR = FALSE,
                 VIP_trim = TRUE,
                 key = NULL,
                 parallel = TRUE,
                 core_choice = "detect",
                 print_progress = TRUE,
                 ...) {

  AABA_Output = list()

  if (isTRUE(print_progress)) {
    pb <- progress::progress_bar$new(
      format = "Processing :what [:bar] :percent ETA: :eta",
      clear = TRUE, total = 100, width = 80)

    pb$tick(5, tokens = list(what = "Starting AABA Analysis"))
  }

  Directory = getwd()

  if(!is.null(key)) {
    if(!dir.exists(paste("AABA Results", "-", key))) {
      dir.create(paste("AABA Results", "-", key))
    }
    ct.dir <- paste("AABA Results", "-", key)
  } else {
    if(!dir.exists("AABA Results")) {
      dir.create("AABA Results")
    }
    ct.dir <- "AABA Results"
  }
  path = paste(Directory, "/", ct.dir, sep = "")
  setwd(path)

  key = NULL
  args = list(...)
  subset = if(!is.null(point_set)) TRUE else FALSE

  dt_params = names(formals(DT))

  dt_parameters = list(dt_method = "procdist")

  dt_parameters = modifyList(dt_parameters,args) %>% .[intersect(names(.),dt_params)]

  #______________________________________________________________________________________________________________________


  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = "Computing: Correlation Analysis"))
  }

  if(isTRUE(cors)) {

    cor_parameters = list(Models = Models,
                          point_set = point_set,
                          paired_subsets = paired_subsets,
                          parallel = parallel,
                          core_choice = core_choice,
                          key = NULL) %>%
                          modifyList(., dt_parameters)

    Cor.Results <- do.call(CSGM::cor.bio, cor_parameters)
    AABA_Output$Cor.Results <- Cor.Results
    setwd(path)
    gc()
    closeAllConnections()
  }

  if (isTRUE(print_progress)) {
    pb$tick(20, tokens = list(what = "Finished: Correlation Analysis"))
  }
  #___________________________________________________________________________________________________________________________

  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = "Computing: VIP & Covariation Analysis"))
  }

  if(isTRUE(vips)) {

    vip_parameters = modifyList(args, list(Models = Models,
                                           point_set = point_set,
                                           subset = subset,
                                           paired_subsets = paired_subsets,
                                           key = NULL))

      VIP.Results <- do.call(CSGM::pls.bio, vip_parameters)
      AABA_Output$VIP.Results <- VIP.Results
      setwd(path)
      gc()
      closeAllConnections()
  }

  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = "Finished: VIP & Covariation Analysis"))
  }
  #___________________________________________________________________________________________________________________________
  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = "Computing: Procrustes Shape Regression"))
  }

  if(isTRUE(PSLR)) {

    PSLR.Results = list()

    if(isTRUE(VIP_trim)) { # decision to trim by vips or not

      Vars = vector("list", length(Models))
      Vars.subset = vector("list", length(Models))
      for(m in seq_along(AABA_Output$VIP.Results)) {
        Vars = vector("list", length(Models[[m]]))
        Vars.subset = vector("list", length(Models[[m]]))
        for(i in seq_along(AABA_Output$VIP.Results[[m]]$plsr.results.all)) {
          if("vip" %in% names(AABA_Output$VIP.Results[[m]]$plsr.results.all[[i]])) {
            Vars[[m]][[i]] = colnames(data.frame(AABA_Output$VIP.Results[[m]]$plsr.results.all[[i]]$vip))
          }
          if(!is.null(point_set)){
            Vars.subset[[m]][[i]] = list()
            for(j in seq_along(AABA_Output$VIP.Results[[m]]$plsr.results[[i]])){
              if("vip" %in% names(AABA_Output$VIP.Results[[m]]$plsr.results[[i]])) {
                Vars.subset[[m]][[i]][[j]] = colnames(data.frame(AABA_Output$VIP.Results[[m]]$plsr.results[[i]]$vip))
              } else {
                Vars.subset[[m]][[i]][[j]] = colnames(data.frame(AABA_Output$VIP.Results[[m]]$plsr.results[[i]][[j]]$vip))
              }
            } # j loop
          } # for subsets
        } # i loop
      } # m loop
    } # for VIPs trimming


    pslr_parameters = list(Models = Models,
                           paired_subsets = paired_subsets,
                           point_set = point_set,
                           all_vips = if(isTRUE(VIP_trim)) Vars else NULL,
                           subset_vips = if(isTRUE(VIP_trim)) Vars.subset else NULL,
                           key = NULL) %>% modifyList(., dt_parameters)

    PSLR.Results <- do.call(CSGM::pslr.bio, pslr_parameters)

    AABA_Output$PSLR.Results <- PSLR.Results
    gc()
    closeAllConnections()
  }

  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = "Finished: Proc Shape Regression"))
  }

  #___________________________________________________________________________________________________________________________

  if (isTRUE(print_progress)) {
    pb$tick(20, tokens = list(what = "Computing: Multivariate Regressions"))
  }

  if(isTRUE(regs)) {


    caret_params = c("method","preProcess","weights","metric","maximize",
                     "trControl","tuneGrid","tuneLength","mtry","replace","sampsize",
                     "nodesize","maxnodes","importance","localImp","nPerm","proximity",
                     "oob.prox","norm.votes","keep.forest","keep.inbag")

    dt_params = names(formals(DT))

    pred_params = names(formals(reg.bio))

    reg_params = c(caret_params, dt_params, pred_params)

    reg_parameters = args %>% .[intersect(names(.),reg_params)]

    if(isTRUE(VIP_trim)) {
      reg_parameters <- modifyList(reg_parameters, list(full_model = FALSE))

      Vars = vector("list", length(Models))
      Vars.subset = vector("list", length(Models))
      for(m in seq_along(AABA_Output$VIP.Results)) {
        Vars = vector("list", length(Models[[m]]))
        Vars.subset = vector("list", length(Models[[m]]))
        for(i in seq_along(AABA_Output$VIP.Results[[m]]$plsr.results.all)) {
          if("vip" %in% names(AABA_Output$VIP.Results[[m]]$plsr.results.all[[i]])) {
            Vars[[m]][[i]] = colnames(data.frame(AABA_Output$VIP.Results[[m]]$plsr.results.all[[i]]$vip))
          }
          if(!is.null(point_set)){
            Vars.subset[[m]][[i]] = list()
            for(j in seq_along(AABA_Output$VIP.Results[[m]]$plsr.results[[i]])){
              if("vip" %in% names(AABA_Output$VIP.Results[[m]]$plsr.results[[i]])) {
                Vars.subset[[m]][[i]][[j]] = colnames(data.frame(AABA_Output$VIP.Results[[m]]$plsr.results[[i]]$vip))
              } else {
                Vars.subset[[m]][[i]][[j]] = colnames(data.frame(AABA_Output$VIP.Results[[m]]$plsr.results[[i]][[j]]$vip))
              }
            } # j loop
          } # for subsets
        } # i loop
      } # m loop
    } # for VIPs trimming

    Reg.Results = vector("list", length(Models))

    for(m in seq_along(Reg.Results)) {
      Reg.Results[[m]]$pred.all <- vector("list", length(Models[[m]]))

      if(isTRUE(subset)) {
        Reg.Results[[m]]$pred.subset <- vector("list", length(Models[[m]]))
      }
      for(i in seq_along(Reg.Results[[m]]$pred.all)) {

        if(isTRUE(subset)){
          if(isTRUE(paired_subsets)) {
            p_subset = expand.grid(seq_along(point_set[[m]][[i]][[1]]),
                                   seq_along(point_set[[m]][[i]][[2]]))

            Reg.Results[[m]]$pred.subset <- vector("list", nrow(p_subset))
          } else {
            p_subset = data.frame(Response = rep(1:length(point_set[[m]][[i]][[1]])),
                                  Predictor = rep(1:length(point_set[[m]][[i]][[2]])))

            Reg.Results[[m]]$pred.subset[[i]] <- vector("list", length(point_set[[m]][[i]][[2]]))
          }

          for(j in 1:nrow(p_subset)) {

            dt_parameters2 <- tryCatch({modifyList(dt_parameters,Models[[m]][[i]][[3]])},error = function(e){dt_parameters})
            dt_parameters2 = modifyList(dt_parameters2, list(Response = Models[[m]][[i]][[1]],
                                                           Predictor = Models[[m]][[i]][[2]],
                                                           Res_point_set = point_set[[m]][[i]][[1]][[p_subset[j,1]]],
                                                           Pred_point_set = point_set[[m]][[i]][[2]][[p_subset[j,2]]]
                                                           )) %>% .[intersect(names(.), dt_params)]


            dt_parameters2 = modifyList(dt_parameters2, list(VIPS = if(isTRUE(VIP_trim)) Vars.subset[[m]][[i]][[j]] else NULL))

            reg_parameters2 = modifyList(reg_parameters, dt_parameters2) %>% .[intersect(names(.),reg_params)]

            Reg.Results[[m]]$pred.subset[[i]][[j]] <- do.call(CSGM::reg.bio, reg_parameters2)
            closeAllConnections()
          } # j loop
        } # Subset

        dt_parameters1 <- tryCatch({modifyList(dt_parameters, Models[[m]][[i]][[3]])},error = function(e){dt_parameters})
        dt_parameters1 = modifyList(dt_parameters1, list(Response = Models[[m]][[i]][[1]],
                                                       Predictor = Models[[m]][[i]][[2]],
                                                       Res_point_set = NULL,
                                                       Pred_point_set = NULL)) %>% .[intersect(names(.), dt_params)]


        dt_parameters1 = modifyList(dt_parameters1, list(VIPS = if(isTRUE(VIP_trim)) Vars[[m]][[i]] else NULL))

        reg_parameters1 = modifyList(reg_parameters, dt_parameters1) %>% .[intersect(names(.),reg_params)]

        Reg.Results[[m]]$pred.all[[i]] <- do.call(CSGM::reg.bio, reg_parameters1)
        closeAllConnections()
      } # i loop
    } # m loop
    names(Reg.Results) <- names(Models)
    AABA_Output$Reg.Results <- Reg.Results
    gc()
    closeAllConnections()
  }

  # summarize results___________________________________________________________________________________________________________________________

  if (isTRUE(print_progress)) {
    pb$tick(25, tokens = list(what = "Finished: Adv. Regression Analysis"))
  }
  setwd(Directory)

  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = "Finished!: AABA_Output Analysis :D"))
  }
  class(AABA_Output) <- "AABA"

  return(AABA_Output)
} # End of AABA_Output function
