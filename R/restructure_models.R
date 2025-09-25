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
