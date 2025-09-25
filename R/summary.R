#' @title Summarize Results
#'
#' @description A function used to extract and summarize the results from
#' regression based functions within the CSGM package.
#'
#' @param object An object resulting from either the \code{reg.bio}, \code{PSLR},
#' and \code{AABA} functions. Otherwise, the function uses the default \pkg{stats} package.
#'
#' @returns A data table consisting of statistics for model performance.
#'
#' @details This function extracts the regression results for the \code{GMA},
#' \code{PSLR} and \code{reg.bio} functions as well as conducts an internal function
#' \code{AABA.summary} specifically designed to extract relevant statistics for
#' the \code{AABA} function.
#'
#' @author Keeling et al., 2025
#'
#' @seealso \code{\link{AABA}} \code{\link{GMA}} \code{\link{reg.bio}}
#'
#' @importFrom stats as.formula
#' @importFrom utils head tail
#' @importFrom writexl write_xlsx
#' @importFrom caret trainControl
#' @importFrom doParallel registerDoParallel
#' @importFrom pls plsr
#' @importFrom purrr map_df
#'
#' @export

summary <- function(object,
                    ...) {

  if("reg.bio" %in% class(object)) {
    print("Summary Statistics from the reg.bio Function:")
    cat("\n")
    cat("\n")
    return(object$summary)

  } else if("GMA" %in% class(object)) {

    print("Summary of the Geometric Morphometric Analysis:")
    cat("\n")
    cat("\n")
    if(grepl("sym_data", names(object))) {
      print("Summary of Bilateral Symmetry")
      cat("\n")
      print(object$sym_data)
      cat("\n")
      cat("\n")
    }
    print("PCA Summary")
    print(head(object$PCA$PCA_var$PCA_var))
    cat("\n")
    cat("\n")
    print("Plot 1: Morphometric Relationships")
    plot(object$morphotree$resid$dendro)
    print("Plot 2: Principal Component Plot for PCs 1-3")
    object$PCA$PCA_plots$PCs_1_2_3
    print("Plot 3: Shape variation across PC 1")
    object$PCA$lollipop_PCs[[1]]
    print("Plot 4: Shape variation across PC 2")
    object$PCA$lollipop_PCs[[2]]
    cat("\n")

  } else if("AABA" %in% class(object)) {

    print("Summary Statistics from AABA:")

    if ("Cor.Results" %in% names(object)) {
      object_cor <- object$Cor.Results
      m_length = length(object_cor)
      Cor_Output = vector("list", length(object_cor))
      names(Cor_Output) = paste(rep("Hypothesis Model", length(object_cor)), rep(1:length(object_cor)),sep = "_")
      for(m in seq_along(object_cor)) {
        i_length = length(object_cor[[m]])
        Cor_Output[[m]] = vector("list", length(object_cor[[m]]$Cor.all.Mantels))
        names(Cor_Output[[m]]) = (paste(rep("Correlation Tests", length(object_cor[[m]]$Cor.all.Mantels)), rep(1:length(object_cor[[m]]$Cor.all.Mantels))))
        for(i in seq_along(object_cor[[m]]$Cor.all.Mantels)) {
          Cor.all.Mantels = data.frame(object_cor[[m]]$Cor.all.Mantels[[i]])
          if(ncol(Cor.all.Mantels) < 4) {
            Cor.all.Mantels = t(Cor.all.Mantels)
          }
          Cor_Output[[m]][[i]] = Cor.all.Mantels
          colnames(Cor_Output[[m]][[i]]) = c("r","p-value","I","Ix","Iy")
          rownames(Cor_Output[[m]][[i]])[1] = paste("Overall_Test_", i, sep = "")

          if("Cor.subset.Mantels" %in% names(object_cor[[m]])) {
            for(j in seq_along(object_cor[[m]]$Cor.subset.Mantels[[i]])) {
              Cor.subset.Mantels = data.frame(object_cor[[m]]$Cor.subset.Mantels[[i]][[j]])
              if(ncol(Cor.subset.Mantels) < 4) {
                Cor.subset.Mantels = t(Cor.subset.Mantels)
              }
              colnames(Cor.subset.Mantels) = c("r","p-value","I","Ix","Iy")
              Cor_Output[[m]][[i]] <- rbind(Cor_Output[[m]][[i]],Cor.subset.Mantels)
              rownames(Cor_Output[[m]][[i]])[j+1] = paste("Test_", i, "_Subset_", j, sep = "")
            }
          }
        }
      }
    } else {
      Cor_Output = NULL
    }

    if ("VIP.Results" %in% names(object)) {
      object_vip <- object$VIP.Results
      m_length = length(object_vip)
      VIP_Output = vector("list", length(object_vip))
      names(VIP_Output) = paste(rep("Hypothesis Model", length(object_vip)), rep(1:length(object_vip)),sep = "_")
      for(m in 1:m_length) {
        i_length = length(object_vip[[m]])
        VIP_Output[[m]] = vector("list", length(object_vip[[m]]$plsr.results.all))
        names(VIP_Output[[m]]) = paste(rep("Covariation Tests", length(object_vip[[m]]$plsr.results)), rep(1:length(object_vip[[m]]$plsr.results)), sep = "_")
        for(i in 1:i_length) {
          VIP_Output[[m]][[i]] = data.frame(object_vip[[m]]$plsr.results.all[[i]]$summary)
          colnames(VIP_Output[[m]][[i]]) <- c("X_var", "R2", "Q2", "RMSEP", "p_RMSEP", "p_Q2")
          rownames(VIP_Output[[m]][[i]])[1] = paste("Overall_Test_", i, sep = "")
          if("Cov.all" %in% names(object_vip[[m]])){
            cov.results <- data.frame(object_vip[[m]]$Cov.all[[i]]$summary, object_vip[[m]]$Cov.all.VIP[[i]]$summary)
            rownames(VIP_Output[[m]][[i]])[1] = paste("Overall_Test_", i, sep = "")
            colnames(VIP_Output[[m]][[i]]) <- c("df", "RV", "p-RV", "VIP_df", "VIP_RV", "VIP_p-RV")
          }
          if("plsr.results" %in% names(object_vip[[m]])){

            j_length <- if("summary" %in% names(object_vip[[m]]$plsr.results[[i]])) 1 else length(object_vip[[m]]$plsr.results[[i]])

            for(j in 1:j_length) {
              VIP_Output[[m]][[i]][j+1,] <- if("Cov.subset" %in% names(object_vip[[m]])) { data.frame(object_vip[[m]]$plsr.results[[i]][[j]]$summary,
                                                                                           object_vip[[m]]$Cov.subset[[i]][[j]]$summary,
                                                                                           object_vip[[m]]$Cov.subset.VIP[[i]][[j]]$summary)
                                             } else {
                                               if("summary" %in% names(object_vip[[m]]$plsr.results[[i]])){
                                                 data.frame(object_vip[[m]]$plsr.results[[i]]$summary)
                                               }else {
                                                 data.frame(object_vip[[m]]$plsr.results[[i]][[j]]$summary)
                                               }
                                             }

              rownames(VIP_Output[[m]][[i]])[j+1] = paste("Test_", i, "_Subset_", j, sep = "")
            } # J loop
          } # if subsets
        } # I loop
      } # M loop
    } else {
      VIP_Output = NULL
    }

    if ("PSLR.Results" %in% names(object)) {
      object_pslr = object$PSLR.Results
      m_length = length(object_pslr)
      PSLR_Output = vector("list", m_length)
      names(PSLR_Output) = paste(rep("Hypothesis Model", m_length), rep(1:m_length),sep = "_")
      for(m in seq_along(object_pslr)) {
        i_length = length(object_pslr[[m]]$pslr.results.all)
        PSLR_Output[[m]] = vector("list", i_length)
        names(PSLR_Output[[m]]) = paste(rep("Procrustes Regression Test", i_length), rep(1:i_length),sep = "_")
        for(i in seq_along(object_pslr[[m]]$pslr.results.all)) {
          PSLR_Output[[m]][[i]] = data.frame(object_pslr[[m]]$pslr.results.all[[i]]$summary)[,1:14]
          rownames(PSLR_Output[[m]][[i]])[[1]] = paste("Overall Test_", i, sep = "")
          colnames(PSLR_Output[[m]][[i]]) <- c("df","SS","MS","Rsq","F","Z","p_value","QR_df","QR_SS","QR_MS","QR_Rsq","QR_F","QR_Z","QR_p_value")
          if("pslr.results" %in% names(object_pslr[[m]])) {
            for(j in seq_along(object_pslr[[m]]$pslr.results[[i]])) {
              if("summary" %in% names(object_pslr[[m]]$pslr.results[[i]])) {
                PSLR_Output[[m]][[i]][j+1,] = data.frame(object_pslr[[m]]$pslr.results[[i]]$summary)[,1:14]
                rownames(PSLR_Output[[m]][[i]])[j+1] = paste("Test_", i, "_Subset_", j, sep = "")
              } else {
                PSLR_Output[[m]][[i]][j+1,] = data.frame(object_pslr[[m]]$pslr.results[[i]][[j]]$summary)[,1:14]
                rownames(PSLR_Output[[m]][[i]])[j+1] = paste("Test_", i, "_Subset_", j, sep = "")
              }
            }
          }
        } # I loop
      } # m loop
    } else {
      PSLR_Output = NULL
    }

    if("Reg.Results" %in% names(object)) {
      object_pred = object$Reg.Results
      m_length = length(object_pred)
      pred_Output = vector("list", length(object_pred))
      names(pred_Output) = paste(rep("Hypothesis Model", m_length), rep(1:m_length),sep = "_")
      for(m in 1:m_length) {
        i_length = length(object_pred[[m]]$pred.all)
        pred_Output[[m]] = vector("list", i_length)
        names(pred_Output[[m]]) = paste(rep("Hypothesis Model", m_length), rep(1:m_length),sep = "_")
        for(i in 1:i_length) {
          pred.all = data.frame(object_pred[[m]]$pred.all[[i]]$summary)
          pred_Output[[m]][[i]] = data.frame(matrix(nrow = 1, ncol = ncol(pred.all)))
          colnames(pred_Output[[m]][[i]]) = colnames(pred.all)
          pred_Output[[m]][[i]][1,] <- pred.all
          rownames(pred_Output[[m]][[i]])[1] = paste("Overall_Test_", i, sep = "")

          if("pred.subset" %in% names(object_pred[[m]])) {
            j_length =  if("summary" %in% object_pred[[m]]$pred.subset[[i]]) 1 else length(object_pred[[m]]$pred.subset[[i]])
            for(j in 1:j_length) {
              if("summary" %in% object_pred[[m]]$pred.subset[[i]]) {
                pred_subset = data.frame(object_pred[[m]]$pred.subset[[i]]$summary)
              } else {
                pred.subset = data.frame(object_pred[[m]]$pred.subset[[i]][[j]]$summary)
              }
              pred_Output[[m]][[i]][j+1,] <- pred.subset
              rownames(pred_Output[[m]][[i]])[j+1] = paste("Test_", i, "_Subset_", j, sep = "")
            } # j loop
          } # subset detection
        } # i loop
      } # m loop
    } else {
      pred_Output = NULL
    }

   Output <- vector("list", m_length)
   names(Output) <- paste(rep("Hypothesis_Model", m_length), 1:m_length, sep="_")

   for(m in 1:m_length) {
     Output[[m]] <- vector("list", i_length)
     names(Output[[m]]) <- paste(rep("Test", i_length), 1:i_length, sep="_")

     for(i in 1:i_length) {

       if(!is.null(Cor_Output) && !is.null(Cor_Output[[m]][[i]])) {
         combined_df <- data.frame(row.names = rownames(Cor_Output[[m]][[i]]))
       } else if(!is.null(VIP_Output) && !is.null(VIP_Output[[m]][[i]])) {
         combined_df <- data.frame(row.names = rownames(VIP_Output[[m]][[i]]))
       } else if(!is.null(PSLR_Output) && !is.null(PSLR_Output[[m]][[i]])) {
         combined_df <- data.frame(row.names = rownames(pslr_Output[[m]][[i]]))
       } else if(!is.null(pred_Output) && !is.null(pred_Output[[m]][[i]])) {
         combined_df <- data.frame(row.names = rownames(pred_Output[[m]][[i]]))
       } else {
         combined_df <- data.frame()
       }

       if(!is.null(Cor_Output) && !is.null(Cor_Output[[m]][[i]])) {

         colnames(Cor_Output[[m]][[i]]) <- paste("Cor", colnames(Cor_Output[[m]][[i]]), sep="_")
         combined_df <- cbind(combined_df, Cor_Output[[m]][[i]])
       }

       if(!is.null(VIP_Output) && !is.null(VIP_Output[[m]][[i]])) {
         colnames(VIP_Output[[m]][[i]]) <- paste("VIP", colnames(VIP_Output[[m]][[i]]), sep="_")
         combined_df <- cbind(combined_df, VIP_Output[[m]][[i]])
       }

       if(!is.null(PSLR_Output) && !is.null(PSLR_Output[[m]][[i]])) {
         colnames(PSLR_Output[[m]][[i]]) <- paste("PSLR", colnames(PSLR_Output[[m]][[i]]), sep="_")
         combined_df <- cbind(combined_df, PSLR_Output[[m]][[i]])
       }

       if(!is.null(pred_Output) && !is.null(pred_Output[[m]][[i]])) {
         colnames(pred_Output[[m]][[i]]) <- paste("Pred", colnames(pred_Output[[m]][[i]]), sep="_")
         combined_df <- cbind(combined_df, pred_Output[[m]][[i]])
       }

       Output[[m]][[i]] <- combined_df
     }
   }
   return(Output)

  } else if ("cor.bio" %in% class(object)) {
    print("Summary Statistics from cor.bio")
    Output = vector("list", length(object))
    names(Output) = paste(rep("Hypothesis Model", length(object)), rep(1:length(object)),sep = "_")
    for(m in seq_along(object)) {
      Output[[m]] = vector("list", length(object[[m]]$Cor.all.Mantels))
      names(Output[[m]]) = (paste(rep("Correlation Tests", length(object[[m]]$Cor.all.Mantels)), rep(1:length(object[[m]]$Cor.all.Mantels))))
      for(i in seq_along(object[[m]]$Cor.all.Mantels)) {
        Output[[m]][[i]] = data.frame(matrix(nrow = 1, ncol = length(object[[m]]$Cor.all.Mantels[[i]])))
        colnames(Output[[m]][[i]]) = c("r","p-value","I","Ix","Iy")
        Output[[m]][[i]][1,] = object[[m]]$Cor.all.Mantels[[i]]
        rownames(Output[[m]][[i]])[1] = paste("Overall_Test_", i, sep = "")

        if("Cor.subset.Mantels" %in% names(object[[m]])) {
          for(j in seq_along(object[[m]]$Cor.subset.Mantels[[i]])) {
            Output[[m]][[i]][j+1,] = object[[m]]$Cor.subset.Mantels[[i]][[j]]
            rownames(Output[[m]][[i]])[j+1] = paste("Test_", i, "_Subset_", j, sep = "")
          }
        }
      }
    }
    return(Output)

  } else if ("Bio.VIP" %in% class(object)) {
    print("Summary Statistics from Bio.VIP")
    Output = vector("list", length(object))
    names(Output) = paste(rep("Hypothesis Model", length(object)), rep(1:length(object)),sep = "_")
    for(m in seq_along(object)) {
      Output[[m]] = vector("list", length(object[[m]]$plsr.results.all))
      names(Output[[m]]) = paste(rep("Covariation Tests", length(object[[m]]$plsr.results)), rep(1:length(object[[m]]$plsr.results)), sep = "_")
      for(i in seq_along(object[[m]]$plsr.results.all)) {
        Output[[m]][[i]] = data.frame(object[[m]]$plsr.results.all[[i]]$summary)
        colnames(Output[[m]][[i]]) = c("X_var", "R2", "Q2", "RMSEP", "p_RMSEP", "p_Q2")
        rownames(Output[[m]][[i]])[1] = paste("Overall_Test_", i, sep = "")
        if("Cov.all" %in% names(object[[m]])){
          Output[[m]][[i]] <- data.frame(Output[[m]][[i]], object[[m]]$Cov.all[[i]]$summary, object[[m]]$Cov.all.VIP[[i]]$summary)
          rownames(Output[[m]][[i]])[1] = paste("Overall_Test_", i, sep = "")
        }
        if("plsr.results" %in% names(object[[m]])){
          for(j in seq_along(object[[m]]$plsr.results[[i]])) {
            if("summary" %in% names(object[[m]]$plsr.results[[i]])){
              plsr_summary = data.frame(object[[m]]$plsr.results[[i]]$summary)
              Cov_summary = data.frame(object[[m]]$Cov.subset[[i]][[1]]$summary)
              cvip_summary = data.frame(object[[m]]$Cov.subset.VIP[[i]][[1]]$summary)

              Output[[m]][[i]][2,] <- if("Cov.subset" %in% names(object[[m]])) {
                                          data.frame(plsr_summary, Cov_summary, cvip_summary)
                                        } else {
                                          data.frame(plsr_summary)
                                        }

            } else {
              Output[[m]][[i]][j+1,] <- if("Cov.subset" %in% names(object[[m]])){
                                          data.frame(data.frame(object[[m]]$plsr.results[[i]][[j]]$summary),
                                          data.frame(object[[m]]$Cov.subset[[i]][[j]]$summary),
                                          data.frame(object[[m]]$Cov.subset.VIP[[i]][[j]]$summary))
                                        } else {
                                          data.frame(object[[m]]$plsr.results[[i]][[j]]$summary)
                                        }
              rownames(Output[[m]][[i]])[j+1] = paste("Test_", i, "_Subset_", j, sep = "")
            } # adjustments
          } # J loop
        } # if subsets
      }  # I loop
    }# M loop
    return(Output)

  } else if ("PSLR" %in% class(object)) {
    print("Summary Statistics from PSLR")
    m_length = length(object)
    Output = vector("list", m_length)
    names(Output) = paste(rep("Hypothesis Model", m_length), rep(1:m_length),sep = "_")
    for(m in seq_along(object)) {
      i_length = length(object[[m]]$pslr.results.all)
      Output[[m]] = vector("list", i_length)
      names(Output[[m]]) = paste(rep("Procrustes Regression Test", i_length), rep(1:i_length),sep = "_")
      for(i in seq_along(object[[m]]$pslr.results.all)) {
        Output[[m]][[i]] = data.frame(object[[m]]$pslr.results.all[[i]]$summary)[,1:14]
        rownames(Output[[m]][[i]])[[1]] = paste("Overall Test_", i, sep = "")
        colnames(Output[[m]][[i]]) <- c("df","SS","MS","Rsq","F","Z","p_value","QR_df","QR_SS","QR_MS","QR_Rsq","QR_F","QR_Z","QR_p_value")
        if("pslr.results" %in% names(object[[m]])) {
          for(j in seq_along(object[[m]]$pslr.results[[i]])) {
            if("summary" %in% names(object[[m]]$pslr.results[[i]])) {
              Output[[m]][[i]][j+1,] = data.frame(object[[m]]$pslr.results[[i]]$summary)[,1:14]
              rownames(Output[[m]][[i]])[j+1] = paste("Test_", i, "_Subset_", j, sep = "")
            } else {
              Output[[m]][[i]][j+1,] = data.frame(object[[m]]$pslr.results[[i]][[j]]$summary)[,1:14]
              rownames(Output[[m]][[i]])[j+1] = paste("Test_", i, "_Subset_", j, sep = "")
            }
          }
        }
      } # I loop
    } # m loop
    return(Output)
  } else {
    base::summary(object)
  }
} # end of function
