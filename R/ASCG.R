#' @title ASCG: Allometric Shape Changes by Groups
#'
#' @description
#' This function takes a landmark configuration, a list of
#' grouping variables, and an allometrically influenced variable to calculate
#' the allometrically influenced shape changes for each group across each grouping
#' vector in the group list. The function uses the visualization package **plotly**
#' to generate interactive figures and compiles each of the allometrically influenced
#' shapes by each grouping vector separately as an html file stored in the
#' working directory.
#'
#' @param GMA_object Either a GMA function output OR a list object
#' with the following structure and components:
#' \itemize{
#'    \item{coords      p x k x n (landmarks, dimensions, observations) landmark array}
#'    \item{consensus   p x k matrix of the center (mean) shape}
#'    \item{Csize       numeric vector of allometric variable or any variable of interest}
#' }
#'
#' @param Groups A list object containing vectors of grouping variables.
#'
#' @param size An optional vector for an allometrically influenced proxy variable.
#' By default, the natural log centroid size is used.
#'
#' @param color Optional vector containing the color values to be used for each group.
#' The number of colors should be equivalent to the longest grouping vector in the
#' Group list.
#' @param coord_name an optional vector containing the names for each landmark
#' in the landmark configuration.
#'
#' @param subset An optional vector specifying the subset or group names associated
#' with each landmark in the landmark configuration. If provided, the resulting
#' interactive figure will allow toggling between these landmark subsets.
#'
#' @param Mesh an optional 3D mesh file (.ply format), preferably the center shape
#' of the landmark array, to be thin plate spline warped into the group mean shapes.
#'
#' @param include_whole a boolean value indicating whether to create visuals
#' for the allometric influences on shape variation for the entire sample.
#'
#' @param save_plots a logical value indicating whether to save the plots to a user-specified directory.
#'
#' @param key an optional character string indicating the name of a folder
#' either to be created or already exists in the working directory where generated
#' files can be stored.
#'
#' @param parallel a boolean value indicating whether to use parallel processing
#' of the shape predictions and mesh generations for each group vector in the
#' Groups object. Parallel processing is strongly recommended for long Group lists.
#'
#' @details
#' This function accepts either a GMA function output or a list object
#' with the following structure:
#' \itemize{
#'    \item{coords       p x k x n (landmarks, dimensions, observations) landmark array}
#'    \item{consensus    p x k matrix of the center (mean) shape}
#'    \item{Csize        numeric vector of allometric variable or any variable of interest}
#' }
#'
#' The coords object is used as the landmark configuration for all observations
#' in the entire sample. The consensus object is the center shape of all observations
#' in the entire sample. The Csize object is the allometrically influenced variable.
#' This function uses a linear Procrustes regression via the *procD.lm* function
#' from the **Geomorph** package. Please see the *procD.lm* function for details
#' on how the regression is performed and the landmark data is subjected to
#' dimensionality reduction.
#'
#' The function first performs a linear Procrustes regression of the landmark array
#' on the allometrically influenced variable for the entire sample if the
#' "include_whole" argument is TRUE. If FALSE, this step will be skipped.
#' Subsequently, and in parallel, each grouping vector in the Groups list will be
#' split from the landmark array as separate, group arrays. For each group value
#' in a grouping vector, these specific grouped observations will be subjected to
#' a linear Procrustes regression on the allometrically influenced variable. For
#' each group the predicted shape for the observation with the highest and smallest
#' values from allometrically influenced variable will be compiled into a single
#' interactive plotly plot html file in the specified directory folder. Each file
#' is automatically labeled with the phrase "Allometric Shape Changes_" followed
#' by the name of the unique grouping value and ended with ".html". This is
#' performed for each group in parallel if the argument "parallel" is TRUE.
#'
#' The "Group" argument should be a list object with one or more vectors. The
#' function will then iterate through each grouping vector in the list and
#' predict the allometrically influenced shapes for each group in the grouping vector.
#' This function allows Grouping vectors to take on character or numeric values.
#' If a vector of numeric values is provided, the predicted allometric shapes
#' for the minimum and maximum values will be used for the *plotly* comparisons.
#'
#' If a 3D mesh is supplied, the landmark array should also be 3D. When a 3D
#' mesh object is supplied, the mesh will be warped using the thin plate spline
#' technique into each allometrically predicted group shape. Specifically, this
#' function calls the warpRefMesh function from the geomorph package to warp
#' mesh objects to the group mean shape.
#'
#' @returns
#' A folder titled "Allometric Influences by Group" in the working
#' directory. This folder will contain plotly objects for each group vector in the
#' Group argument list. If a mesh object is provided, the function will also store
#' warped meshes for each group value in the grouping vector.
#' \itemize{
#'    \item{Group            list object each containing a plotly plot object}
#'    \item{Subsetted_Group  list object each containing a plotly plot object with subsetted landmarks}
#' }
#'
#' @author Brian Anthony Keeling
#'
#' @seealso \code{\link{plotly}} \code{\link{gpagen}} \code{\link{GMA}}
#'
#' @importFrom geomorph gpagen arrayspecs mshape warpRefMesh
#' @importFrom Morpho mesh2ply
#' @importFrom grDevices png dev.off
#' @importFrom htmlwidgets saveWidget
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Scenario: Identify the allometrically influenced shape changes of the mandibular
#' # corpus by skeletal collection and mean cortical area
#'
#' # Import the sample dataset Corpus_Land
#'
#' data(Corpus_Land)
#'
#' # Import the sample dataset Biomechanics
#'
#' data(Biomechanics)
#'
#' # Create a list of grouping vectors for the variables Collection and Cortical Area (CA)
#'
#' Group = list(Collection = Biomechanics[1:119,]$Collection,
#'               'Cortical Area' = Biomechanics %>%
#'                                   group_by(ID) %>%
#'                                     summarise(CA = mean(CA))
#'
#'               )
#'
#' # Option 1: Perform generalized procrustes analysis on the landmark array
#'
#' GMA_object = gpagen(Corpus_Land)
#'
#' # Option 2: Perform the GMA analysis
#'
#' GMA_object = GMA(Corpus_Land, Groups = Group)
#'
#' # Conduct the SCG analysis
#'
#' ASCG(GMA_object # gpagen or GMA object output
#'     Group # list of grouping vectors
#'     ) # all other arguments are optional. Colors are generated automatically.
#' }

ASCG <- function(GMA_object,
                 Groups,
                 color = NULL,
                 coord_name = NULL,
                 subset = NULL,
                 Mesh = NULL,
                 include_whole = TRUE,
                 save_plots = TRUE,
                 key = NULL,
                 parallel = TRUE) {

  if(isTRUE(save_plots)) {

    Dir <- getwd()
    if(!is.null(key)) {
      if(!dir.exists(paste("Group Allometry", key))) {
        dir.create(paste("Group Allometry", key))
        temp.dir = paste("Group Allometry", key)
      } else {
        temp.dir = paste("Group Allometry", key)
      }
    } else {
      if(!dir.exists(paste("Group Allometry"))) {
        dir.create(paste("Group Allometry"))
        temp.dir = paste("Group Allometry")
      } else {
        temp.dir = paste("Group Allometry")
      }
    }
    path = file.path(Dir, temp.dir)
    setwd(path)
  }
  #_______________________________________________________________________

  if(class(GMA_object) == "GMA") {
    coords = GMA_object$GPA$coords
    Center = GMA_object$GPA$consensus
    Size = GMA_object$GPA$Csize
    k = as.numeric(dim(coords)[2])
  } else {
    coords = GMA_object$coords
    Center = GMA_object$consensus
    Size = GMA_object$Csize
    k = as.numeric(dim(coords)[2])
  }

  if(is.null(coord_name)) {
    coord_name = rep(1:dim(coords)[1])
  }
  #___________________________________________________________________________
  # Allo shape without group influence
  if(isTRUE(include_whole)) {
    regressions = list()
    Data = list(coords, Size)
    reg <- summary(procD.lm(coords ~ log(Size), data = Data))
    regressions$Allo.size.reg <- reg
    reg_data <- procD.lm(coords ~ log(Size), data = Data)
    if(!is.null(subset)) {
      for(s in seq_along(unique(subset))) {
        subdex = which(subset == unique(subset)[[s]])
        sub_data = list(coords = coords[subdex,,], Size = Size)
        sub_reg_results <- summary(procD.lm(coords ~ log(Size), data = sub_data))

        regressions$Subset_Reg[[s]] <- sub_reg_results[[1]]
        names(regressions$Subset_Reg)[s] <- paste("Total.Sample", ".Subset.", unique(subset)[s], sep = "")
      }
    }
    coords.fitted = arrayspecs(reg_data$fitted, dim(coords)[1], k)
    big = which(data.frame(reg_data$X)[2] == max(data.frame(reg_data$X)[2]))
    small = which(data.frame(reg_data$X)[2] == min(data.frame(reg_data$X)[2]))
    allo.shape = list(Big = coords.fitted[,,big], Small = coords.fitted[,,small])

    if(!is.null(Mesh)) {
      big.shape = warpRefMesh(Center, allo.shape$Big, mesh = Mesh)
      mesh2ply(big.shape, filename = paste("Large Ind", "_Allometric Shape Change", ".ply", sep = ""))
      graphics.off()
      rm(big.shape)
      small.shape = warpRefMesh(Center, allo.shape$Small, mesh = Mesh)
      mesh2ply(small.shape, filename = paste("Small Ind", "_Allometric Shape Change", ".ply", sep = ""))
      graphics.off()
      rm(small.shape)
    }

    allo.shape = lapply(allo.shape, function(x) { array(x, dim = c(dim(x)[1], dim(x)[2], 1)) })

    Size.Shape = lolshape(allo.shape, ID = c(paste("Large Ind:", big), paste("Small Ind:", small)), color = color,
                            coord_name = coord_name, title = "Allometrically Influenced Shape Changes",
                            subset = subset)
    if(isTRUE(save_plots)) {
      saveWidget(Size.Shape, file = "Allometrically Influenced Shape Changes.html", selfcontained = TRUE)
    }
  } # end of if including whole allometric analysis
  #___________________________________________________________________________
  # Allometry by group
  if(isTRUE(parallel)) {
    doParallel::registerDoParallel(cores = detectCores() - 2)
  }

  allo.parallel = foreach(j = 1:length(Groups), .combine = 'rbind', .export = c("lolshape"), .packages =
                            c("geomorph", "Morpho", "doParallel", "tidyverse", "RColorBrewer", "plotly", "htmlwidgets")) %dopar% {

    if("character" %in% class(Groups[[j]][[1]]) || "factor" %in% class(Groups[[j]][[1]]) ||
       "logical" %in% class(Groups[[j]][[1]])) {
      reg_by_group = vector("list", length(na.omit(unique(Groups[[j]]))))
      reg_results = list()
      r.pred = list()

      for(g in seq_along(na.omit(unique(Groups[[j]])))) {
        index = which(Groups[[j]] == na.omit(unique(Groups[[j]]))[g])
        name = na.omit(unique(Groups[[j]]))[g]
        if(class(GMA_object) == "GMA") {
          coords = GMA_object$GPA$coords[,,index]
          Center = GMA_object$GPA$consensus
          Size = GMA_object$GPA$Csize[index]
        } else {
          coords = GMA_object$coords[,,index]
          Center = GMA_object$consensus
          Size = GMA_object$Csize[index]
        }
        Data = list(coords, Size)
        reg_by_group[[g]] <- summary(procD.lm(coords ~ log(Size), data = Data))
        names(reg_by_group[[g]]) <- paste(names(Groups)[j], ".",na.omit(unique(Groups[[j]]))[g], sep = "")
        r.pred[[g]] <- procD.lm(coords ~ log(Size), data = Data)
        if(!is.null(subset)) {
          subset_reg_by_group = vector("list", length(unique(subset)))
          for(s in seq_along(unique(subset))) {
            subdex = which(subset == unique(subset)[[s]])
            sub_data = list(coords = coords[subdex,,], Size = Size)
            sub_reg_results <- summary(procD.lm(coords ~ log(Size), data = sub_data))
            subset_reg_by_group[[s]][[g]] <- sub_reg_results
            names(subset_reg_by_group[[s]][[g]]) <- paste(names(Groups)[j], ".",
                                                          na.omit(unique(Groups[[j]]))[g], ".",
                                                          unique(subset)[s], sep = "")
          }
        }
        coords.fit = arrayspecs(r.pred[[g]]$fitted, dim(coords)[1], k)
        big = which(data.frame(r.pred[[g]]$X)[2] == max(data.frame(r.pred[[g]]$X)[2]))
        small = which(data.frame(r.pred[[g]]$X)[2] == min(data.frame(r.pred[[g]]$X)[2]))
        allo.list = list(big = coords.fit[,,big], small = coords.fit[,,small])

        if(!is.null(Mesh)) {
          big.shape = warpRefMesh(Mesh, Center, allo.list$big)
          mesh2ply(big.shape, filename = paste("Large Ind", "_Allo Shape_", names(Groups)[j],"_",name,".ply", sep = ""))
          graphics.off()
          small.shape = warpRefMesh(Mesh, Center, allo.list$small)
          mesh2ply(small.shape, filename = paste("Small Ind", "_Allo Shape_", names(Groups)[j],"_",name,".ply", sep = ""))
          graphics.off()
        }

        allo.shape = lapply(allo.list, function(x) { array(x, dim = c(dim(x)[1], dim(x)[2], 1)) })

        allo.shape = lolshape(allo.shape, ID = c(paste("Large Ind:", big), paste("Small Ind:", small)),
                                title = paste("Allometric Shape Changes", names(Groups[[j]]), "_",
                                              name),  subset = subset, color = color)
        if(isTRUE(save_plots)) {
          saveWidget(allo.shape, paste("Allometric Shape Changes", names(Groups)[j], "-",
                     name, ".html", sep = ""), selfcontained = TRUE)
        }
      }
    } else {
      reg_by_group = vector("list", 1)
      index = which(Groups[[j]][[1]] %in% na.omit(Groups[[j]][[1]]))
      if(class(GMA_object) == "GMA") {
        coords = GMA_object$GPA$coords[,,index]
        Center = GMA_object$GPA$consensus
        Size = GMA_object$GPA$Csize[index]
      } else {
        coords = GMA_object$coords[,,index]
        Center = GMA_object$consensus
        Size = GMA_object$Csize[index]
      }
      Data = list(coords = coords, Size = Size, Group = na.omit(Groups[[j]][[1]][index]))
      reg_by_group[[1]] <- summary(procD.lm(coords ~ log(Size) * Group, data = Data))
      names(reg_by_group) <- paste(names(Groups)[j], sep = "")
      r.pred <- procD.lm(coords ~ log(Size) * Group, data = Data)
      coords.fit = arrayspecs(r.pred$fitted, dim(coords)[1], k)
      if(!is.null(subset)) {
        subset_reg_by_group = vector("list", length(unique(subset)))
        for(s in seq_along(unique(subset))) {
          subdex = which(subset == unique(subset)[[s]])
          sub_data = list(coords = coords[subdex,,], Size,
                          Group = na.omit(Groups[[j]])[[1]][index])
          subset_reg_by_group[[s]] <- summary(procD.lm(coords ~ log(Size), data = sub_data))
          names(subset_reg_by_group)[s] <- paste(names(Groups)[j], ".", unique(subset)[s], sep = "")
        }
      }
      big = which(data.frame(r.pred$X)[4] == max(data.frame(r.pred$X)[4]))
      small = which(data.frame(r.pred$X)[4] == min(data.frame(r.pred$X)[4]))
      allo.list = list(big = coords.fit[,,big], small = coords.fit[,,small])

      if(!is.null(Mesh)) {
        big.shape = warpRefMesh(Mesh, Center, allo.list$big)
        mesh2ply(big.shape, filename = paste("Large Ind", "_Allo Shape_", names(Groups)[j], ".ply", sep = ""))
        graphics.off()
        small.shape = warpRefMesh(Mesh, Center, allo.list$small)
        mesh2ply(small.shape, filename = paste("Small Ind", "_Allo Shape_", names(Groups)[j], ".ply", sep = ""))
        graphics.off()
      }

      allo.shape = lapply(allo.list, function(x) { array(x, dim = c(dim(x)[1], dim(x)[2], 1)) })

      allo.shape = lolshape(allo.shape, ID = c(paste("Large Ind:", big), paste("Small Ind:", small)), color = color,
                              title = paste("Allometric Shape Changes", ".", names(Groups[[j]])), subset = subset)

      if(isTRUE(save_plots)) {

        saveWidget(allo.shape, paste("Allometric Shape Changes_", names(Groups)[j], ".html", sep = ""), selfcontained = TRUE)
      }
    }
    if(!is.null(subset)) {
      return(list(reg_by_group = reg_by_group, subset_reg_by_group = subset_reg_by_group))
    } else {
      return(list(reg_by_group = reg_by_group))
    }
  } # end of j foreach parallel loop
  closeAllConnections()
  #___________________________________________________________________________
  place = list()
  for(a in 1:nrow(allo.parallel)) {
    place[[a]] <- allo.parallel[a,1]
    #place[[a]] <- unlist(place[[a]], recursive = FALSE)
    place[[a]] <- Filter(function(x) !is.null(x), place[[a]])
    place[[a]][[1]] <- Filter(function(x) !is.null(x), place[[a]][[1]])
    #place[[a]][[1]] <- unlist(place[[a]][[1]], recursive = FALSE)
    # place[[a]][[1]] <- place[[a]][[1]][class(place[[a]][[1]]) %in% "anova"]
    place[[a]] <- place[[a]][[1]]
    names(place)[a] <- names(Groups[[a]])
  }
  for(a in seq_along(place)) {
    for(b in seq_along(place[[a]])) {
      names(place[[a]])[b] <- names(place[[a]][[b]])[1]
      place[[a]][[b]] <- place[[a]][[b]][[1]]
    }
  }

  names(place) <- names(Groups)

  if(!is.null(subset)) {

    place2 = list()
    for(a in 1:nrow(allo.parallel)) {
      place2[[a]] <- allo.parallel[a,2]
      #place2[[a]] <- unlist(place2[[a]], recursive = FALSE)
      place2[[a]] <- Filter(function(x) !is.null(x), place2[[a]])
      place2[[a]][[1]] <- Filter(function(x) !is.null(x), place2[[a]][[1]])
      #place2[[a]][[1]] <- unlist(place2[[a]][[1]], recursive = FALSE)
      # place2[[a]][[1]] <- place2[[a]][[1]][class(place2[[a]][[1]]) %in% "anova"]
      place2[[a]] <- place2[[a]][[1]]
      names(place2)[a] <- names(Groups[[a]])

    }
    for(a in seq_along(place2)) {
      for(b in seq_along(place2[[a]])) {
        place2[[a]][[b]] <- Filter(function(x) !is.null(x), place2[[a]][[b]])
        names(place2[[a]])[b] <- names(place2[[a]][[b]])[1]
        place2[[a]][[b]] <- place2[[a]][[b]][[1]]
      }
    }
    for(a in seq_along(place2)) {
      for(b in seq_along(place2[[a]])) {
        names(place2[[a]])[b] <- names(place2[[a]][[b]])[1]
        place2[[a]][[b]] <- place2[[a]][[b]][[1]]
      }
    }

    names(place2) <- names(Groups)

    for(j in seq_along(Groups)) {
      for(s in seq_along(subset)) {
        names(place2)[j] <- na.omit(unique(Groups[[j]]))[s]
      }
    }
  }
  regressions$Group <- place
  regressions$Subsetted_Group <- place2
  #_____________________________________________________________________________
  setwd(Dir)
  return(regressions)
}
