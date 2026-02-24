#' @title SCG: Shape Changes by Group
#'
#' @description
#' This function takes a landmark configuration and a list of
#' grouping variables, calculating each group mean shape for each grouping vector
#' in the list. Then, the function automatically compiles each of the group mean
#' shapes together into a single interactive visualization as an html file stored
#' in the working directory.
#'
#' @param Data Either a p x k x n (landmarks, dimensions, observations) array or
#' a GMA object.
#' @param Groups A list object containing vectors of grouping variables.
#' @param color Optional vector containing the color values to be used for each group.
#' The number of colors should be equivalent to the longest grouping vector in the
#' Groups list.
#' @param coord_name an optional vector containing the names for each landmark
#' in the landmark configuration.
#' @param subset An optional vector specifying the subset or group names associated
#' with each landmark in the landmark configuration. If provided, the resulting
#' interactive figure will allow toggling between these landmark subsets.
#' @param key an optional character string indicating the name of a folder
#' either to be created or already exists in the working directory where generated
#' files can be stored.
#' @param Mesh an optional 3D mesh file (.ply format), preferably the center shape
#' of the landmark array, to be thin plate spline warped into the group mean shapes.
#'
#' @details This function accepts either a p x k x n (landmarks, dimensions, observations)
#' array or a GMA function output. The "Groups" argument should be a list object
#' with one or more vectors. The function will then iterate through each grouping
#' vector in the list and calculate the group mean shapes for each group in the
#' grouping vector.
#'
#' This function allows Grouping vectors to take on character or numeric values.
#' If a vector of numeric values is provided, the minimum and maximum values
#' will be calculated as the "group mean shape".
#'
#' If a 3D mesh is supplied, the landmark array should also be 3D. When a 3D
#' mesh object is supplied, the mesh will be warped using the thin plate spline
#' technique into each group mean shape. Specifically, this function calls the
#' warpRefMesh function from the geomorph package to warp mesh objects to the
#' group mean shape.
#'
#' Subsequently, the function combines each of the group shape means into a single,
#' interactive plot for each grouping variable in the Groups list. These plots
#' are then saved as HTML files in the specified directory. If no string is
#' provided to the argument "key," the function will store these plots in a
#' folder within the working directory called "Shape Differences Between Groups".
#'
#' Interactive visualization plots are generated using the plotly function
#' from the plotly package.
#'
#' @returns A folder titled "Shape Differences Between Groups" in the working
#' directory. This folder will contain plotly objects for each group vector in the
#' Groups argument list. If a mesh object is provided, the function will also store
#' warped meshes for each group value in the grouping vector.
#' \itemize{
#'    \item{Plotly object for grouping vector 1}
#'    \item{Plotly object for grouping vector 2}
#'    \item{Plotly object for grouping vector...i}
#'    \item{Optional PLY Mesh object for grouping vector 1}
#'    \item{Optional PLY Mesh object for grouping vector 2}
#' }
#'
#' @seealso \code{plotly} \code{warpRefMesh}
#'
#' @author Keeling et al., 2025
#'
#' @importFrom geomorph mshape warpRefMesh arrayspecs
#' @importFrom htmlwidgets saveWidget
#' @importFrom dplyr filter
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Scenario:
#' # Identify the mandibular corpus shape differences by skeletal collection.
#' # Identify the mandibular corpus shape differences by mean cortical area
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
#' Groups = list(Collection = Biomechanics[1:119,]$Collection,
#'               'Cortical Area' = Biomechanics %>%
#'                                   group_by(ID) %>%
#'                                     summarise(CA = mean(CA))
#'
#'               )
#'
#' # Perform generalized procrustes analysis on the landmark array
#'
#' gpa_coords = gpagen(Corpus_Land)
#'
#' # Conduct the SCG analysis
#'
#' SCG(gpa_coords # landmark array
#'     Groups # list of grouping vectors
#'     ) # all other arguments are optional. Colors are generated automatically.
#' }
#'

SCG <- function(Data,
                Groups,
                color = NULL,
                coord_name = NULL,
                subset = NULL,
                key = NULL,
                Mesh = NULL,
                save_plots = FALSE) {

  if(isTRUE(save_plots)) {

    Dir <- getwd()
    if(!is.null(key)) {
      if(!dir.exists(paste("Shape Differences Between Groups", key))) {
        dir.create(paste("Shape Differences Between Groups", key))
        temp.dir = paste("Shape Differences Between Groups", key)
      } else {
        temp.dir = paste("Shape Differences Between Groups", key)
      }
    } else {
      if(!dir.exists(paste("Shape Differences Between Groups"))) {
        dir.create(paste("Shape Differences Between Groups"))
        temp.dir = paste("Shape Differences Between Groups")
      } else {
        temp.dir = paste("Shape Differences Between Groups")
      }
    }
    path = file.path(Dir, temp.dir)
    setwd(path)

    on.exit(setwd(Dir))

  }

  if(!is.list(Groups)) {
    Groups = list(Group = Groups)
  }

  if(class(Data) == "GMA") {
    coords = Data$GPA$coords
    Center = Data$GPA$consensus
    k = as.numeric(dim(coords)[2])
  } else {
    coords = Data
    Center = apply(coords, c(1, 2), mean)
    k = as.numeric(dim(coords)[2])
  }

  if(is.null(coord_name)) {
    coord_name = rep(1:dim(coords)[1])
  }

  Group_shapes = vector("list", length(Groups))

  for (j in seq_along(Groups)) {
    if("character" %in% class(Groups[[j]][[1]]) || "factor" %in% class(Groups[[j]][[1]]) ||
       "logical" %in% class(Groups[[j]][[1]])) {
      group_name <- names(Groups)[j]
      if(is.null(group_name) || group_name == "") {
        group_name <- paste0("Group_", j)
      }
      dif_groups <- unlist(na.omit(unique(Groups[[j]])))
      split_groups = list()
      for(g in seq_along(unique(Groups[[j]]))) {
        index = which(Groups[[j]] == unique(Groups[[j]])[g])
        split_groups[[g]] <- coords[,,index,drop = FALSE]
      }
      split_groups <- Filter(function(x) !(is.array(x) && length(x) == 0), split_groups)
      split_groups <- lapply(na.omit(split_groups), geomorph::mshape)
      for(g in seq_along(split_groups)){
        split_groups[[g]] <- array(split_groups[[g]], dim = c(dim(split_groups[[g]])[1], dim(split_groups[[g]])[2], 1))
        if(!is.null(Mesh)) {
          shape = geomorph::warpRefMesh(Center, split_groups[[g]], mesh = Mesh)
          Morpho::mesh2ply(shape, filename = paste(dif_groups[g], "_Shape_", names(Group[[j]]), ".ply", sep = ""))
        }
      }
    } else {
      group_name <- names(Groups)[j]
      if(is.null(group_name) || group_name == "") {
        group_name <- paste0("Group_", j)
      }
      dif_groups <- c("Min", "Max")
      Min = which(Groups[[j]] == min(Groups[[j]], na.rm = TRUE))
      Max = which(Groups[[j]] == max(Groups[[j]], na.rm = TRUE))
      split_groups <- list(Min = coords[,,Min,drop = FALSE], Max = coords[,,Maxdrop = FALSE])
      split_groups <- lapply(split_groups, function(coords) { array(coords, dim = c(dim(split_groups[[1]])[1],
                                                                                    dim(split_groups[[2]])[2], 1))
      })

      split_groups <- lapply(split_groups, function(x) {
        if(dim(x)[3] > 1) geomorph::mshape(x) else x[,,1]
      })

      split_groups <- lapply(split_groups, function(coords) {
        array(coords, dim = c(dim(coords)[1], dim(coords)[2], 1))
      })


      if(!is.null(Mesh)) {
        min.shape = geomorph::warpRefMesh(Center, split_groups$Min, mesh = Mesh)
        Morpho::mesh2ply(min.shape, filename = paste("Min", "_Shape_", group_name, ".ply", sep = ""))
        max.shape = geomorph::warpRefMesh(Center, split_groups$Max, mesh = Mesh)
        Morpho::mesh2ply(max.shape, filename = paste("Max", "_Shape_", group_name, ".ply", sep = ""))
      }
    }
    Group_shapes[[j]] <- lolshape(split_groups, ID = dif_groups, coord_name = coord_name, color = color,
                             title = paste("Shape Differences Between Means of Group:", group_name), subset = subset)

    if(isTRUE(save_plots)){
      saveWidget(Group_shapes[[j]], file.path(path, paste0(group_name, ".html", collapse = "")), selfcontained = TRUE)
    }
  } # end of j for loop

  return(Group_shapes)
} # end of SCG function
