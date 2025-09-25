#' @title lolshape - Am interactive tool to Visualize shape data
#'
#' @description
#' This function allows the interactive comparison of shape based variation between
#' two or more landmark configurations and/or meshes. Interactive graphics are made
#' using functions within the \pkg{plotly} package. This function is best for
#' comparing inter/intra-individual shape comparisons, mean group shape comparisons,
#' and compare shape predictions with the original. This function provides a versatile
#' graphic that permits the visibility of shape configurations and landmark subsetting.
#'
#' @param shape.list a list object containing p x k x n arrays and/or 3D meshes to be visualized.
#'
#' @param ID an optional vector containing the identification for each observation in the landmarked sample.
#'
#' @param title an optional vector to create a custom title for the interactive graphic.
#'
#' @param coord_name a vector listing the names of the landmark coordinates.
#'
#' @param color an optional vector containing the color of each observation in shape.list.
#'
#' @param subset An optional vector specifying the subset or group names associated
#'  with each landmark in the landmark configuration. If provided, the resulting
#'  interactive figure will allow toggling between these landmark subsets.
#'
#' @param mesh a logical value indicating whether shape.list contains at least
#'  one mesh. If mesh is set to TRUE, the function will generate interactive visuals
#'  for the meshes and landmark configurations separately and also as a single combined
#'  interactive visualization.
#'
#' @details
#' This function requires a list of p x k x n arrays and/or 3D mesh objects. Both
#' arrays and meshes can be provided in this list but the arrays must be only either
#' 2D or 3D. Function takes these shape arrays and/or meshes and transforms them
#' into a single, interactive \pkg{plotly} graphic to make interobservation comparisons.
#' If meshes are added, argument mesh must be set to TRUE. If mesh = TRUE, the function
#' will create three outputs: two to separate meshes from landmark arrays and another
#' which combines the two together. Optional arguments are added to allow for landmark
#' subsetting, object colors, plot title, and adding landmark names.
#'
#'
#' @returns a list object containing the following:
#'
#' \itemize{
#'    \item *Graph*:   a plotly object containing only the landmarks listed in shape.list.
#'    \item *Meshes*:  a plotly object containing an interactive visual for all of
#'    the meshes listed in shape.list. Only applicable if mesh = TRUE.
#'    \item *Combined*:  a plotly object containing an interactive visual that
#'    shows both 3D meshes and landmarks in the same shape space. Only applicable
#'    if mesh = TRUE.
#' }
#'
#' @author Keeling et al., 2025
#'
#' @seealso [plotly()]
#'
#' @export
#'
#' @importFrom plotly plot_ly add_trace layout add_mesh
#' @importFrom htmlwidgets saveWidget
#' @importFrom dplyr filter
#' @importFrom grDevices png dev.off
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#'
#' # Import the sample dataset Corpus_Land
#' data(Corpus_Land)
#'
#' # Conduct a Generalized Procrustes Analysis
#' gpa_land = gpagen(Corpus_Land)
#'
#' # Create visualization for the shape comparison between individuals 1212 and 6220
#' output = lolshape(shape.list = list('10' = landmarks[,,"10"],
#'                                     '20' = landmarks[,,"20"]),
#'                   ID = c('10', '20'),
#'                   title = "Shape Comparison Between individuals 10 and 20",
#'                   coord_name = dimnames(Corpus_Land)[[1]],
#'                   color = c("orange", "blue"),
#'                   subset = rep(c("LM1M2", "LP3P4", "Symphysis", "RP3P4", "RM1M2"), each = 120),
#'                   mesh = FALSE)
#'  # Visualize output
#'  output
#'
#'  # Save interactive graphic to computer as html
#'  saveWidget(output, "output.htnl")
#'
#' }

lolshape <- function(shape.list, # a list object containing p x k x n arrays to be visualized
                     ID = NULL, # the observation identification for each landmark configuration
                     title = NULL, # Create a custom title to identify these comparisons.
                     coord_name = NULL, # Add identifiers to each landmark in the configuration.
                     color = NULL, # an optional vector containing the color of each observation in shape.list.
                     subset = NULL, # an optional vector to define subsets of landmarks. Vector should be the same length as the number of total landmarks indicating which
                     mesh = FALSE # a logical value indicating whether a 3D mesh of the landmark configuraiton should be generated.
                     ) {

  coords = list()

  # define k dimensions for shape.list
  array_indices <- which(sapply(shape.list, is.array))
  first_array <- shape.list[[array_indices[1]]]
  k <- dim(first_array)[2]

  # Validate k is 2 or 3 dimensions.
  if (!k %in% c(2, 3)) {
    k = NULL
  }

  if (is.null(coord_name)) {
    n_point = dim(shape.list[[1]])[1]
    coord_name = rep(1:n_point)
    coord_name = as.character(coord_name)
  }
  if(missing(ID)) {
    ID <- c(rep(1:length(shape.list)))
  }
  if(is.null(color)) {
    palette <- suppressWarnings(brewer.pal(length(ID), "Set2"))
    if(length(palette) > length(ID)) {
      palette = if(length(ID) == 1) palette[1] else palette[1:length(ID)]
    } else if (length(ID) <= 65) {
      palette = c(brewer.pal(8, "Set2"),
                  brewer.pal(12, "Set3"),
                  brewer.pal(9, "Pastel1"),
                  brewer.pal(8, "Pastel2"),
                  brewer.pal(12, "Paired"),
                  brewer.pal(8, "Dark2"),
                  brewer.pal(8, "Set1"))

    } else if (length(ID) > 65) {
      n_colors <- length(ID)
      hues <- seq(0, 360 * (n_colors - 1)/n_colors, length.out = n_colors)
      palette <- hcl(h = hues,
                     c = 50,
                     l = 50)
    }
  } else {
    palette = color
  }


  for (i in 1:length(ID)) {

    if(length(dim(shape.list[[i]])) == 2) {
      k = if(dim(shape.list[[i]])[2] == 3) 3 else 2
      shape.list[[i]] <- array(shape.list[[i]], dim = c(dim(shape.list[[i]])[1], dim, 1))
    }

    if(!is.null(subset)) {

      if(k == 3) {
        coords[[i]] <- data.frame(X = shape.list[[i]][,1,], Y = shape.list[[i]][,2,], Z = shape.list[[i]][,3,],
                                  coord_name = coord_name, subset = subset)
      } else {
        coords[[i]] <- data.frame(X = shape.list[[i]][,1,], Y = shape.list[[i]][,2,],
                                  coord_name = coord_name, subset = subset)
      }

    } else {
      if(k == 3) {
        coords[[i]] <- data.frame(X = shape.list[[i]][,1,], Y = shape.list[[i]][,2,], Z = shape.list[[i]][,3,],
                                  coord_name = coord_name)
      } else {
        coords[[i]] <- data.frame(X = shape.list[[i]][,1,], Y = shape.list[[i]][,2,],
                                  coord_name = coord_name)
      }
    }
  } # end of ID loop


  if(is.null(k)) {
    for (i in 1:length(ID)) {
      x = 0
      mesh.list = vector("list", length = sum(x, for(i in seq_along(shape.list)) {x[i] = isTRUE("vb" %in% names(shape.list[[i]]))}))
      for(i in seq_along(shape.list)) {mesh.list[[i]] <- shape.list[[i]][isTRUE("vb" %in% names(shape.list[[i]]))]}
      for(i in seq_along(shape.list)) {names(mesh.list)[i] <- names(shape.list)[i][isTRUE("vb" %in% names(shape.list[[i]]))]}
      for(i in seq_along(mesh.list)) {mesh.list[[i]]$vb <- t(mesh.list[[i]]$vb[1:3, ])
      mesh.list[[i]]$it <- t(mesh.list[[i]]$it[1:3, ]) - 1
      mesh.list[[i]]$Graph <- plot_ly()}
      for(i in seq_along(mesh.list)) {
        mesh.list[[i]]$Graph = mesh.list[[i]]$Graph %>% add_trace(data = mesh.list[[i]],
                                                                  type = "mesh3d", name = names(mesh.list)[i],
                                                                  x = ~vb[, 1],
                                                                  y = ~vb[, 2],
                                                                  z = ~vb[, 3],
                                                                  i = ~it[, 1],
                                                                  j = ~it[, 2],
                                                                  k = ~it[, 3],
                                                                  facecolor = rep("#e8dac9", length(mesh.list[[i]]$it[1,])),
                                                                  opacity = 1)
        mesh.list[[i]]$Graph <- mesh.list[[i]]$Graph %>% layout(title = paste(title), scene = list(
          xaxis = list(showgrid = FALSE, zeroline = FALSE,
                       showticklabels = FALSE, title = list(text = "")),
          yaxis = list(showgrid = FALSE, zeroline = FALSE,
                       showticklabels = FALSE, title = list(text = "")),
          zaxis = list(showgrid = FALSE, zeroline = FALSE,
                       showticklabels = FALSE, title = list(text = "")),
          dragmode = "turntable"))
      }
    }
  } else if (k == 3) {

    if(is.null(subset)) {
      mean.shape = data.frame(X = rowMeans(sapply(coords, function(coords) coords[, 1])),
                              Y = rowMeans(sapply(coords, function(coords) coords[, 2])),
                              Z = rowMeans(sapply(coords, function(coords) coords[, 3])),
                              coord_name = coord_name)
    } else {
      mean.shape = data.frame(X = rowMeans(sapply(coords, function(coords) coords[, 1])),
                              Y = rowMeans(sapply(coords, function(coords) coords[, 2])),
                              Z = rowMeans(sapply(coords, function(coords) coords[, 3])),
                              coord_name = coord_name, subset = subset)
    }

    if(is.null(subset)) {
      Graph = plot_ly()
      Graph <- Graph %>% add_trace(data = mean.shape,
                                   x = ~X, y = ~Y, z = ~Z, text = ~ coord_name,
                                   type = 'scatter3d', mode = 'markers', name = "Mean Shape",
                                   marker = list(color = "black", size = 6))
      for(i in seq_along(ID)) {
        Graph <- Graph %>%
          add_trace(data = coords[[i]],
                    x = ~X, y = ~Y, z = ~Z, text = ~coord_name,
                    type = 'scatter3d', mode = 'markers', name = ID[i], color = palette[i],
                    marker = list(size = 6))
      }
      if(length(shape.list) == 2) {
        x_c <- c(rbind(coords[[1]]$X, coords[[2]]$X))
        y_c <- c(rbind(coords[[1]]$Y, coords[[2]]$Y))
        z_c <- c(rbind(coords[[1]]$Z, coords[[2]]$Z))

        line_name = rep(coord_name, times = 2)

        Graph <- Graph %>%
          add_trace(x = x_c, y = y_c, z = z_c, text = line_name,
                    type = 'scatter3d', mode = 'lines', name = "Mean Shape Differences",
                    line = list(color = "black", width = 2))
      }

      if(length(shape.list) >= 3) {

        for (i in seq_along(coords)) {
          x_c <- c(rbind(coords[[i]]$X, mean.shape$X))
          y_c <- c(rbind(coords[[i]]$Y, mean.shape$Y))
          z_c <- c(rbind(coords[[i]]$Z, mean.shape$Z))
          line_name = rep(coord_name, times = 2)

          Graph <- Graph %>%
            add_trace(x = x_c, y = y_c, z = z_c,
                      type = 'scatter3d', mode = 'lines', name = ID[i], text = line_name,
                      color = palette[i], line = list(width = 2))
        }
      }
      Graph <- Graph %>% layout(title = paste(title), scene = list(
        xaxis = list(showgrid = FALSE, zeroline = FALSE,
                     showticklabels = FALSE, title = list(text = "")),
        yaxis = list(showgrid = FALSE, zeroline = FALSE,
                     showticklabels = FALSE, title = list(text = "")),
        zaxis = list(showgrid = FALSE, zeroline = FALSE,
                     showticklabels = FALSE, title = list(text = "")),
        dragmode = "turntable"))

    } else { # if there are subsets

      subset_markers <- function(graph, data, title = NULL, group_name = "Group",
                                 base_color = "black", extra = NULL) {

        if(is.null(title)) {
          title = "Subsetted Shapes"
        }
        regions <- unique(data$subset)
        for (i in seq_along(regions)) {
          region_data <- filter(data, subset == regions[i])
          graph <- graph %>%
            add_trace(data = region_data,
                      x = ~X, y = ~Y, z = ~Z,
                      type = 'scatter3d', mode = 'markers', text = ~coord_name,
                      marker = list(size = 4, color = base_color),
                      name = paste(group_name, regions[i]))

          if(!is.null(extra)) {
            line_data <- filter(extra, subset == regions[i])
            graph <- graph %>% add_trace(data = line_data,
                                         x = ~X, y = ~Y, z = ~Z,
                                         type = 'scatter3d', mode = 'lines', text = ~coord_name,
                                         line = list(width = 2, color = base_color),
                                         name = paste(group_name, regions[i], "lines"))
          }
        }
        graph <- graph %>% layout(title = paste(title), scene = list(
          xaxis = list(showgrid = FALSE, zeroline = FALSE,
                       showticklabels = FALSE, title = list(text = "")),
          yaxis = list(showgrid = FALSE, zeroline = FALSE,
                       showticklabels = FALSE, title = list(text = "")),
          zaxis = list(showgrid = FALSE, zeroline = FALSE,
                       showticklabels = FALSE, title = list(text = "")),
          dragmode = "turntable"))
        return(graph)
      } # end of subset function

      Graph = plot_ly()
      mean.connection = list()
      for (i in 1:length(ID)) {
        mean.connection[[i]] <- data.frame(X = c(rbind(coords[[i]]$X, mean.shape$X)), Y = c(rbind(coords[[i]]$Y, mean.shape$Y)),
                                           Z = c(rbind(coords[[i]]$Z, mean.shape$Z)), coord_name = rep(coord_name, times = 2),
                                           subset = rep(subset, times = 2))
      } # end of i for loop

      Graph <- subset_markers(Graph, mean.shape, title = NULL, group_name = "Mean Shape", base_color = "black")

      for(i in 1:length(ID)) {
        Graph <- subset_markers(Graph, coords[[i]], title = title, group_name = ID[i],
                                base_color = palette[i], extra = mean.connection[[i]])
      } # end of i for loop

    } # end of subsets if statement

  } else { # End of 3D shape comparisons and beginning of 2D shape comparisons

    if(is.null(subset)) {
      coords[[i]] <- data.frame(X = shape.list[[i]][,1,], Y = shape.list[[i]][,2,], coord_name = coord_name)

      Graph = plot_ly()
      for(i in seq_along(shape.list)) {
        Graph <- Graph %>%
          add_trace(data = coords[[i]] ,
                    x = ~X, y = ~Y, text = ~ coord_name,
                    type = 'scatter', mode = 'markers',
                    marker = list(size = 6, color = palette[i]),
                    name = ID[i])
      }
      if(length(shape.list) == 2){
        x_c <- c(rbind(coords[[1]]$X, coords[[2]]$X))
        y_c <- c(rbind(coords[[1]]$Y, coords[[2]]$Y))

        Graph <- Graph %>%
          add_trace(x = x_pos_vals, y = y_pos_vals, text = coord_name,
                    type = 'scatter', mode = 'lines', name = "Shape Differences",
                    line = list(color = "black", width = 2))
      }
      if(length(shape.list) >= 3){
        mean.shape = data.frame(X = rowMeans(sapply(coords, function(coords) coords[, 1])),
                                Y = rowMeans(sapply(coords, function(coords) coords[, 2])))

        for (i in seq_along(shape.list)){

          x_c <- c(rbind(coords[[i]]$X, mean.shape$X))
          y_c <- c(rbind(coords[[i]]$Y, mean.shape$Y))

          Graph <- Graph %>%
            add_trace(x = x_c, y = y_c, text = coord_name,
                      type = 'scatter', mode = 'lines', name = ID[i],
                      line = list(color = palette[i], width = 2))
        }
      }
      Graph <- Graph %>% layout(title = paste(title), scene = list(
        xaxis = list(showgrid = FALSE, zeroline = FALSE,
                     showticklabels = FALSE, title = list(text = "")),
        yaxis = list(showgrid = FALSE, zeroline = FALSE,
                     showticklabels = FALSE, title = list(text = ""))))

    } else { # if there are subsets

      if(is.null(color)) {
        palette <- c(brewer.pal(length(ID), "Set2"))
      } else {
        palette = color
      }
      if(missing(ID)) {
        ID <- c(rep(1:length(shape.list)))
      }

      subset_markers <- function(graph, data, title = NULL, group_name = "Group",
                                 base_color = "black", extra = NULL) {

        if(is.null(title)) {
          title = "Subsetted Data"
        }
        regions <- unique(data$subset)
        for (i in seq_along(regions)) {
          region_data <- filter(data, subset == regions[i])
          graph <- graph %>%
            add_trace(data = region_data,
                      x = ~X, y = ~Y,
                      type = 'scatter', mode = 'markers', text = ~coord_name,
                      marker = list(size = 4, color = base_color),
                      name = paste(group_name, regions[i]))

          if(!is.null(extra)) {
            line_data <- filter(extra, subset == regions[i])
            graph <- graph %>% add_trace(data = line_data,
                                         x = ~X, y = ~Y,
                                         type = 'scatter', mode = 'lines', text = ~coord_name,
                                         line = list(width = 2, color = base_color),
                                         name = paste(group_name, regions[i], "lines"))
          }
        }
        graph <- graph %>% layout(title = paste(title), scene = list(
          xaxis = list(showgrid = FALSE, zeroline = FALSE,
                       showticklabels = FALSE, title = list(text = "")),
          yaxis = list(showgrid = FALSE, zeroline = FALSE,
                       showticklabels = FALSE, title = list(text = ""))))

        return(graph)
      } # end of subset function

      Graph = plot_ly()

      for (i in 1:length(ID)) {
        coords[[i]] <- data.frame(X = shape.list[[i]][,1,], Y = shape.list[[i]][,2,],
                                  coord_name = coord_name, subset = subset)
      } # end of i for loop

      mean.shape = data.frame(X = rowMeans(sapply(coords, function(coords) coords[, 1])),
                              Y = rowMeans(sapply(coords, function(coords) coords[, 2])),
                              coord_name = coord_name, subset = subset)
      mean.connection = list()
      mean.connection[[i]] <- data.frame(X = c(rbind(coords[[i]]$X, mean.shape$X)), Y = c(rbind(coords[[i]]$Y, mean.shape$Y)))

      Graph <- subset_markers(Graph, mean.shape, "Mean Shape", "black")
      for(i in 1:length(ID)) {
        Graph <- subset_markers(Graph, coords[[i]], title = title, group_name = ID[i], base_color = palette[i])
        Graph <- Graph %>% add_markers(data = mean.connection[[i]], x = ~X, y = ~Y, type = "scatter",
                                       mode = "lines", text = ~coord_name, line = list(width = 2, color = palette[i]))
      } # end of i for loop
    } # end of subsets if statement
  } # End of 2D comparisons
  if(isTRUE(mesh)) {
    Lolshape = list()
    Lolshape$Meshes <- mesh.list
    Lolshape$Graph <- Graph

    for(m in seq_along(mesh.list)) {
      Combined = Graph
      Combined = Combined %>% add_mesh(mesh.list[[m]])
    }
    Lolshape$Combined <- Combined
  } else {
    Lolshape = Graph
  }
  return(Lolshape)
} # End of lolshape function
