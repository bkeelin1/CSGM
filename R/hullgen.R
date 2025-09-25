#' @title hullgen: an interactive group convex hull generator for Principal Component Analysis
#'
#' @description
#' This function generates a convex hull to visualize group dynamics in a
#' principal component, or similar, space. hullgen utilizes the workhorse
#' package "plotly" to generate interactive plots with the addition of 3D and 2D convex hulls.
#'
#' @param Data_PCs a data table containing principal component scores for each observation.
#'
#' @param select_PC a vector indicating which principal components are to be compared.
#'
#' @param Group a vector containing the grouping values for each individual.
#'
#' @param color a vector containing the color values for each group.
#'
#' @param ID an optional vector indicating the individual names for each specimen in the sample
#'
#' @param PCA_var an optional data table containing three columns: the principal
#' component number, the explained percentage of variance for each principal component,
#' and the cumulative explained percentage of variance. This is automatically
#' produced by the *pca_variances* function.
#'
#' @param sym a boolean indicating whether groups should take on different symbols
#' in the interactive plot.
#'
#' @details
#' This function generates an interactive convex hull plot to highlight group
#' differences between variables or principal components. This plot uses the
#' *plotly* function from the **plotly** package to generate the visualization.
#' To generate the indices used to construct the convex hull, this function uses
#' the *chull* function from the **grDevices** package for 2D convex hulls and the
#' *convhulln* function from the **geometry** package for 3D convex hulls.
#'
#' @returns a plotly object containing a plot with group convex hulls.
#'
#' @author Keeling et al., 2025
#'
#' @seealso \code{plotly} \code{chull} \code{convhulln}
#'
#' @importFrom plotly plot_ly add_trace layout add_polygons
#' @importFrom grDevices chull
#' @importFrom geometry convhulln
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr filter
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' # Example 1: Convex hull generation of principal component data from mandibular corpus shape
#'
#' # Import the sample 3D p x k x n array *Corpus_Land*
#' data(Corpus_Land)
#'
#' # Import grouping variable Collection from the Biomechanics dataset
#' data(Biomechanics)
#' Group = Biomechanics[1:119,]$Collection
#' ID = unique(Biomechanics$ID)
#'
#' # Conduct a Principal Component Analysis on the 3D array using *gm.prcomp* from **geomorph**
#' Data_PCA = gm.prcomp(Corpus_Land)
#' Data_PCs = Data_PCA$x # extract only principal component scores
#'
#' # Summarize the results with the pca_variances function
#' PCA_var = pca_variances(Data_PCA) # input object from the PCA
#'
#' # Generate a PCA plot of PCs 1-2 with a convex hull based on Collection
#' output = hullgen(Data_PCs,
#'                  Group,
#'                  color = NULL, # allow color autofill for each group
#'                  ID,
#'                  PCA_var,
#'                  sym = FALSE # make all groups share the same symbol (circle)
#'                  )
#'
#' # visualize the plot
#' output
#'
#' # Save the plot
#' saveWidget(output, file = "output.html")
#' }

hullgen <- function(Data_PCs,
                    select_PC,
                    Group,
                    color = NULL,
                    ID = NULL,
                    PCA_var = NULL,
                    sym = FALSE) {

  if(length(select_PC) == 3) {
    Data_save = Data_PCs
    Data_PCs = Data_PCs[,select_PC]
    colnames(Data_PCs) <- c("x", "y", "z")
  } else {
    Data_save = Data_PCs
    Data_PCs = Data_PCs[,select_PC]
    colnames(Data_PCs) <- c("x", "y")
  }
  Data_PCs <- data.frame(Data_PCs, Group, ID)

  if(is.null(color)) {
    palette <- c(brewer.pal(length(unique(Group)), "Set2"))
    if(length(unique(Group)) <= 2){
      palette = palette[-3]
    }
    if(length(unique(Group)) >= 2){
      palette = c(palette, palette, palette, palette, palette)
      palette = palette[1:length(unique(Group))]
    }
    group_val = unique(Group)
    group_ind = match(group_val, Group)
    Group_names = group_val[order(group_ind)]
    color <- palette[match(as.factor(Group), Group_names)]
  } else {
    palette = color
    group_val = unique(Group)
    group_ind = match(group_val, Group)
    Group_names = group_val[order(group_ind)]
    color <- palette[match(as.factor(Group), Group_names)]
  }
  if(is.null(ID)) {
    ID <- c(rep(1:nrow(Data_PCs)))
  }
  if(isTRUE(sym)) {
    if(length(unique(Group)) == 1) {
      sym = "circle"
    }
    if(length(unique(Group)) == 2) {
      sym = c("circle", "square")
    }
    if(length(unique(Group)) == 3) {
      sym = c("circle", "square", "diamond")
    }
    if(length(unique(Group)) == 4) {
      sym = c("circle", "square", "diamond", "cross")
    }
    if(length(unique(Group)) == 5) {
      sym = c("circle", "square", "diamond", "cross", "x")
    }
    if(length(unique(Group)) >= 6) {

      sym = c("circle", "square", "diamond", "cross", "x",
              "circle", "square", "diamond", "cross", "x",
              "circle", "square", "diamond", "cross", "x",
              "circle", "square", "diamond", "cross", "x",
              "circle", "square", "diamond", "cross", "x",
              "circle", "square", "diamond", "cross", "x",
              "circle", "square", "diamond", "cross", "x",
              "circle", "square", "diamond", "cross", "x",
              "circle", "square", "diamond", "cross", "x",
              "circle", "square", "diamond", "cross", "x")

      sym = sym[1:length(unique(Group))]

      print("Too many groups for autogeneration. Repeats will be used. Please stop analysis and use sym argument to select your shapes.")
    }
  }
  convex_groups = list()

  if ("z" %in% colnames(Data_PCs)) {

    plot = plot_ly()
    for (i in 1:length(unique(Group))) {

      data = Data_PCs[Data_PCs$Group == Group_names[i],]
      color = rep(palette[i], times = nrow(data))
      symbol = rep(sym[i], times = nrow(data))

      if(isTRUE(sym) || !isFALSE(sym)) {
        symbol = rep(sym[i], times = nrow(data))
        plot = plot %>% add_trace(data = data, x = ~x, y = ~y, z = ~z,
                                  type = "scatter3d", mode = "markers", name = Group_names[i], text = ~ID,
                                  marker = list(size = 10, color = color, symbol = symbol))
      } else {
        plot = plot %>% add_trace(data = data, x = ~x, y = ~y, z = ~z,
                                  type = "scatter3d", mode = "markers", name = Group_names[i], text = ~ID,
                                  marker = list(size = 10, color = color))
      }


      palette2 = data.frame(palette, Group = Group_names)
      current_group2 <- palette2[i,]$Group
      Group_data2 <- Data_PCs[Group == current_group2, ]
      row.names(Group_data2) <- rep(1:nrow(Group_data2))
      colnames(Group_data2) <- c("x","y","z","Group","ID")
      indices <- tryCatch({geometry::convhulln(Group_data2[,c("x","y","z")])}, error = function(e) {
        tryCatch({geometry::convhulln(Data_save[,1:4])}, error = function(e) {
          geometry::convhulln(Data_save[,1:5])
        })
      })
      if(ncol(indices) > 3) {
        indices = indices[,1:3]
      }

      hull <- data.frame(Group_data2[indices,])
      index <- data.frame(x = indices[,1], y = indices[,2], z = indices[,3])
      #palette[match(as.factor(Group), Group_names)]
      #color_hull = palette[match(hull, Group_names)]

      plot <- plot %>% add_mesh(x = hull$x, y = hull$y, z = hull$z,
                                name = paste("Hull", current_group2),
                                showlegend = TRUE,
                                type = 'mesh3d',
                                alphahull = 0,
                                inherit = FALSE,
                                opacity = 0.4,
                                facecolor = rep(paste(palette2[i,1]), times = nrow(hull)))

    } # end of i for loop
    if (!is.null(PCA_var)) {
      PC_total = PCA_var[select_PC[1],2] + PCA_var[select_PC[2],2] + PCA_var[select_PC[3],2]
      plot = plot %>% layout(title = paste("PCA Plot of Principal Components",
                                           select_PC[1], "-",
                                           select_PC[3],"(",
                                           round(PC_total, 4), "%)"), scene = list(
                                             xaxis = list(title = paste("PC", select_PC[1], "(",
                                                                        round(PCA_var[select_PC[1],2],4), "%)")),
                                             yaxis = list(title = paste("PC", select_PC[2] ,"(",
                                                                        round(PCA_var[select_PC[2],2],4), "%)")),
                                             zaxis = list(title = paste("PC", select_PC[3], "(",
                                                                        round(PCA_var[select_PC[3],2],4), "%)"))))
    } else {
      plot = plot %>% layout(title = paste("PCA Plot of Principal Components",
                                           select_PC[1], "-", select_PC[3]),
                             scene = list(
                               xaxis = list(title = paste("PC", select_PC[1])),
                               yaxis = list(title = paste("PC", select_PC[2])),
                               zaxis = list(title = paste("PC", select_PC[3]))))
    }

  } else { # for 2D convex hulls

    plot = plot_ly()

    for (i in 1:length(unique(Group))) {

      data = Data_PCs[Data_PCs$Group == Group_names[i],]
      color = rep(palette[i], times = nrow(data))

      if(isTRUE(sym) || !isFALSE(sym)) {
        symbol = rep(sym[i], times = nrow(data))
        plot = plot %>% add_trace(data = data, x = ~x, y = ~y,
                                  type = "scatter", mode = "markers", name = Group_names[i], text = ~ID,
                                  marker = list(size = 10, color = color, symbol = symbol))
      } else {
        plot = plot %>% add_trace(data = data, x = ~x, y = ~y,
                                  type = "scatter", mode = "markers", name = Group_names[i], text = ~ID,
                                  marker = list(size = 10, color = color))
      }
      palette2 = data.frame(palette, Group = Group_names)
      current_group2 <- palette2[i,2]
      Group_data2 <- Data_PCs[Group == current_group2, ]
      row.names(Group_data2) <- rep(1:nrow(Group_data2))
      colnames(Group_data2) <- c("x","y","Group","ID")
      indices <- chull(Group_data2[, c("x", "y")])
      hull <- data.frame(Group_data2[indices,])

      colnames(hull) <- c("x","y")

      # Add density contour for the group
      plot <- plot %>% add_polygons(plot, x = hull[,1], y = hull[,2], fillcolor = palette2[i,1],
                                    fill = "toself",
                                    name = paste("Hull", current_group2),
                                    opacity = 0.4,
                                    line = list(color = palette2[i,1], width = 2))


    }

    if (!is.null(PCA_var)) {
      PC_total = PCA_var[select_PC[1],2] + PCA_var[select_PC[2],2]
      plot <- plot %>% layout(title = paste("PCA Plot of Principal Components",select_PC[1], "-",
                                            select_PC[2],"(", round(PC_total, 4), "%)"),
                              xaxis = list(title = paste("PC", select_PC[1], "(",
                                                         round(PCA_var[select_PC[1],2],4), "%)")),
                              yaxis = list(title = paste("PC", select_PC[2] ,"(",
                                                         round(PCA_var[select_PC[2],2],4), "%)")))
    } else {
      plot <- plot %>% layout(title = paste("PCA Plot of Principal Components",
                                            select_PC[1], "-", select_PC[2]),
                              xaxis = list(title = paste("PC", select_PC[1])),
                              yaxis = list(title = paste("PC", select_PC[2])))
    }




  } # end of 2D convex construction
  return(plot)
}# End of hull_gen function
