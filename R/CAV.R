#' @title CAV: Cluster Analysis Visualizations
#'
#' @description This function conducts both an agglomerative hierarchical
#' cluster analysis and an k-medoids cluster analysis, generating
#' visually interactive figures to evaluate relationships among observations. Both
#' of these cluster analyses are unsupervised, that is, they aren't given group labels,
#' and instead use distance matrices to identify clusters of individual observations.
#' This is useful in contexts, such as geometric morphometrics, in identifying
#' taxonomic or even inter/intra-population differences based on their relationships
#' to one another.
#'
#' @param Data a data table or data frame object.
#' @param Group a list object containing a vector of grouping variables.
#' @param ID an optional vector containing the observation identifications.
#' @param nboot the number of bootstraps to be run for the agglomerative cluster analysis.
#' @param Dir_name an optional character string indicating the name of a folder
#' either to be created or already exists in the working directory where generated
#' files can be stored.
#' @param parallel a boolean value indicating whether to use parallel processing.
#' Parallel processing is strongly recommended in cluster analysis.
#' @param method a character string indicating the distance matrix method that
#' should be applied on the data. Options include: "manhattan",
#' "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower",
#' "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao",
#' "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison".
#' Please see the *vegdist* "method" argument from the **vegan** package for more
#' details.
#'
#' @param ... Optional arguments to be passed into the Predictor arguments for the
#' \code{DT} function. This is used to create distance matrices for the analysis.
#'
#' @details This function accepts a data table and from it generates a distance matrix.
#' A variety of standardized methods for distance matrix calculation are avaialble
#' including mahalanobis distances. For a full list, please view the method argument
#' from the *vegdist* function in the **vegan** package. Mahalanobis distances
#' are calculated based on the *mahalanobis* function from the **stats** package.
#' Subsequently, this distance matrix is subject to a bootstrapped, agglomerative
#' hierarchical cluster analysis via the *pvclust* function from the **pvclust** package.
#' The agglomerative cluster analysis linkage method for bottom-up clustering
#' is automatically identified using the *agnes* function from the **cluster**
#' package by computing yhe agglomerative coefficient (AC), as a measure how well
#' the given data can be clustered for each possible method ("average", "single",
#' "complete", "ward", "weighted"). The method with the highest AC is used. Additionally,
#' a partitioning cluster analysis is conducted through k-medoids using the *pam* function
#' from the **cluster** package. Partitioning Around Medoids (PAM), as opposed to
#' k-means clustering, is more robust to data outliers and uses medoids in its clusterings.
#' First, the number of possible k clusters are iteratively performed using the *pam*
#' function. Then the average silhouette width, an average measure of how well the data
#' fits the cluster, is calculated for each k and the k clusters with the highest
#' average silhouette width is selected as k. Then, the function is run again with
#' the selected parameters.
#'
#' Furthermore, this function publishes each cluster plot into an automatically
#' generated folder whose directory name can be defined through Dir_name. Cluster
#' data from both the agglomerative and pam clustering is output from the function.
#'
#' @returns a list object of the following:
#' \itemize{
#'    \item{k_medoid:   a pam object containing the results of the pam cluster analysis.}
#'    \item{k_medoid_plot:    an integer indicating the number of k-medoid clusters selected.}
#'    \item{dendro:   a phylo object containing the dendrogram from the agglomerative cluster analysis.}
#'    \item{agg_plot:   a list object containing fan plots of the agglomerative clusters for each group}
#' }
#'
#' @author
#' Keeling et al., 2025
#'
#' @seealso [\code{pam} \code{pvclust} \code{agnes}
#'
#' @importFrom vegan vegdist
#' @importFrom stats mahalanobis cov dist
#' @importFrom cluster pam agnes
#' @importFrom pvclust pvclust
#' @importFrom RColorBrewer brewer.pal
#' @importFrom factoextra fviz_cluster
#' @importFrom grDevices png dev.off recordPlot
#' @importFrom ape as.phylo plot.phylo Ntip Nnode
#' @importFrom geomorph two.d.array gm.prcomp
#' @importFrom ggpubr theme_pubclean
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # Example grouping mandibular corpus landmarks by Collection (for agglomerative clustering)
#'
#' # Import the sample 3D p x k x n array *Corpus_Land*
#' data(Corpus_Land)
#'
#' # conduct a generalized procrustes analysis of the data
#'
#' gpa_data = gpagen(Corpus_Land)
#'
#' # Transform the 3D array as a data table using the geomorph package
#' table_corpus_land = data.frame(two.d.array(gpa_data))
#'
#' # Import grouping variable Collection from the Biomechanics dataset
#' data(Biomechanics)
#' Group = Biomechanics[1:119,]$Collection
#' ID = unique(Biomechanics$ID)
#'
#' # Run the cluster analysis vizualization tool
#'
#' output = CAV(Data = table_corpus_land,
#'              Group = Group,
#'              ID = ID,
#'              nboot = 1000,
#'              Dir_name = "Preferred Directory Folder Name",
#'              parallel = TRUE,
#'              method = "mahalanobis"
#'              )
#'
#' output
#' }
#'
CAV <- function(Data,
                Group = NULL,
                ID = NULL,
                nboot = 1000,
                Dir_name = NULL,
                parallel = TRUE,
                method = "mahalanobis",
                ...) {

  Morphotree = list()

  args = list(...)

  if(is.null(nboot)) {
    nboot = 1000
  }

  if(is.null(ID)) {
    ID = c(1:nrow(Data))
  }

  if(is.null(Group)) {
    Group = data.frame(Group = rep("No Group Data", times = nrow(Data)))
  }

  if(!is.null(Dir_name)) {
    original_wd <- getwd()
    if (!dir.exists(paste("Morphograms_", Dir_name, sep = " "))) {
      Directory <- paste("Morphograms_", Dir_name, sep = " ")
      path = file.path(original_wd, Directory)
      dir.create(path)
      path = file.path(original_wd, Directory)
    } else {
      Directory <- paste("Morphograms_", Dir_name, sep = " ")
    }
  } else {
    original_wd <- getwd()
    if (!dir.exists(paste("Morphograms"))) {
      Directory <- paste("Morphograms")
      path = file.path(original_wd, Directory)
      dir.create(path)
      path = file.path(original_wd, Directory)
    } else {
      Directory <- paste("Morphograms")
    }
  }

  if("dist" %in% class(Data)) {
    Euc_dist = Data
  } else {
    dt_params = names(formals(DT))
    dt_parameters = modifyList(list(Predictor = Data,
                                    Pred_transform = "procdist",
                                    dt_method = method), args) %>%
                          .[intersect(names(.), dt_params)]
    Euc_dist = do.call(DT,dt_parameters)$Predictor
  }

  Euc_dist = as.dist(Euc_dist)

  # Hierarchical Cluster Generation
  m <- c( "average", "single", "complete", "ward", "weighted")
  names(m) <- c("average", "single", "complete", "ward.D2", "mcquitty")
  ac <- function(x) {cluster::agnes(Euc_dist, method = x)$ac}
  models = purrr::map_dbl(m, ac)
  choice <- names(which.max(models))
  clust = pvclust::pvclust(as.matrix(Euc_dist), method.hclust = choice, weight = TRUE,
                  nboot = nboot, parallel=parallel, iseed=1, quiet=TRUE)

  # K-Medoids Cluster Generation
  tryCatch({
    if(length(ID) == 10) {
      kchoice = numeric(9)
    } else if (length(ID) < 10) {
      kchoice = numeric(length(ID))
    } else {
      kchoice = numeric(10)
    }

    for(k in 2:length(kchoice)) {
      pam_clust = cluster::pam(Euc_dist, k = k, metric = method, nstart = 25, pamonce = TRUE, diss = TRUE)
      pam_clust$data <- as.matrix(Euc_dist)
      kchoice[k] = pam_clust$silinfo$avg.width
    }
    kchoice = which.max(kchoice)
    pam_clust = cluster::pam(Euc_dist, k = kchoice, metric = method, nstart = 25, pamonce = TRUE, diss = TRUE)
    pam_clust$data <- as.matrix(Euc_dist)

    hclust = clust$hclust
    hclust$labels <- ID
    dendro <- ape::as.phylo(hclust)

    k_plot <- factoextra::fviz_cluster(pam_clust, pointsize = 5, repel = TRUE, alpha = 0.6, label = 8, main = "Results of K-Medoids Cluster Analysis", ggtheme = ggpubr::theme_pubclean())

    if (length(dev.list()) > 0) {
      graphics.off()
    }
    cluster_plot_path <- paste(Directory,"/plot_cluster_", "K-medoids", ".png", sep = "")
    grDevices::png(cluster_plot_path, res = 600, width = 1920, height = 1080)
    # label used to be 8
    print(factoextra::fviz_cluster(pam_clust, pointsize = 5, repel = TRUE, alpha = 0.6, label = 8, main = "Results of K-Medoids Cluster Analysis", ggtheme = ggpubr::theme_pubclean()))
    dev.off()
  }, error = function(e) {
    print("PAM analysis cannot determine the number of clusters. Consider using another distance matrix and review the nature of your data.")
  })

  # Plot Generation

  if (class(Group) == "data.frame" || "data.frame" %in% class(Group)) {
    Group_list = list()
    for(i in seq_along(Group)){
      Group_list[[i]] <- as.factor(Group[,i])
    }
    names(Group_list) <- colnames(Group)
    Group = Group_list
  } else if ("numeric" %in% class(Group) || "character" %in% class(Group) || "factor" %in% class(Group)) {
    Group = list(Group = as.factor(Group))
    Group
  }

  if (class(Group) == "list") {

    #agg_clad_plot = list()
    agg_fan_plot = list()

    for (i in 1:length(Group)) {

      group_name <- names(Group)[i]
      group_val = unique(Group[[i]])
      group_ind = match(group_val, Group[[i]])
      group_val = group_val[order(group_ind)]
      palette <- suppressWarnings(brewer.pal(length(group_val), "Set2"))
      if(length(palette) > length(group_val)) {
        palette = if(length(group_val) == 1) palette[1] else palette[1:length(group_val)]
      } else if (length(group_val) <= 65) {
        palette = c(brewer.pal(8, "Set2"),
                    brewer.pal(12, "Set3"),
                    brewer.pal(9, "Pastel1"),
                    brewer.pal(8, "Pastel2"),
                    brewer.pal(12, "Paired"),
                    brewer.pal(8, "Dark2"),
                    brewer.pal(8, "Set1"))

      } else if (length(group_val) > 65) {
        n_colors <- length(group_val)
        hues <- seq(0, 360 * (n_colors - 1)/n_colors, length.out = n_colors)
        palette <- hcl(h = hues,
                       c = 50,
                       l = 50)
      }

      Group_col <- palette[match(as.factor(Group[[i]]), group_val)]


      edge_colors <- rep("black", nrow(dendro$edge))
      for (i in seq_along(Group_col)) {
        edge_colors[dendro$edge[,2] == i] <- Group_col[i]
      }

      getAllDescendantTips <- function(node, tree) {
        if (node <= Ntip(tree)) {
          return(node)  # If it's a tip, return itself
        }
        children <- tree$edge[, 2][tree$edge[, 1] == node]
        unlist(lapply(children, function(child) getAllDescendantTips(child, tree)))
      }

      getMostCommonColor <- function(colors) {
        valid_colors <- colors[colors != "black"]
        if (length(valid_colors) == 0) return("black")  # Return black if no valid colors
        names(which.max(table(valid_colors)))  # Most common color
      }

      colortree2 <- function(tree, Group_colors) {
        # Initialize edge colors based on tip colors provided by Group_colors
        edge_colors <- rep("black", nrow(tree$edge))
        tip_ids <- 1:Ntip(tree)
        edge_colors[tip_ids] <- Group_colors  # Assuming Group_colors is named by tip labels

        # Apply colors to internal nodes based on the most common descendant tip color
        for (node in (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))) {
          descendant_tips <- getAllDescendantTips(node, tree)
          tip_colors <- edge_colors[descendant_tips]
          most_common_color <- getMostCommonColor(tip_colors)

          # Apply the most common color to the edges leading to this node
          parent_edge <- which(tree$edge[, 2] == node)
          edge_colors[parent_edge] <- most_common_color
        }
        for (i in seq_along(Group_colors)) {
          edge_colors[dendro$edge[,2] == i] <- Group_col[i]
        }
        return(edge_colors)
      }  # end of color tree function

      colortree <- function(tree, Group_colors) {
        # Get tip labels
        tip_labels <- tree$tip.label

        # Initialize edge colors based on tip colors provided by Group_colors
        edge_colors <- rep("black", nrow(tree$edge))

        # Map colors to tips based on their labels
        for (i in 1:Ntip(tree)) {
          if (i <= length(Group_colors)) {
            # Find edges leading to tips
            tip_edges <- which(tree$edge[,2] == i)
            if (length(tip_edges) > 0) {
              edge_colors[tip_edges] <- Group_colors[i]
            }
          }
        }

        # Apply colors to internal nodes based on descendant tips
        for (node in (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))) {
          descendant_tips <- getAllDescendantTips(node, tree)
          if (length(descendant_tips) > 0) {
            tip_colors <- Group_colors[descendant_tips]
            if (length(tip_colors) > 0) {
              most_common_color <- names(which.max(table(tip_colors)))
              parent_edges <- which(tree$edge[,2] == node)
              if (length(parent_edges) > 0) {
                edge_colors[parent_edges] <- most_common_color
              }
            }
          }
        }

        return(edge_colors)
      }

      edge_colors <- colortree(dendro, Group_col)

      if (length(dev.list()) > 0) {
        graphics.off()
      }

      agg_fan_plot[[i]] <- list(x = dendro,
                            edge.color = edge_colors,
                            tip.color = Group_col,
                            edge.width = 3,
                            cex = 0.7,
                            type = "fan",
                            main = "Fan of Morphometric Shape Distances")

      fan_plot_path <- paste0(Directory,"/plot_fan_", group_name, ".png", sep = "")

      grDevices::png(fan_plot_path, res = 600, width = 1920, height = 1080)
        par(mar = c(1, 4, 1, 1))
        plot(dendro,
             edge.color = edge_colors,
             tip.color = Group_col,
             edge.width = 3,
             cex = 0.7,
             type = "fan",
             main = "Fan of Morphometric Shape Distances")
        legend("topright", legend = group_val, col = palette, pch = 19, cex = 0.4)
      dev.off()

      if (length(dev.list()) > 0) {
        graphics.off()
      }

      #cladogram_plot_path <- paste0(Directory,"/plot_clado_", group_name, ".png", sep = "")

      #grDevices::png(cladogram_plot_path, res = 600, width = 3840, height = 2160)
      #par(mar = c(1, 0.5, 1, 1))
      #agg_clad_plot[[i]] <-plot(dendro, edge.color = Group_col, tip.color = Group_col, edge.width = 3, cex = 0.4, type = "phylogram",
      #     main = "Dendrogram of Morphometric Shape Distances")
      #agg_clad_plot[[i]]
      #legend("topleft", legend = group_val, col = palette, pch = 19, cex = 0.4)

      #dev.off()

      #if (length(dev.list()) > 0) {
      #  graphics.off()
      #}
    }
    setwd(original_wd)
  }
  Morphotree$k_medoid <- pam_clust
  Morphotree$k_medoid_plot <- k_plot
  Morphotree$dendro <- dendro
  #Morphotree$agg_clad_plot <- agg_clad_plot
  agg_fan_plot <- Filter(Negate(is.null), agg_fan_plot)
  names(agg_fan_plot) <- names(Group)

  Morphotree$agg_plot <- agg_fan_plot
  return(Morphotree)
} # End of CAV function
