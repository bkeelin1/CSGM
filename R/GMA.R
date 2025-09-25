#' @title Geometric Morphometric Analysis
#'
#' @description This functions provides a way to conduct a geometric
#' morphometric analysis and produce both informative results and graphic
#' visualizations with a single function. This function flexibly takes a wide
#' range of landmark configurations and subjects it to a Generalized Procrustes
#' Analysis (GPA). Additionally, this function can perform an analysis of bilateral
#' symmetry, permuted agglomerative or K-Medoids cluster analysis to view individual
#' relationships, and a weighted or unweighted principal component analysis. This
#' function's greatest strength is in its flexibility as well as new interactive
#' graphics to visualize shape relationships leveraging the \pkg{\link{plotly}} package.
#'
#' \itemize{\strong{GMA can perform the following:}
#'        \item{Generalized Procrustes Analysis}
#'        \item{bilateral (as)symmetry test with visualizations}
#'        \item{agglomerative and K-Medoids cluster analysis with graphics}
#'        \item{Test for overall allometric influence on shape and outlier analysis}
#'        \item{phylogenetically and/or morphometrically weighted Principal
#'              Component Analysis}
#'        \item{Principal Component plots and extreme shape visualizations}
#' }
#'
#' @param Landmarks Either a 2D array (p x k x n) (landmarks, dimensions, observations),
#' 3D array (p x k x n), n x p data.frame, n x p data matrix, p x k data.frame
#' (also tbl.df or tbl) object. Data should be the original landmark configuration
#' prior to any Procrustes Superimposition or data transformation.
#' @param k number of dimensions of the landmark data
#' @param ID An optional vector or data.frame column indicating a string of names
#' for each specimen within your dataset. This will be used as tip labels for
#' hierarchical clustering as well as in PCA shape visualizations. Default is
#' set to NULL.
#' @param Group Optional data.frame or list indicating the grouping variable(s)
#' the user wants to include for each specimen in the dataset, array, or object.
#' This variable will be used in principal component visualizations. Default is
#' set to NULL.
#' @param coord_name Optional vector or data.frame column containing the names
#' of each landmark and semilandmark in the configuration
#' @param phy A phylogenetic tree of class *phylo*. This tree will be used to
#' tailor PCA configuration and weight specimen shape variation by their
#' phylogenetic relationships. Default is set to NULL.
#' @param residual_test Logical TRUE/FALSE statement indicating whether user
#' desires to have Procrustes residuals calculated to be separately performed on
#' tests and PCA. Default is set to FALSE.
#' @param w_morpho_pca A logical TRUE/FALSE statement to indicate whether a
#' permutated hierarchical cluster analysis should be performed and subsequently
#' applied as an apriori weighting of individuals based on morphological shape
#' dissimilarity in a principal component analysis. Default is set to FALSE.
#' @param pca_centering character ("OLS", "OLS+align", "GLS", "GLS+align")
#' selecting which method of centering technique the user would prefer a principal
#' component analysis should center the Procrustes landmark (and/or Procrustes residuals)
#' configurations. The options are through Ordinary Least Squares (OLS; a common
#' approach assuming equal variances and weighting all landmarks equally), General
#' Least Squares (GLS; an approach assuming unequal variances and weighting
#' landmarks with little variance much greater than those with larger variances),
#' OLS centering with alignment of the PCA by the phylogeny or morphometric
#' cluster, or GLS centering with alignment of the PCA by the phylogeny or
#' morphometric cluster. When a phylogeny is absent and w_morpho_pca is set to
#' default (= NULL), alignment or GLS centering makes little sense and subsequently
#' will result in an error. Default option is set to "OLS".
#' @param bilat_sym A logical TRUE/FALSE statement indicating whether or not the
#' user desires an analysis of bilateral symmetry to conduct tests to evaluate
#' the degree which shape variation is influenced by factors of directional and/or
#' fluctuating asymmetry. If this argument is set to TRUE, user must add additional
#' parameters required by the \code{\link{bilat.symmetry}} function from the
#' \pkg{geomorph} package. Particularly relevant is the arguments "side", "replicate",
#' and "object.sym". See \code{\link{bilat.symmetry}} for more details.
#' @param subset An optional group vector which indicates if there are
#' subdivisions in the landmark configuration. This will be used to generate
#' extreme shape graphics using plotly to toggle different landmark configuration
#' groups. This vector should only name the landmark number by their respective
#' group and not of their specific axes coordinates (x,y,z).
#' @param curves An optional curve matrix via the \code{gen_curves} function to
#' ensure equidistant semilandmark points.
#' @param side An optional vector which EITHER indicates which individuals
#' represent the left and right sides OR indicates which landmarks in the
#' configuration represents the left and right sides. Only necessary when
#' argument 'bilat_sym = TRUE' and the type of symmetry is set to "matching".
#' Default is set to NULL. Please see the function \code{bilat.symmetry} from
#' the \pkg{geomorph} package for more information.
#' @param replicate An optional vector which indicates which specimens in the
#' sample are replicates of each other. Only necessary when argument
#' 'bilat_sym = TRUE' AND the type of symmetry is set to "matching". Default is
#' set to NULL. Please see the function \code{bilat.symmetry} from the \pkg{geomorph}
#' package for more information.
#' @param key Character string. Identifier for output files and directories. Default is NULL.
#' @param parallel Logical value indicating whether to use parallel processing
#' @param print_progress An logical value indicating whether the function should
#' generate a progress bar in the Console with updates on the function process.
#' @param ... Any additional arguments from the \code{bilat.symmetry}, \code{gpagen},
#' \code{hullgen}, or \code{CAV} functions from \pkg{CSGM} and \pkg{geomorph}
#' can be passed through here.
#'
#' @returns A list of class "GMA" containing:
#' \itemize{
#'   \item{GPA_data:    Generalized Procrustes Analysis results:
#'     \itemize{
#'       \item{coords:    Procrustes superimposed coordinates}
#'       \item{Csize:   Centroid sizes}
#'       \item{consensus:   Mean shape configuration}
#'       \item{Outliers:    Results of outlier analysis}
#'     }
#'   }
#'   \item{sym_data:  If bilat_sym=TRUE, bilateral symmetry analysis results
#'     \itemize{
#'       \item{symm.shape:    Symmetric component of shape variation}
#'       \item{asymm.shape:   Asymmetric component of shape variation}
#'       \item{Mantel_tests:    Results of symmetry significance tests}
#'     }
#'   }
#'   \item{morphotree:    Morphological clustering results via CAV:
#'     \itemize{
#'       \item{resid:   Clustering results and plots for Procrustes data}
#'       \item{asym:    Optional clustering and plots of asymmetric component}
#'       \item{sym:   Optional clustering and plots of symmetric component}
#'     }
#'   }
#'   \item{PCA:   Principal Component Analysis results
#'     \itemize{
#'       \item{PCA_var:   Principal component summary statistics via \code{pca_variances}}
#'       \item{PCs:   Identified significant principal components}
#'       \item{scores:    PC scores for each specimen}
#'       \item{loadings:    PC loadings for shape variables}
#'       \item{PCA_plots:   Interactive \pkg{plotly} Principal Component plots for each meaningful PC comparison}
#'       \item{PCA_groups:    Optional when Group != NULL. Interactive \pkg{plotly} Principal Component plots by Group}
#'       \item{lollipop_PCs:    Visualize shape variation for each meaningful Principal Component}
#'     }
#'   }
#'   \item{PCA_sym:   When bilat.sym = TRUE. Principal Component Analysis results for symmetric shape variation
#'     \itemize{
#'       \item{PCA_var:   Principal component summary statistics via \code{pca_variances}}
#'       \item{PCs:   Identified significant principal components}
#'       \item{scores:    PC scores for each specimen}
#'       \item{loadings:    PC loadings for shape variables}
#'       \item{PCA_plots_sym:   Interactive \pkg{plotly} Principal Component plots for each meaningful PC comparison}
#'       \item{PCA_groups_sym:    Optional when Group != NULL. Interactive \pkg{plotly} Principal Component plots by Group}
#'       \item{lollipop_PCs_sym:    Visualize shape variation for each meaningful Principal Component}
#'     }
#'   }
#'   \item{PCA_asym:    When bilat.sym = TRUE. Principal Component Analysis results for asymmetric shape variation:
#'     \itemize{
#'       \item{PCA_var:   Principal component summary statistics via \code{pca_variances}}
#'       \item{PCs:   Identified significant principal components}
#'       \item{scores:    PC scores for each specimen}
#'       \item{loadings:    PC loadings for shape variables}
#'       \item{PCA_plots_asym:    Interactive \pkg{plotly} Principal Component plots for each meaningful PC comparison}
#'       \item{PCA_groups_asym:   Optional when Group != NULL. Interactive \pkg{plotly} Principal Component plots by Group}
#'       \item{lollipop_PCs_asym:   Visualize shape variation for each meaningful Principal Component}
#'     }
#'   }
#' }
#'
#' @details
#' GMA takes a landmark configuration and conducts a thorough geometric morphometric
#' analysis which includes Generalized Procrustes Analysis (GPA), bilateral symmetry
#' analysis, an interactive Principal Component Analysis (PCA), and an optional
#' unsupervised K-Medoids and Agglomerative Hierarchcical cluster analysis. This
#' function was designed to allow a user to quick and efficient way to conduct an
#' exploratory data analysis with shape data and provides interactive graphics
#' useful for subsequent down-stream analyses.
#'
#' This function accepts a few standard input data formats for the argument Landmarks,
#' including a "p x k x n" 2D or 3D array, an "n x p" data frame or matrix, a
#' "p x k" data frame or matrix, as well as table objects of the same dimensions.
#' Data should be the original landmark configuration prior to any data transformation.
#' Landmark configurations are subjected to Generalized Procrustes Analysis (GPA)
#' using the \code{gpagen} from the \pkg{geomorph} package. The argument (...)
#' allows the user to add in additional arguments which can be passed to \code{gpagen}.
#' When semilandmarks are present and a curves matrix are provided, the function
#' implements sliding semilandmark procedures to minimize bending energy and provide
#' a more improved overall fit.
#'
#' For bilateral structures, the function implements comprehensive symmetry analysis through
#' either matching symmetry (comparing separate structures) or object symmetry (analyzing
#' symmetric structures). This analysis decomposes shape variation into symmetric and
#' asymmetric components, enabling investigation of both directional and fluctuating
#' asymmetry. This function takes the approach following Lazic et al. (2008) and
#' additionally provides mantel significance tests to evaluate whether the same
#' factors are influencing the symmetric and asymmetric shape variation on the
#' left and right sides.
#'
#' Morphological clustering is performed through \code{CAV}, which implements both supervised
#' agglomerative hierarchical clustering and unsupervised k-medoids clustering. The
#' agglomerative clustering method is automatically selected based on the highest
#' agglomerative coefficient among available methods (average, single, complete, ward,
#' weighted). For k-medoids clustering, the optimal number of clusters is determined
#' through silhouette width analysis (see \code{\link{PMA}} and \code{\link{pvclust}}
#' for more details on how each analysis is performed). Landmarks are transformed
#' to Mahalanobis distances with parallel processing. Additional arguments can be passed
#' to CAV. See \code{\link{CAV}} for more details on argument names and choices.
#'
#' Principal Component Analysis is conducted using either standard or phylogenetically/
#' morphometrically weighted approaches. The function determines significant PCs through
#' a principal component based factor analysis via \code{fa.parallel}. For visualizations,
#' the function generates several interactive plots using the workhorse package
#' \pkg{\link{plotly}}.
#'
#' \itemize{Generated Plots include:
#'    \item{2D and 3D Principal Component plots for each combination of meaningful components}
#'    \item{When Group != NULL, 2D and 3D Princiapl Component plots with convex hulls overlaid}
#'    \item{2D or 3D graphics showing Principal Component shape variation for each meaningful component}
#' }
#'
#' When analyzing bilateral symmetry, both the symmetric and asymmetric shape components undergo
#' separate clustering and PCA analyses. This enables detailed investigation of patterns
#' in both the shared (symmetric) and unique (asymmetric) aspects of shape variation.
#' The function generates comprehensive visualization plots for both components.
#'
#' For shape data with semilandmarks, the function automatically handles sliding while
#' maintaining correct geometric relationships. Curve definitions must be provided to
#' enable proper sliding along curves as argument (curves). The sliding process
#' minimizes either bending energy or Procrustes distances based on the selected criterion.
#'
#' Cluster analysis results are organized in the "GMA Results" folder directory,
#' Interactive visualization plots are stored in function output but can saved to
#' the directory by using the \code{saveWidget} function. GMA supports parallel processing
#' to optimize computational efficiency for large datasets as well as prints function
#' progress.
#'
#' @author Keeling et al., 2025
#'
#' @seealso \pkg{\link{cluster}} \pkg{\link{geomorph}} \pkg{\link{plotly}} \pkg{\link{psych}} \pkg{\link{pvclust}}
#'
#' @importFrom geomorph gpagen plotOutliers procD.lm arrayspecs bilat.symmetry two.d.array
#' @importFrom progress progress_bar
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Example 1: Exploration of mandibular corpus shape variation
#' #            in a human study sample
#'
#' # First load in example landmark data
#' data(Corpus_Land)
#'
#' # Create a semilandmark curves matrix for equidistant spacing between
#' # landmarks during Procrustes superimposition
#'
#'
#' # load in helper data to generate a curves matrix
#'
#' data(anchor_indices)
#' data(semiland_indices)
#'
#' #Use the gen_curves function
#'
#' curves = gen_curves(10, # There are 10 total curves (5 external; 5 internal)
#'                     3, # There are six anchor points (lingual x2, basal x2, buccal x2) for each curve
#'                     590, # There are 590 total landmarks
#'                     def_points = anchor_indices
#'                     def_semi = semiland_indices
#'                    )
#'
#'
#' # Create grouping variables wih example data Biomechanics
#'
#' data(Corpus_Prop)
#' ID = dimnames(Corpus_Land)[[3]]
#' Collection = Biomechanics[1:50,]$Collection
#'
#' # Conduct the GMA test
#'
#' output <- GMA(Corpus_Land,  # Complete landmark configuration
#'               k = 3,  # 3D data
#'               ID,  # Specimen identifiers
#'               Group = Collection,  # Grouping variable
#'               coord_name = paste("Landmark", 1:590, sep="_"),  # Landmark names
#'               subset = paste(rep(c("LM1M2", "LP3P4", "Symphysis", "RP3P4", "RM1M2"), each = 118)),
#'               bilat_sym = FALSE,  # No bilat symmetry
#'               key = "Mandibular_Corpus",
#'               curves = curves # add in curves matrix
#' ) # keep the rest defaults
#'
#'
#' # Statistics ________________________________________________________________
#'
#' # Check top outliers
#' GMA$GPA_stats$outliers
#'
#' # Check percentage of variation in shape explained by log centroid size (allometry)
#' GMA$GPA_stats$Csize
#'
#' # Agglomerative Cluster Analysis ____________________________________________
#'
#' #view shape relationships through agglomerative clustering
#' # NOTE K-Medoids plot and cluster plot is in file directory
#'  plot(results_standard$morphotree$cluster$dendro)
#'
#' # PCA Summary _______________________________________________________________
#'
#' output$PCA$PCA_var$PCA_var
#'
#' # PCA Plots _________________________________________________________________
#' # View 2D Principal Component plot for PCs 1-2
#' output$PCA_plots$PCs_1_2
#'
#' # View *Group* 2D Principal Component plot for PCs 1-2
#'
#' output$PCA_groups$Collection$PCs_1_2
#'
#'
#' # View 3D Principal Component plot for PCs 1-3
#' output$PCA_plots$PCs_1_2
#'
#' # View *Group* 3D Principal Component plot for PCs 1-3
#' output$PCA_groups$Collection$PCs_1_2_3
#'
#'
#' # PCA Shape Variation _______________________________________________________
#'
#' # View shape variation for PC 1 via \pkg{plotly}
#'
#' output$PCA$lollipops_PCs[[1]]
#'
#'
#'
#' # Example 2: Bilateral symmetry analysis of the left and right mandibular corpus
#'
#' # First load in the original mandibular corpus landmark configuration
#' data(Corpus_Land)
#'
#' # Create a semilandmark curves matrix for equidistant spacing between landmarks
#' # during Procrustes superimposition (important for semilandmarks).
#'
#' #Use the gen_curves function
#'
#' curves_sym = gen_curves(8, # There are 8 total curves (4 external; 4 internal)
#'                     3, # There are three fixed points (lingual, basal, buccal) for each curve
#'                     480 # There are 600 total landmarks
#'                     def_points = c(1, 31, 60, 61, 91, 120,       #LM1M2
#'                                    121, 151, 180, 181, 211, 240, #LP3P4
#'                                    241, 271, 300, 301, 331, 360, #RP3P4
#'                                    361, 391, 420, 421, 451, 480) #RM1M2
#'                     semi_points = c(2:30,32:59, 62:90, 92:119,
#'                                     122:150, 152:179, 182:210, 212:239,
#'                                     242:270, 272:299, 302:330, 332:359,
#'                                     362:390, 392:419, 422:450, 452:479)
#'                    )
#'
#'  # Create grouping variables
#'
#' ID = dimnames(Corpus_Land)[[3]]
#' Collection = data(Biomechanics)$Collection
#'
#' # Create a data table identifying the right and left sides for the analysis
#'
#' Paired_sym = cbind(left = c(1:120, 121:240),
#'                   right = c(361:480, 241:360))
#'
#' # Conduct the GMA test
#'
#' output <- GMA(Corpus_Land,  # Complete landmark configuration
#'               k = 3,  # 3D data
#'               ID,  # Specimen identifiers
#'               Group = Collection,  # Grouping variable
#'               coord_name = paste("LM", 1:600, sep=""),  # Landmark names
#'               subset = paste(rep(c("LM1M2", "LP3P4", "Symphysis", "RP3P4", "RM1M2"), each = 120)),
#'               bilat_sym = TRUE,  # No bilat symmetry
#'               curves = curves_sym, # curves matrix for semilandmarks
#'               side = Paired_sym,  # identifying the landmarks associated with each side of the mandibular corpus
#'               sym_type = "object", # the type of bilateral symmetry analysis
#'               key = "Mandibular_Corpus" # defines the folder directory under "GMA Results - Mandibular_Corpus"
#'              ) # keep the rest defaults
#'
#' # Statistics ________________________________________________________________
#'
#' # View Bilateral Symmetry Results
#'
#'  GMA$sym_data
#'
#' # View Mantel associations between the left and right fluctuating asymmetry
#'
#'  GMA$sym_data$Mantel_tests
#'
#' # Check top outliers
#' GMA$GPA_stats$outliers
#'
#' # Check percentage of variation in shape explained by log centroid size (allometry)
#' GMA$GPA_stats$Csize
#'
#' # Agglomerative Cluster Analysis ____________________________________________
#'
#' #view shape relationships through agglomerative clustering
#' # NOTE K-Medoids plot and cluster plot is in file directory
#'  plot(output$morphotree$cluster$dendro)
#'
#'  # View symmetric shape component relationships
#'  plot(output$morphotree$cluster_sym$dendro)
#'
#'  # View asymmetric shape component relationships
#'  plot(output$morphotree$cluster_asym$dendro)
#'
#'
#' # PCA Summary _______________________________________________________________
#'
#' output$PCA$PCA_var$PCA_var
#'
#' output$PCA$PCA_var$PCs # number of meaningful PCs detected
#'
#' # PCA Plots _________________________________________________________________
#' # View 2D Principal Component plot for PCs 1-2
#' output$PCA$PCA_plots$PCs_1_2
#'
#' # View *Group* 2D Principal Component plot for PCs 1-2
#' output$PCA$PCA_groups$Collection$PCs_1_2
#'
#' #View for symmetric shape components
#' output$PCA_sym$PCA_groups_sym$Collection$PCs_1_2
#'
#' #View for asymmetric shape components
#' output$PCA_asym$PCA_groups_asym$Collection$PCs_1_2
#'
#'
#'
#'
#' # View 3D Principal Component plot for PCs 1-3
#' output$PCA$PCA_plots$PCs_1_2
#'
#' # View *Group* 3D Principal Component plot for PCs 1-3
#' output$PCA_groups$Collection$PCs_1_2_3
#'
#' #View for symmetric shape components
#' output$PCA_sym$PCA_groups_sym$Collection$PCs_1_2_3
#'
#' #View for asymmetric shape components
#' output$PCA_asym$PCA_groups_asym$Collection$PCs_1_2_3
#'
#'
#' # PCA Shape Variation _______________________________________________________
#'
#' # View shape variation for PC 1 via \pkg{plotly}
#' output$PCA$lollipops_PCs[[1]]
#'
#' # View symmetric shape variation for PC 1 via \pkg{plotly}
#' output$PCA_sym$lollipop_PCs_sym[[1]]
#'
#' # View asymmetric shape variation for PC 1 via \pkg{plotly}
#' output$PCA_asym$lollipop_PCs_asym[[1]]
#' }

GMA <- function(Landmarks,
                k,
                ID = NULL,
                Group = NULL,
                coord_name = NULL,
                subset = NULL,
                phy = NULL,
                bilat_sym = FALSE,
                w_morpho_pca = TRUE,
                pca_centering = "OLS",
                key = NULL,
                print_progress = TRUE,
                ...) {

  if (isTRUE(print_progress)) {
    pb <- progress_bar$new(
      format = "Processing :what [:bar] :percent ETA: :eta",
      clear = FALSE, total = 100, width = 60)
  }

  Directory = getwd()

  if(!is.null(key)) {
    if(!dir.exists(paste("GMA Results", "-", key))) {
      dir.create(paste("GMA Results", "-", key))
    }
    ct.dir = paste("GMA Results", "-", key)
  } else {
    if(!dir.exists("GMA Results")) {
      dir.create("GMA Results")
    }
    ct.dir = "GMA Results"
  }
  path = paste(Directory, "/", ct.dir, sep = "")
  setwd(path)

  # Data wrangling for Landmarks

  GMA_Output <- list() # output data

  args = list(...)

  gpagen_params = names(formals(geomorph::gpagen))

  gpa_parameters = list(A = NULL,
                        curves = NULL,
                        approxBE = FALSE,
                        verbose = FALSE,
                        print.progress = FALSE
                      ) %>% modifyList(., args) %>%
                            .[intersect(names(.), gpagen_params)]

  if (isTRUE(print_progress)) {
    pb$tick(5, tokens = list(what = "GPA generation"))
  }

  if (class(Landmarks) == "array") {
    land_origin = Landmarks
    gpa_parameters = modifyList(gpa_parameters, list(A = Landmarks))
    GPA_data = do.call(gpagen, gpa_parameters)
    Data = data.frame(two.d.array(GPA_data$coords))
    GMA_Output$GPA <- GPA_data
  }

  if (any(class(Landmarks) %in% c("data.frame", "tbl_df", "tbl"))) {

    if(ncol(Landmarks) == 3) {
      if(is.null(ID)){
        stop("Error: Cannot identify how many individuals are in the sample
                       to calculate number of landmarks per individual. Please add vector to argument 'ID'.")
      } else{
        n <- c(rep(1:length(ID), each = (nrow(Landmarks)/n)))
        p <- c(rep(1:(nrow(Landmarks)/length(ID)), times = length(ID)))
        matrix1 <- matrix(Landmarks[,1], nrow = length(ID), ncol = length(Landmarks)/length(ID))
        matrix2 <- matrix(Landmarks[,2], nrow = length(ID), ncol = length(Landmarks)/length(ID))
        matrix3 <- matrix(Landmarks[,3], nrow = length(ID), ncol = length(Landmarks)/length(ID))
        array_3d <- array(dim = c(length(ID), k, nrow(Landmarks)/length(ID)))


        array_3d[,1,] <- matrix1
        array_3d[,2,] <- matrix2
        array_3d[,3,] <- matrix3
      }
    } else if(ncol(Landmarks) == 2) {
      if(is.null(ID)){
        stop("Error: Cannot identify how many individuals are in the sample
                       to calculate number of landmarks per individual. Please add vector to argument 'ID'.")
      } else{
        n <- c(rep(1:length(ID), each = length(ncol(Landmarks))))
        p <- c(rep(1:length(ncol(Landmarks)), times = length(ID)))
        matrix1 <- matrix(Landmarks[,1], nrow = length(ID), ncol = length(Landmarks),
                          dimnames = list(unique(ID), unique(p)))
        matrix2 <- matrix(Landmarks[,2], nrow = length(ID), ncol = length(Landmarks),
                          dimnames = list(unique(ID), unique(p)))
        array_2d <- array(dim = c(length(ID), k, nrow(Landmarks)/length(ID)))


        array_3d[,1,] <- matrix1
        array_3d[,2,] <- matrix2
      }
    }
    Landmarks = arrayspecs(Landmarks, length(Landmarks)/k, k)
    land_origin = Landmarks
    gpa_parameters = modifyList(gpa_parameters, list(A = Landmarks))
    GPA_data = do.call(gpagen, gpa_parameters)
    Data = data.frame(two.d.array(GPA_data$coords))
    GMA_Output$GPA_data <- GPA_data
  } # End of data frame transformation

  if(is.null(ID)) {
    if(k == 3) {
      ID <- c(rep(1:dim(GPA_data$coords)[3]))
    } else {
      ID <- c(rep(1:dim(GPA_data$coords)[2]))
    }
  }

  if(is.data.frame(ID)) {
    ID <- ID %>% pull(1)
  }

  if (!is.character(ID)) {
    ID <- as.character(ID)
  }

  if(is.null(Group)) {
    Group = c(rep("NA", times = dim(GPA_data$coords)[3]))
    Group = data.frame(Group = Group)
  }

  if(isTRUE(print_progress)){
    pb$tick(15, tokens = list(what = "Symmetry analysis"))
  }
  # Bilateral symmetry analysis

  if (bilat_sym == TRUE) { #To designate whether an analysis of symmetry is desired

      bilat_params = names(formals(geomorph::bilat.symmetry))

      bilat_parameters = args %>%
                           modifyList(., list(A = GPA_data,
                                              print.progress = FALSE)) %>%
                              .[intersect(names(.), bilat_params)]
      sym_data = do.call(bilat.symmetry,bilat_parameters)
      # Mantel test calculation
      side = if(!is.null(bilat_parameters$side)) bilat_parameters$side else bilat_parameters$land.pairs
      Mantel_sym = mantel_sym(sym_data, side, perm = 1000)
      sym_data$Mantel_tests <- Mantel_sym
      # Dataset generation for PCA
      Data_asym = data.frame(two.d.array(sym_data$asymm.shape))
      GPA_data_asym = sym_data$asymm.shape
      Data_sym = data.frame(two.d.array(sym_data$symm.shape))
      GPA_data_sym = sym_data$symm.shape

      # Procrustes Residuals
      residuals = Evomorph::GpaResiduals(land_origin, GPA_data$coords)
      residuals = geomorph::arrayspecs(data.frame(residuals$resid), ncol(data.frame(residuals$resid))/k, k)
      sym_resid_data = modifyList(bilat_parameters, list(A = GPA_data)) %>% do.call(bilat.symmetry,.)

      sym_resid_data$Csize <- GPA_data$Csize
      sym_data$Csize <- GPA_data$Csize
      sym_data$residuals <- residuals

      resid_asym <- sym_resid_data$asymm.shape
      resid_sym <- sym_resid_data$symm.shape
      GMA_Output$sym_resid_data <- sym_resid_data

      GMA_Output$sym_data <- sym_data
  }

  if(isTRUE(print_progress)){
    pb$tick(20, tokens = list(what = "Create Morpho Clusts"))
  }

  morphotree = list()
  GMA_Output$morphotree <- morphotree

  CAV_params = names(formals(CAV))
  CAV_parameters = modifyList(args, list(Data = Data, Group = Group, ID = ID,
                                         nboot = 1000, Dir_name = "Residuals")) %>%
                        .[intersect(names(.), CAV_params)]

  Pred_ncomp = PCA_variances(Data_PCA = geomorph::gm.prcomp(A = GPA_data$coords))$PCs


  dimnames(GPA_data$coords)[[3]] <- ID
  resid_data = Evomorph::GpaResiduals(land_origin, GPA_data$coords)
  GMA_Output$residuals = arrayspecs(resid_data$resid, dim(GPA_data$coords)[1], k)

  CAV_parameters = modifyList(CAV_parameters, list(Data = data.frame(resid_data$resid), Dir_name = "Procrustes Residuals", Pred_ncomp = Pred_ncomp))

  cluster = tryCatch({do.call(CAV, CAV_parameters)},
             error = function (e) {
               tryCatch({
                 CAV_parameters = modifyList(CAV_parameters, list(Data = data.frame(geomorph::two.d.array(GPA_data$coords)),
                                                                  Dir_name = "Procrustes Residuals",
                                                                  Pred_ncomp = Pred_ncomp, method = "euclidean")
                                             )
                 do.call(CAV,CAV_parameters)
               }, error = function (e) {
                print("CAV Function Failed. Cannot perform cluster analysis with data that has a standard deviation of 0. Select different method")
               })
            })


  dendro_resid = cluster$dendro
  GMA_Output$morphotree$resid <- cluster

  if (bilat_sym == TRUE) { # Hierarchical cluster generation when bilateral symmetry analysis is desired

    CAV_parameters = modifyList(CAV_parameters, list(Data = data.frame(two.d.array(resid_asym)), Dir_name = "Asymmetry"))
    cluster = do.call(CAV, CAV_parameters)
    dendro_asym = cluster$dendro
    GMA_Output$morphotree$asym <- cluster

    CAV_parameters = modifyList(CAV_parameters, list(Data = data.frame(two.d.array(resid_sym)), Dir_name = "Symmetry"))
    cluster = do.call(CAV, CAV_parameters)
    dendro_sym = cluster$dendro
    GMA_Output$morphotree$sym <- cluster
  }

  if(isTRUE(print_progress)){
    pb$tick(20, tokens = list(what = "Compute PCAs"))
  }

  if (bilat_sym == TRUE) {
    # Calculation of the Principal Component Analysis for the Symmetric and Asymmetric Shape Components
    dimnames(GPA_data_asym)[[3]] <- ID
    dimnames(GPA_data_sym)[[3]] <- ID

    align.to.phy = if(grepl("align", pca_centering)) TRUE else FALSE
    GLS = if(grepl("GLS", pca_centering)) TRUE else FALSE
    transform = if(grepl("align", pca_centering)) TRUE else FALSE

    if (!is.null(phy) | isTRUE(w_morpho_pca)) {
      if(!is.null(phy)) {
        phy_asym = phy
        phy_sym = phy
      } else {
        phy = dendro_resid
        phy_asym = dendro_asym
        phy_sym = dendro_sym
      }
    } else {
      phy = NULL
      phy_asym = NULL
      phy_sym = NULL
    }

    pca_params = names(formals(geomorph::gm.prcomp))
    pca_parameters = modifyList(gpa_parameters, list(A = GPA_data$coords,
                                                     phy = phy,
                                                     align.to.phy = align.to.phy,
                                                     GLS = GLS, transform = transform)) %>%
                     .[intersect(names(.), pca_params)]
    Data_PCA = do.call(gm.prcomp, pca_parameters)
    pca_parameters = modifyList(pca_parameters, list(A = GPA_data_asym, phy = phy_asym))
    Data_PCA_asym = do.call(gm.prcomp, pca_parameters)
    pca_parameters = modifyList(pca_parameters, list(A = GPA_data_sym, phy = phy_sym))
    Data_PCA_sym = do.call(gm.prcomp, pca_parameters)

  } else {

    align.to.phy = if(grepl("align", pca_centering)) TRUE else FALSE
    GLS = if(grepl("GLS", pca_centering)) TRUE else FALSE
    transform = if(grepl("align", pca_centering)) TRUE else FALSE

    phy = if(isTRUE(w_morpho_pca) & is.null(phy)) dendro_resid

    pca_params = names(formals(geomorph::gm.prcomp))
    pca_parameters = modifyList(gpa_parameters, list(A = GPA_data$coords,
                                                     phy = phy,
                                                     align.to.phy = align.to.phy,
                                                     GLS = GLS, transform = transform)) %>%
                     .[intersect(names(.), pca_params)]
    Data_PCA = do.call(gm.prcomp, pca_parameters)

  } # End of non-bilateral symmetric analysis

  if(isTRUE(print_progress)){
    pb$tick(10, tokens = list(what = "Analyzing PCA Data"))
  }

  #PCA results of original component
  PCA_variance = PCA_variances(Data_PCA = Data_PCA)
  PCA_var = PCA_variance$PCA_var
  PCs = PCA_variance$PCs
  GMA_Output$PCA <- PCA_variance

  if(bilat_sym == TRUE) { # If a PCA of bilateral symmetry was computed to summarize the a/symmetric variances
    # PCA results of asymmetric component
    PCA_variance = PCA_variances(Data_PCA_asym, PCs = FALSE)
    PCA_var_asym = PCA_variance$PCA_var
    PCs_asym = PCs
    GMA_Output$PCA_asym <- PCA_variance

    # PCA results of symmetric component
    PCA_variance = PCA_variances(Data_PCA_sym, PCs = FALSE)
    PCA_var_sym = PCA_variance$PCA_var
    PCs_sym = PCs
    GMA_Output$PCA_sym <- Data_PCA_sym
  }

  if(isTRUE(print_progress)){
    pb$tick(20, tokens = list(what = "Creating PCA Visuals"))
  }

  lollipops <- lolligen(Data_PCA = Data_PCA, PCs = PCs, ID = ID, Group = Group, coord_name = coord_name, subset = subset, PCA_var = PCA_var, k = k)
  GMA_Output$PCA$PCA_plots <- lollipops$PCA_list
  GMA_Output$PCA$PCA_groups <- lollipops$PCA_group
  GMA_Output$PCA$lollipop_PCs <- lollipops$lollipop_PCs

  # Apply lollipop graphing function by whether bilateral symmetry is warranted or not
  if (bilat_sym == TRUE) {
    lollipops_asym <- lolligen(Data_PCA_asym, PCs = PCs_asym, ID = ID, Group = Group, coord_name = coord_name, subset = subset, PCA_var = PCA_var_asym, k = k)
    GMA_Output$PCA_asym$asym_PCA_plots <- lollipops_asym$PCA_list
    GMA_Output$PCA_asym$asym_PCA_groups <- lollipops_asym$PCA_group
    GMA_Output$PCA_asym$asym_lollipop_PCs <- lollipops_asym$lollipop_PCs

    lollipops_sym <- lolligen(Data_PCA_sym, PCs = PCs_sym, ID = ID, Group = Group, coord_name = coord_name, subset = subset, PCA_var = PCA_var_sym, k = k)
    GMA_Output$PCA_sym$sym_PCA_plots <- lollipops_sym$PCA_list
    GMA_Output$PCA_sym$sym_PCA_groups <- lollipops_sym$PCA_group
    GMA_Output$PCA_sym$sym_lollipop_PCs <- lollipops_sym$lollipop_PCs
  }

  if(isTRUE(print_progress)){
    pb$tick(9, tokens = list(what = "Making Final Adjustments"))
  }

  calregs <- function(GPA_data) {
    output <- list()
    tryCatch({
      GPA_data <<- GPA_data
      size_formula <- as.formula("coords ~ log(Csize)")
      output$Csize <- procD.lm(size_formula, data = GPA_data)
      output$outliers <- unique(plotOutliers(GPA_data$coords, inspect.outliers = FALSE))
      output$dimensions <- dim(GPA_data$coords)
      }, error = function (e) {print("Centroid Size is the same in all individuals after Procrustes Registration. Please check the nature of the data.")}
    )
    return(output)
  }

  GMA_Output$GPA_stats <- calregs(GPA_data)

  if(isTRUE(print_progress)){
    pb$tick(1, tokens = list(what = "FINISHED!"))
  }

  class(GMA_Output) = "GMA"
  return(GMA_Output)

} # end GMA function
