% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ASCG.R
\name{ASCG}
\alias{ASCG}
\title{ASCG: Allometric Shape Changes by Groups}
\usage{
ASCG(
  GMA_object,
  Groups,
  color = NULL,
  coord_name = NULL,
  coord_subset = NULL,
  subset = FALSE,
  Mesh = NULL,
  include_whole = TRUE,
  save_plots = TRUE,
  key = NULL,
  parallel = TRUE
)
}
\arguments{
\item{GMA_object}{Either a GMA function output OR a list object
with the following structure and components:
\itemize{
\item{coords      p x k x n (landmarks, dimensions, observations) landmark array}
\item{consensus   p x k matrix of the center (mean) shape}
\item{Csize       numeric vector of allometric variable or any variable of interest}
}}

\item{Groups}{A list object containing vectors of grouping variables.}

\item{color}{Optional vector containing the color values to be used for each group.
The number of colors should be equivalent to the longest grouping vector in the
Group list.}

\item{coord_name}{an optional vector containing the names for each landmark
in the landmark configuration.}

\item{coord_subset}{An optional vector specifying the subset or group names associated
with each landmark in the landmark configuration.If provided, the resulting
interactive figure will allow toggling between these landmark subsets.}

\item{subset}{A logical value to individually conduct regression tests on each
subsetted landmark group in coord_subset. Default is NULL.}

\item{Mesh}{an optional 3D mesh file (.ply format), preferably the center shape
of the landmark array, to be thin plate spline warped into the group mean shapes.}

\item{include_whole}{a boolean value indicating whether to create visuals
for the allometric influences on shape variation for the entire sample.}

\item{save_plots}{a logical value indicating whether to save the plots to a user-specified directory.}

\item{key}{an optional character string indicating the name of a folder
either to be created or already exists in the working directory where generated
files can be stored.}

\item{parallel}{a boolean value indicating whether to use parallel processing
of the shape predictions and mesh generations for each group vector in the
Groups object. Parallel processing is strongly recommended for long Group lists.}

\item{size}{An optional vector for an allometrically influenced proxy variable.
By default, the natural log centroid size is used.}
}
\value{
A folder titled "Allometric Influences by Group" in the working
directory. This folder will contain plotly objects for each group vector in the
Group argument list. If a mesh object is provided, the function will also store
warped meshes for each group value in the grouping vector.
\itemize{
\item{Group_Shape         list object each containing a plotly plot object for each Group}
\item{Group_stats         list object each summary statistics of allometric influence on shape by group}
\item{Group_subset_stats  list object each summary statistics of allometric influence on shape by group by subsetted region}
\item{whole_Shape         list object each containing a plotly plot object of the allometric influence on all individauls irregardless of group}
}
}
\description{
This function takes a landmark configuration, a list of
grouping variables, and an allometrically influenced variable to calculate
the allometrically influenced shape changes for each group across each grouping
vector in the group list. The function uses the visualization package \strong{plotly}
to generate interactive figures and compiles each of the allometrically influenced
shapes by each grouping vector separately as an html file stored in the
working directory.
}
\details{
This function accepts either a GMA function output or a list object
with the following structure:
\itemize{
\item{coords       p x k x n (landmarks, dimensions, observations) landmark array}
\item{consensus    p x k matrix of the center (mean) shape}
\item{Csize        numeric vector of allometric variable or any variable of interest}
}

The coords object is used as the landmark configuration for all observations
in the entire sample. The consensus object is the center shape of all observations
in the entire sample. The Csize object is the allometrically influenced variable.
This function uses a linear Procrustes regression via the \emph{procD.lm} function
from the \strong{Geomorph} package. Please see the \emph{procD.lm} function for details
on how the regression is performed and the landmark data is subjected to
dimensionality reduction.

The function first performs a linear Procrustes regression of the landmark array
on the allometrically influenced variable for the entire sample if the
"include_whole" argument is TRUE. If FALSE, this step will be skipped.
Subsequently, and in parallel, each grouping vector in the Groups list will be
split from the landmark array as separate, group arrays. For each group value
in a grouping vector, these specific grouped observations will be subjected to
a linear Procrustes regression on the allometrically influenced variable. For
each group the predicted shape for the observation with the highest and smallest
values from allometrically influenced variable will be compiled into a single
interactive plotly plot html file in the specified directory folder. Each file
is automatically labeled with the phrase "Allometric Shape Changes_" followed
by the name of the unique grouping value and ended with ".html". This is
performed for each group in parallel if the argument "parallel" is TRUE.

The "Group" argument should be a list object with one or more vectors. The
function will then iterate through each grouping vector in the list and
predict the allometrically influenced shapes for each group in the grouping vector.
This function allows Grouping vectors to take on character or numeric values.
If a vector of numeric values is provided, the predicted allometric shapes
for the minimum and maximum values will be used for the \emph{plotly} comparisons.

If a 3D mesh is supplied, the landmark array should also be 3D. When a 3D
mesh object is supplied, the mesh will be warped using the thin plate spline
technique into each allometrically predicted group shape. Specifically, this
function calls the warpRefMesh function from the geomorph package to warp
mesh objects to the group mean shape.
}
\examples{
\dontrun{

# Scenario: Identify the allometrically influenced shape changes of the mandibular
# corpus by skeletal collection and mean cortical area

# Import the sample dataset Corpus_Land

data(Corpus_Land)

# Import the sample dataset Biomechanics

data(Biomechanics)

# Create a list of grouping vectors for the variables Collection and Cortical Area (CA)

Group = list(Collection = Biomechanics[1:119,]$Collection,
              'Cortical Area' = Biomechanics \%>\%
                                  group_by(ID) \%>\%
                                    summarise(CA = mean(CA))

              )

# Option 1: Perform generalized procrustes analysis on the landmark array

GMA_object = gpagen(Corpus_Land)

# Option 2: Perform the GMA analysis

GMA_object = GMA(Corpus_Land, Groups = Group)

# Conduct the SCG analysis

ASCG(GMA_object # gpagen or GMA object output
    Group # list of grouping vectors
    ) # all other arguments are optional. Colors are generated automatically.
}
}
\seealso{
\code{\link{plotly}} \code{\link{gpagen}} \code{\link{GMA}}
}
\author{
Brian Anthony Keeling
}
