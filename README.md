# *CSGM* - A Toolkit to Conduct Cross-Sectional Geometry(ics) Morphometric Analyses

## Description

Cross-sectional Geometry and Morphometrics (CSGM) analysis consists of two popular interdisciplinary methods in biological anthropology to study the relationships between shape and biomechanical function in contexts such as developmental biology, behavioral plasticity, and evolutionary morphology. As such, the integration of these two research methods involve the comparisons of complex, multidimensional shape data from multiple skeletal regions to multivariate, and often multicollinear, biomechanical variables. Thus, we developed the CSGM package in R to provide a statistically holistic framework and rigorous nested hypothesis model testing approach to facilitate a complete statistical data analysis. 

The CSGM package offers an innovative and comprehensive toolkit to automate an entire geometric morphometric and inferential statistical analysis from start to finish. This package addresses the data complexity issues within CSGM studies which often require repetitive analyses and inter/intra-bone analyses by providing a holistic statistical approach to functionally test hypotheses and analyze complex multivariate data with only a few functions. This multimodal approach combines geometric morphometric analysis, correlation, covariation, and regression analyses to offer novel insights into the relationships between complex geometric and biological data. By applying a hypothesis model approach, the CSGM package allows researchers to organize multiple hypothesis testing across multivariate dataset comparisons. The functions within this package also offer accessible and interactive visualizations to easily interpret data relationships. Therefore, this package offers numerous ways to accessibly study geometric morphometrics and interpret complex data relationships needed to answer complex questions in human biology.

# Installation Instructions:

------------------------------------------------------------------------

## Install the `devtools` package to install Github package repositories like CSGM

```{r}
install.packages("devtools") # to install GitHub packages
```

## Now Install the `CSGM` package from GitHub!

```{r}
devtools::install_github("bkeelin1/CSGM",
                         build_vignettes = TRUE, 
                         dependencies = TRUE) 
```

## Follow along in RStudio or Posit Cloud and use this README through the vignette CSGM-Workflow

```{r}
library(CSGM)
# CSGM function to extract vignette in file directory
get_vignette("CSGM-Workflow")
file.edit("CSGM-Workflow.Rmd")
```

Note: When the Vignette opens, please look in the Script (top left) Panel and click the button "Visual" to view this vignette in a more visually appealing mode.

------------------------------------------------------------------------

# Visualizing the CSGM Workflow

![](images/clipboard-1066059136.png)

# Step 1: Importing Example Data

### Load the landmark array, and helper indices

```{r Load_Landmarks,echo=TRUE,eval=FALSE}
data("Corpus_Land")
data("anchor_indices")
data("semiland_indices")
```

## Load in Corpus Cross-sectional Properties

```{r, Load_Properties,echo=TRUE,eval=FALSE}
data("Corpus_Prop")
```

## Create a curves matrix of the semilandmark array

```{r gencurves, echo=TRUE,eval=FALSE}
curves = curves = gen_curves(10, # There are 10 total curves (5 external; 5 internal)
                    6, # There are six anchor points (lingual x2, basal x2, buccal x2) for each curve
                    590, # There are 590 total landmarks
                    def_points = anchor_indices, # Landmark numbers which are anchor points
                    def_semi = semiland_indices, # Landmark numbers which are semilandmark points
                   )
```

# Step 2: Exploratory Data Analysis

## Conduct a Geometric Morphometric Analysis (GMA)

```{r GMA,echo=TRUE,eval=FALSE}
# Let's make ID and Group arguments for the GMA analysis

ID = paste(1:50)

Group = list(Collection = Corpus_Prop[1:50,]$Collection) # Skeletal Collection

# Conduct the analysis

output <- GMA(Corpus_Land,  # Complete landmark configuration
              k = 3,  # 3D data
              ID,  # Specimen identifiers
              Group = Group,  # Grouping variable
              coord_name = dimnames(Corpus_Land)[[1]],  # Landmark names
              subset = paste(rep(c("LM1M2", 
                                   "LP3P4", 
                                   "Symphysis", 
                                   "RP3P4", 
                                   "RM1M2"), 
                                 each = 118)),
              phy = NULL, # We do not have a phylogeny for this data
              bilat_sym = TRUE,  # We will also conduct a bilateral symmmetry test
              land.pairs = cbind(left = c(1:236), right = c(355:590)), # separate the left and right cross-sections
              object.sym = TRUE, # needed for object symmetry (we are looking at a mandible)
              ind = ID, # needed for bilateral symmetry
              w_morpho_pca = TRUE, # We will weight our PCA analysis by individual Mahalanobis distances
              key = "Mandibular_Corpus",
              curves = curves # add in curves matrix
              ) 


# NOTE: Procrustes Outliers were individually checked and approved as representing real shape variation

```

## Procrustes Shape Variation

### View morphoplots on Procrustes Residuals (in Folder as well)

```{r Morphotrees,echo=TRUE,eval=FALSE}

# View the K-Medoids plot
output$morphotree$resid$k_medoid_plot

#view the agglomerative cluster plot
do.call(ape::plot.phylo, output$morphotree$resid$agg_plot$Collection)
```

### View Weighted PCA plots of Procrustes Landmarks

```{r,echo=TRUE,eval=FALSE}

# View 3D PCA without groups 
output$PCA$PCA_plots$PC_1_2_3  
# View 3D PCA by skeletal collection 
output$PCA$PCA_groups$Collection$PC_1_2_3
```

Notes:

1.  Large degree of overlap in skeletal collection supported by the cluster analysis.

### View Weighted PCA shape extremes of Procrustes Landmarks

```{r Extreme_Shapes,echo=TRUE,eval=FALSE}
# PC 1 shape extremes
output$PCA$lollipop_PCs[[1]]

# PC 2 shape extremes
output$PCA$lollipop_PCs[[2]]

# PC 3 shape extremes
output$PCA$lollipop_PCs[[2]]
```

Notes:

1.  Positive PC1 shapes represent mandibles with generally smaller corpus and symphysis heights. The symphysis shows a more pronounced anteroposterior curvature, a more anteriorly projecting and superiorly placed pogonion on the external aspect, and a deeper genioglossal fossa, and more posteriorly projecting mental spine on the internal aspect. The right corpus at the P~3~/P~4~ is vertical and narrow towards the alveolar margin, but projects buccally towards the basal margin. The corpus also shows a pronounced depression and slightly thicker cortical bone on the lingual aspect and a superoinferior thickening of the basal corpus. The right corpus at the M~1~/M~2~ is also vertical and narrow towards the alveolar margin and projects buccally towards the basal margin, with a pronounced depression on the lingual basal aspect and a modest superoinferior thickening of the basal corpus. The left corpus at the P~3~/P~4~ is relatively vertical along the buccal side but convex inferolingually, and the cortical bone is thicker along the lingual corpus, widens at the mid-corpus and is widest at the basal corpus. The left M~1~/M~2~ takes on a general ‘U’ shape with a widening at the mid-corpus.
2.  Negative PC1 shapes represent mandibles with a taller mandibular corpus and symphysis. The symphysis also shows a less pronounced anteroposterior curvature and a less projecting and more inferiorly placed pogonion on the external aspect. Internally, the symphysis shows a more pronounced lingual projection (*planum alveolare*), a weaker mental spine and considerably thicker cortical bone. The right corpus at the P~3~/P~4~ shows a fairly vertical buccal side and greater curvature and thicker cortical bone on the lingual aspect. The right M~1~/M~2~ also shows a fairly vertical buccal side, a more pronounced depression on the lingual side and a modestly thickened basal margin. The left corpus at the P~3~/P~4~ and at the M~1~/M~2~ shows vertical buccal and lingual sides and modest thickening at the basal margin.
3.  Negative PC2 shapes resemble those of the positive loadings in PC1, and individuals with the most pronounced chins are clustered in positive PC1 and negative PC2 shape space. Negative shapes along PC2 represent more convex shapes at the regions with the greatest cortical bone thickness.
4.  Positive PC2 shapes differ from the negative shape loadings in PC1 in that the thickest distribution of cortical bone is more inferiorly placed along the anteramal and mylohyoid shelves, in the basal corpus, and at the symphysis in the regions of the mental spine and pogonion.

### View Allometric shape changes by group

```{r ASCG,echo=TRUE,eval=FALSE}
allo = ASCG(output, # GMA_object
            Group, # list object of grouping vectors
            coord_name = dimnames(Corpus_Land)[[1]], # names for the landmarks
            coord_subset = paste(rep(c("LM1M2", "LP3P4", "Symphysis", 
                                       "RP3P4", "RM1M2"), each = 118)), 
            key = "Allometric Shapes") 

# View overall allometric influence on all individuals in the sample
## Statistics
allo$Allo.size.reg$table
## Shape
allo$whole_shape

# View allometric influence of study sample by skeletal collection

# allo is the output
# Group_stats is the statistics
# the first [[1]] is for Group 1 in the list
# the second bracket is the specific group in the group you want

## Statistics
# St. Gregory's Priory
allo$Group_stats[[1]][[1]]$`Collection.St. Gregrory Priory`$table
#SMU
allo$Group_stats[[1]][[2]]$Collection.SMU$table
#UP 
allo$Group_stats[[1]][[3]]$Collection.UP$table

# all groups interactively
allo$Group_shape


```

### View Morphological shape changes by Collection

```{r Group_Shape,echo=TRUE,eval=FALSE}
Group_shapes = SCG(output,
                   Group,
                   coord_name = dimnames(Corpus_Land)[[1]],
                   subset = paste(rep(c("LM1M2", "LP3P4",
                                        "Symphysis", 
                                        "RP3P4", "RM1M2"), 
                                      each = 118)),
                   key = "Group Shapes")

# View the Group Differences

Group_shapes[[1]]

```

## Procrustes Shape (A)symmetry Analysis

### Obtain (A)symmetry results for directional and fluctuating asymmetry

```{r Bilat_Symmetry,echo=TRUE,eval=FALSE}
output$sym_data$shape.anova

output$sym_resid_data$shape.anova
```

Notes:

1.  Directional side asymmetry can explain over 90% of the shape symmetry between the left and right sides indicating a very strong presence of asymmetry.
2.  This signal is still present even when accounting for Procrustes Residuals.

### See how much of the directional and fluctuating asymmetry are influenced by similar factors between the right and left sides.

```{r Mantel_Symmetry,echo=TRUE,eval=FALSE}

output$sym_data$Mantel_tests$combined
```

Notes:

1.  The test is significant but only has a correlation coefficient of 0.10 indicating different factors are driving random variation versus the directional asymmetry.

### Visualize the asymmetric shape variation

```{r Asymmetric_PCA,echo=TRUE,eval=FALSE}

# 3D PCA Plot

output$PCA_asym$asym_PCA_plots$PC_1_2_3

# 3D PCA Plot By skeletal collection

output$PCA_asym$asym_PCA_groups$Collection$PC_1_2_3
```

Notes:

1.  There is a lot of overlap in skeletal collection, also noted in the cluster analysis.

### Visualize the asymmetric shape variation across PCs 1,2, and 3

```{r Asymmetric_Shapes,echo=TRUE,eval=FALSE}
#PC 1
output$PCA_asym$asym_lollipop_PCs[[1]]

#PC 2
output$PCA_asym$asym_lollipop_PCs[[2]]

#PC 3
output$PCA_asym$asym_lollipop_PCs[[3]]

```

Notes: The major influences of asymmetric shape variation detected by bilateral symmetry analysis is as follows

1.  An anterior migration of the posterior corpus (PC 1)
2.  A buccobasal projection of the right posterior corpus (PC 2)
3.  A basal elongation of the posterior corpus (PC 3)

### Does this side asymmetry influence the mandibular symphysis?

#### Need to run a simple Procrustes Linear regression of symphyseal shape by the unsigned asymmetry index using Geomorph's procD.lm

```{r Symph_Influence,echo=TRUE,eval=FALSE}
summary(geomorph::procD.lm(output$GPA$coords[c(237:354),,] ~ output$sym_data$unsigned.AI))
```

Notes:

1.  Yes, there is a **Minor** influence of side asymmetry on the mandibular symphysis.
2.  This asymmetry is largely independent of the asymmetry as seen in the PCA of the Procrustes Landmarks (this shape variation is not captured in total shape PCA).

# Step 3: Hypothesis Testing

## AABA - Analysis of Associations Between Attributes

### Create Hypothesis Models

![Figure 2. Workflow for this hypothesis test that the posterior mandibular corpus is biomechanically integrated with the symphysis.](vignettes/images/CSGM_Models.png)

### Create the Hypothesis Model following the structure in Figure 2.

```{r Data_Wrangling,echo=FALSE,eval=FALSE}
library(tidyverse)

#Data wrangle biomechanically relevant data for hypothesis testing

# Create a combined Posterior Corpus
Corpus <- Corpus_Prop[,-c(2,4)] %>% # Get rid of variables Collection and BM3
              filter(Region != "Symphysis") %>% # get rid of the symphysis
              pivot_wider( # make redundant rows into columns
                  id_cols = ID, # use specimen ID for data wrangling
                  names_from = Region, # we want to separate columns by region
                  values_from = -c(ID,Region)) %>%  # the values should be by specimen and region
              .[,-1] # get rid of the ID column

# Create the Left Posterior Corpus
Left_Corpus <- Corpus_Prop[,-c(2,4)] %>% # Get rid of variables Collection and BM3
                  filter(Region %in% c("LM1M2", "LP3P4")) %>% # only left side
                  pivot_wider( # make redundant rows into columns
                      id_cols = ID, # use specimen ID for data wrangling
                      names_from = Region, # we want to separate columns by region
                      values_from = -c(ID,Region)) %>%  # the values should be by specimen and region
                  .[,-1] # get rid of the ID column

# Create the Right Posterior Corpus
Right_Corpus <- Corpus_Prop[,-c(2,4)] %>% # Get rid of variables Collection and BM3
                  filter(Region %in% c("RM1M2", "RP3P4")) %>% # only right side
                  pivot_wider( # make redundant rows into columns
                      id_cols = ID, # use specimen ID for data wrangling
                      names_from = Region, # we want to separate columns by region
                      values_from = -c(ID,Region)) %>%  # the values should be by specimen and region
                  .[,-1] # get rid of the ID column

# Create the Hypothesis Model

Models <- 
  list(
    'Symphysis_Shape ~ Posterior_Corpus_Properties' = # Hypothesis Model 1
        list(
            'Symphyseal_Shape ~ Left_Corpus_Properties' = # Hypothesis Test 2
                list(
                     Symphyseal_Shape = output$GPA$coords[c(237:354),,],
                     Left_Corpus_Prop = Left_Corpus,
                     dt_parameters = list(Res_transform = "procdist+pca", 
                                          dt_method = "procdist", # Procrustes Distances
                                          Res_ncomp = 8, # we saw only 8 important PCs in factor analysis
                                          Pred_transform = "NULL"
                                          )
                    ),
              
              'Symphyseal_Shape ~ Left_Corpus_Shape' = # Hypothesis Test 2
                list(
                     Symphyseal_Shape = output$GPA$coords[c(237:354),,],
                     Left_Corpus_Shape = output$GPA$coords[c(1:236),,],
                     dt_parameters = list(Res_transform = "procdist+pca", 
                                          dt_method = "procdist", # Procrustes Distances
                                          Res_ncomp = 8, # we saw only 8 important PCs in factor analysis
                                          Pred_transform = "procdist+pca",
                                          Pred_ncomp = 8
                                          )
                    ),
              
              'Symphyseal_Shape ~ Right_Corpus_Properties' =  # Hypothesis Test 3
                list(
                     Symphyseal_Shape = output$GPA$coords[c(237:354),,],
                     Right_Corpus_Prop = Right_Corpus,
                     dt_parameters = list(Res_transform = "procdist+pca", 
                                          dt_method = "procdist", # Procrustes Distances
                                          Res_ncomp = 8, # we saw only 8 important PCs in factor analysis
                                          Pred_transform = "NULL"
                                          )
                    ),
            
              'Symphyseal_Shape ~ Right_Corpus_Shape' = # Hypothesis Test 2
                  list(
                       Symphyseal_Shape = output$GPA$coords[c(237:354),,],
                       Right_Corpus_Shape = output$GPA$coords[c(355:590),,],
                       dt_parameters = list(Res_transform = "procdist+pca", 
                                            dt_method = "procdist", # Procrustes Distances
                                            Res_ncomp = 8, # we saw only 8 important PCs in factor analysis
                                            Pred_transform = "procdist+pca",
                                            Pred_ncomp = 8
                                            )
                      )
              )
    )

```

### Create a subsetted Hypothesis Model

```{r Hypothesis_Model,echo=TRUE,eval=FALSE}
point_set <- 
  list(
    'Symphysis_Shape ~ Posterior_Corpus_Properties' = # Hypothesis Model 1
        list(
            'Symphyseal_Shape ~ Left_Corpus_Properties' = # Hypothesis Test 2
                list(
                     Symphyseal_Shape = list(Set1 = 1:118, 
                                     Set2 = 1:118, 
                                     Set3 = 1:118, 
                                     Set4 = 1:118),
                     Left_Corpus_Prop = list(
                                      Dimensions = c(1:6,9:16,19:30,17:18,47:48), 
                                      Biomechanics = c(37:46,49:58),
                                      Lingual_Thickness = c(59:118),
                                      Buccal_Thickness = c(117:176))
                    ),
              'Symphyseal_Shape ~ Left_Corpus_Shape' = # Hypothesis Test 2
                list(
                     Symphyseal_Shape = list(Set1 = 1:118),
                     Left_Corpus_Shape = list(Set1 = 1:118)
                    ),
              'Symphyseal_Shape ~ Right_Corpus_Properties' =  # Hypothesis Test 3
                list(
                     Symphyseal_Shape = list(Set1 = 1:118, 
                                     Set2 = 1:118, 
                                     Set3 = 1:118, 
                                     Set4 = 1:118),
                     Right_Corpus_Prop = list(
                                      Dimensions = c(1:6,9:16,19:30,17:18,47:48), 
                                      Biomechanics = c(37:46,49:58),
                                      Lingual_Thickness = c(59:118),
                                      Buccal_Thickness = c(117:176))
                    ),
              'Symphyseal_Shape ~ Right_Corpus_Shape' = # Hypothesis Test 2
                  list(
                       Symphyseal_Shape = list(Set1 = 1:118),
                       Right_Corpus_Shape = list(Set1 = 1:118)
                      )
              )
    )
```

Note: Symphyseal subsets don't change. Therefore, we shouldn't use paired subsets to conduct a test for every combination isn't needed. Only sequential pairings is computationally ideal for these tests.

### Conduct the Hypothesis Testing with AABA using correlation, covariation, and prediction methods

Note: Due to the robust operations within AABA, this may take a significant amount of time depending on model complexity. The complete, subsetted analysis will take approximately 10-15 minutes to complete.

```{r AABA,echo=TRUE,eval=FALSE}
AABA_results = AABA(Models,
                   point_set, 
                   paired_subsets = FALSE, # Don't do every combination
                   cors = TRUE, # Correlation analysis
                   vips = TRUE, #Partial Leaast Squares (PLS) analysis
                   covs = TRUE, # Shape specific 2B-PLS
                   PSLR = TRUE, # Procrustes Shape Linear Regression
                   regs = TRUE, # Conduct MARS nonlinear regression
                   VIP_trim = TRUE #trim regressions by VIPs
                   ) # Keep defaults or change them. See '?AABA'

```

#### Summarize all of the results of the AABA function

```{r AABA_Summary,echo=TRUE,eval=FALSE}
CSGM::summary(AABA_results)
```

## Individual Tests:

### In this section, we will consider each of the functions separately to demonstrate how they can be used outside of AABA and compare the results of each study in relation to the data.

### *cor.bio*: Conduct a correlation analysis between the posterior corpus and symphysis

```{r Correlation,echo=TRUE,eval=FALSE}
cor_hypothesize = cor.bio(Models, #Your hypothesis models
                          point_set, # Your subsets
                          paired_subsets = FALSE #Not every subset combination
                          )
```

Notes: Visualize the results of the correlograms in your directory under the folder "Correlations_Bio".

#### Obtain summary statistics for cor.bio function

```{r Correlation_Summary,echo=TRUE,eval=FALSE}
CSGM::summary(cor_hypothesize)
```

Notes: The results of the correlation tests demonstrate a few important aspects of the relationships between shape and cross-sectional properties across the mandibular symphysis.

1.  There is a significantly strong correlation between posterior corpus and the symphyseal shapes as demonstrated in the PCA. However, this is not supported by the canonical associations. However, this may also be due to the fact that PC data is orthogonal and thus canonical correlations do not reveal covarying trends.
2.  Overall, the correlations of the right and left side cross-sectional properties were moderately associated with symphyseal shape but not significantly so. While this may reveal that cross–sectional properties are not as important, it is important to consider that this may be dampened by irrelevant orthogonal axes of the PCA. Therefore, a may robust statistical appraoch is needed to confirm these relationships.
3.  Bending and breaking strengths of the posterior corpus do not seem to be strongly associated with the first or second principal components of symphyseal shape except for the minimum breaking strength in the premolars being associated with anterosuperiorly projected chins. While this would demonstrate a lack of association, the following analyses of covariation and prediction will reveal non-linear trends in these variables.

#### *pls.bio*: Conduct a covariation analysis between the posterior corpus and symphysis

```{r Covariation, echo=TRUE,eval=FALSE}
VIP_hypothesize = pls.bio(Models, # Your hypothesis models
                          point_set, # Your subsets
                          paired_subsets = FALSE, # Not every subset combo
                          vips = TRUE, # VIP plot generation
                          vip_method = "spls", # Sparce PLS method
                          lv_method = "mean", # VIPS > mean of all vips
                          cv_type = "CV", # repeated k-fold cross-validation
                          covs = TRUE # spls variable feature selection
                          )
```

#### Obtain summary statistics for Bio.VIP function

```{r Covariation_summary,echo=TRUE,eval=FALSE}
CSGM::summary(VIP_hypothesize)
 
```

Notes: The results of the `pls.bio` analysis are quite informative. Keep in mind our structure

Tests 1 and 3 - Symphysis Shape \~ Corpus Properties

-   Subset 1 - Corpus Dimensions

-   Subset 2 - Corpus Biomechanics

-   Subset 3 - Lingual Thicknesses

-   Subset 4 - Buccal Thicknesses

Tests 2 and 4 - Symphysis shape \~ Posterior Corpus shape

-   Subset 1 - All landmarks

The results reveal more nuanced relationships that couldn't be captured through canonical correlation analysis. For example, where the biomechanical properties of the posterior corpus revealed few correlations with symphyseal shape, the biomechanical properties of the left side significantly explains nearly 60% of symphyseal shape variation as well as the rest of the cross–sectional properties. The total squared correlation coefficient (R2) uses the diagonal of the correlation matrix, demonstrating that the influences of the cross-sectional properties on symphyseal shape follows a non-linear pattern. When examining the VIP selected variables, the maximum and minimum breaking strengths as well as mandibular torsion were the most relevant in predicting symphyseal shape. This is also supported in the correlation analysis, particularly in the premolars as well as in Keeling et al., 2025.

Furthermore, both sides of the mandibular corpus are strongly integrated with symphyseal shape, but the right side shares a stronger integration (Q2 = 0.96; RV = 0.82). This is supported by the results obtained in Keeling et al., 2025 and make sense given the asymmetric corpus results in the PCA.

While these results suggest a biomechanical and shape integration of the posterior corpus, and particularly the right side, with symphyseal shape, its important to further evaluate this using regression.

#### Conduct a Procrustes Shape Linear Regression between the posterior corpus and symphysis with `pslr.bio`

```{r PSLR,echo=TRUE,eval=FALSE}
PSLR_hypothesize = pslr.bio(Models, # Hypothesis model
                        point_set = point_set, # subsets
                        paired_subsets = FALSE, # not every combination
                        all_vips = NULL, # we dont want full model vip trimming
                        subset_vips = NULL # We don't want subset vip trimming
                        )
```

#### Obtain summary statistics for PSLR function

```{r PSLR_Summary,echo=TRUE,eval=FALSE}
CSGM::summary(PSLR_hypothesize)
```

Notes: The results of the PSLR regression test demonstrate the limits that linear regression has when handling multivariate data. As can be seen, the unbelievably high and low squared correlation coefficient (R2) is the result of model overfitting (for variables over 30) and linearity (as we already established there are non-linear relationships involved). However, PSLR was designed to work well with 3D shape data, as indicated in hypothesis tests 2 and 4. Here, while the VIP analysis performed better, we can still see that the right side is contributing much more variation than the left side.

#### Conduct a Multivariate Adaptive Regression Splines (MARS) regression between the posterior corpus and symphysis

Note: `reg.bio` is the only hypothesis testing function without the hypothesis model approach. This is primarily due to the computationally intense nature of `reg.bio` and to permit user flexibility in quick testing without a hypothesis model approach. However, this hypothesis model approach is builtin the function "AABA" to allow for seamless hypothesis testing with `reg.bio`.

Due to its rigorous statistical framework, we can use `reg.bio` without the need for subsetting, due to its automated regularization techniques, to obtain a complete idea about the relationships between the posterior corpus and symphysis.

```{r MARS_Regression,echo=TRUE,eval=FALSE}
pred_tests = list()

for(i in seq_along(Models[[1]])) { # Since we know we have 1 hypothesis model
  pred_tests[[i]] <- reg.bio(Response = Models[[1]][[i]][[1]], # Response Data
                          Predictor = Models[[1]][[i]][[2]], # Predictor Data
                          type = "pred", # Prediction not classification
                          p_calc = "gcv", # generalized cross-validation stat
                          method = "gcvEarth", # Pred or Class model
                          n_boot = 25, # 50 rounds of repeated CV
                          perm = 500, #500 possible rounds of permutation
                          CV = TRUE, # repeated k-fold cross-validation
                          full_model = TRUE, # Trim variables 
                          which_classify = NULL, # for classification
                          VIPS = NULL, # no vip selection
                          Res_transform = Models[[1]][[i]][[3]]$Res_transform,
                          Res_ncomp = Models[[1]][[i]][[3]]$Res_ncomp,
                          Pred_transform = Models[[1]][[i]][[3]]$Pred_transform,
                          dt_method = Models[[1]][[i]][[3]]$dt_method)
                        
  
}
```

#### Obtain summary statistics of the multivariate adaptive regression splines

```{r,echo=TRUE,eval=FALSE}

# Hypothesis test 1 
CSGM::summary(pred_tests[[1]])  
# Hypothesis test 2 
CSGM::summary(pred_tests[[2]])  
# Hypothesis test 3 
CSGM::summary(pred_tests[[3]])  
# Hypothesis test 4 
CSGM::summary(pred_tests[[4]])
```

Notes: As we can see, the MARS predictions outperformed the linear predictions, revealing the predictive strengths noted in the PLS analysis. Two very different methods converging on the same answer is a strong sign that the predictive strength of the model is accurately capturing real shape variation.

## End of *CSGM* Package Demonstration
