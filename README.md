# SINATRA 

Sub-Image Analysis using Topological Summary Statistics.

## Introduction

The sub-image selection problem is to identify physical regions that most explain the variation between two classes of three dimensional shapes.

SINATRA is a statistical pipeline for carrying out sub-image analyses using topological summary statistics. The algorithm follows four key steps:

1. 3D shapes (represented as triangular meshes) are summarized by a collection of vectors (or curves) detailing their topology (e.g. [Smooth Euler Characteristic Transform](https://arxiv.org/abs/1611.06818)). 
2. A statistical model is used to classify the shapes based on their topological summaries.
3. After itting the model, an association measure is computed for each topological feature (e.g. centrality measures, posterior inclusion probabilities, p-values, etc.).
4. Association measures are mapped back onto the original shapes via a reconstruction algorithm — thus, highlighting evidence of the physical (spatial) locations that best explain the variation between the two groups.

Through detailed simulations, we assess the power of our algorithm as a function of its inputs. Lastly, as an application of our pipeline, we conduct feature selection on a dataset consisting of mandibular molars from different genera of New World Monkeys and examine the physical properties of their teeth that summarize their phylogeny. 

## Package Details and Requirements

While our pipeline is flexible and any feature selection algorithm can be implemented, we use the relative centrality criterion ([RATE](https://github.com/lorinanthony/RATE)), which is a variable selection measure for nonlinear and nonparametric statistical methods (Crawford et al. 2019, Ish-Horowicz et al. 2019). Alternative methods implemented include the [elastic net](https://cran.r-project.org/web/packages/elasticnet/elasticnet.pdf) (Zou and Hastie 2003) and Bayesian variable selection using variational inference ([varbvs](https://cran.r-project.org/web/packages/varbvs/index.html)) (Carbonetto, Zhou, and Stephens 2017). 

The SINATRA software also requires the use of the package `rgl`, which requires X11 or XQuartz on Mac.

## R Package Download

To download the package, install [devtools](https://cran.r-project.org/web/packages/devtools/index.html) and run the command

	devtools::install('software') 
	
Next, to load the package, use the command

	library(sinatra)

## Code Usage

Details of our implementation choices for the SINATRA algorithm are detailed below.

### Topological Summary Statistics for 3D Shapes

We convert our dataset into a matrix by measuring topological summary statistics of the shapes in the dataset. After picking a set of direction on which to measure the Euler characteristic of our shapes, we run the function

	compute_standardized_ec_curve.
	
The outputted EC curve can then be transformed - either smoothened or differentiated - by the function 

	update_ec_curve.
	
For every shape in the dataset, amalgamate the EC curves measured in every direction in our set of directions. Stack these curves into a concatenated data matrix, where each row represents a different shape in our dataset.

###  Inference using a Gaussian Process model and RATE for variable selection.

We are interested with classifying two classes of shapes. For this purpose we set up a Gaussian Process model relating the EC curve matrix to the class labels. The EC curve matrix is comprised of *features* which are exactly the Euler characteristic values of sublevel sets in various directions.

To perform variable selection on this model, we use a measure of centrality called RATE (Crawford, 2019) -- each feature is assigned a normalized value, representing its importance.

Thresholding the features above a certain RATE value then determines a feature selection method.

The inference step and RATE step are done within the function

	find_rate_variables_with_other_sampling_methods

### Reconstruction of the selected subimage.

We then map back the selected features via the function

	compute_selected_vertices_cones.
	
This generates the subimage on a given shape in our dataset.

Alternatively, the function

	reconstruct_vertices_on_shape
	
generates a heatmap on the shape which can be interpreted as the importance of each subset of the shape with respect to the outcome, in this case the class label.
 
### Vignettes
Usage of the code is best understood by viewing the examples in `software/vignettes`. We provide examples of running the full pipeline on the following cases:

- Simulated shapes as in the manuscript.
- Imported OFF files of very different shapes.
- Generating power curves for simulated shapes.
- Datasets of primate molars.



## Relevant Citations

