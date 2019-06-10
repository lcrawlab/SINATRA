# SINATRA 
Sub-Image Analysis using Topological Summary Statistics.

## Introduction
The Sub-Image Selection problem is to identify the regions of a collection of three dimensional objects that are most causal for predicting a given response.

SINATRA is a pipeline for analyzing the sub image problem using topological summary statistics. The method can be decomposed into the following steps:

1. Represent the object in question by a collection of vectors, by topological summary statistics such as the [Smooth Euler Characteristic Transform](https://arxiv.org/abs/1611.06818). 
2. Perform a variable selection algorithm on the transformed shape.
3. Reconstruct the regions corresponding to the selected features from the second step.


As an application of our method, we study the problem of classifying molars of New World Monkeys by diet; using SINATRA, we identify the most important regions on these teeth for predicting the diet of these New World Monkeys. The dataset was first curated by Doug M. Boyer. 

## Package Requirements

While any feature selection algorithm can be used, we use the procedure called [RATE](https://github.com/lorinanthony/RATE), a variable selection algorithm for nonlinear models using Relative Centrality. Other alternatives considered in the text include [Elastic Net](https://cran.r-project.org/web/packages/elasticnet/elasticnet.pdf) and Bayesian variable selection using variational inference ([varbvs](https://cran.r-project.org/web/packages/varbvs/index.html)). 

The package also requires the use of the package `rgl`, which requires X11 or XQuartz on Mac.

## Download

To download the package, checkout the repo and within the directory, run the command

	devtools::install('software') 
	
and load the package with

	library(sinatra).

## Code Usage
The SINATRA pipeline can be divided into several steps:

### Data Vectorization
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

