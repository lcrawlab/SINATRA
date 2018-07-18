# SINATRA 
Sub-Image Analysis using Topological Summary Statistics.

## Introduction
The Sub-Image Selection problem is to identify the regions of a collection of three dimensional objects that are most causal for predicting a given response.

SINATRA is a pipeline for analyzing the sub image problem using topological summary statistics. The method can be decomposed into the following steps:

1. Represent the object in question by a collection of vectors, by topological summary statistics such as the [Smooth Euler Characteristic Transform](https://arxiv.org/abs/1611.06818). 
2. Perform a variable selection algorithm on the transformed shape, such as Lasso.
3. Reconstruct the regions corresponding to the selected features from the second step.


As an application of our method, we study the problem of classifying molars of New World Monkeys by diet; using SINATRA, we identify the most important regions on these teeth for predicting the diet of these New World Monkeys. The dataset was first curated by Doug M. Boyer and can be found [here](). 

## Package Requirements

While any feature selection algorithm can be used, we use the procedure called [RATE](https://github.com/lorinanthony/RATE), a variable selection algorithm for nonlinear models using Relative Centrality. Other alternatives considered in the text include [Elastic Net](https://cran.r-project.org/web/packages/elasticnet/elasticnet.pdf) and Bayesian variable selection using variational inference ([varbvs](https://cran.r-project.org/web/packages/varbvs/index.html)). 

## Code Usage
While the usage of the code will be explained in the simulation scripts, we provide an overview of the general pipeline here.
### Finding Critical Points

### Feature Selection

### Reconstruction

### Plotting 

## Relevant Citations

