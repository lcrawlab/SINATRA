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
While the usage of the code will be explained in the simulation scripts, we provide an overview of the general pipeline here. We demonstrate the first steps of the pipeline on a single shape. Let's first import a shape to work with (which in this case is a bean). The file can be found in the data folder. See also the [Princeton Shape Benchmark](http://shape.cs.princeton.edu/benchmark/index.cgi) for more files. The following code takes in an OFF file an converts it to a simplicial complex, a R List data structure. 

	bean <- process_off_file_v3('bean_off_file')
	
`bean$Vertices` will return the coordinates of the vertices of the bean; `bean$Faces, bean$Edges` return the indices of the vertices that make up the face or edge respectively. We also need to create a plot of the object for visualization purposes.

	bean_plot <- vcgImport('bean_off_file')

To open a 3D plot, we can use:

### Finding Critical Points
The SECT transforms a shape, represented as a simplicial complex, into a collection of vectors indexed by the input directions. That is for each of these directions we compute the Euler Characteristic curves (EC curve) of the shape in that direction. We first can specify a set of directions manually or by calling a function in the package:

	directions <- matrix(c(1,0,0),ncol=3, byrow = TRUE)
	directions <- generate_equidistributed_directions(N)

We define the critical points as the points of the complex at which the EC curve changes value, for one of the specified directions above. We get these with the functions

	critical_points <- find_critical_points_multiple_directions(directions,bean,curve_length,phi)
	critical_point_coordinates <- get_all_critical_points(critical_points)

The second line is a simple function; it just extracts the coordinates of the critical points.


As parameters to the reconstruction function, we have `phi` and `curve_length`. Curve length dictates the length of the ec curve. Phi does not have to be specified, but is the angle between an input direction and its perturbed directions (which are used to reconstruct the critical points).

### Feature Selection

### Reconstruction

### Plotting 

## Relevant Citations

