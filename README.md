# SINATRA 

Sub-Image Analysis using Topological Summary Statistics.

## Introduction

The sub-image selection problem is to identify physical regions that most explain the variation between two classes of three dimensional shapes. SINATRA is a statistical pipeline for carrying out sub-image analyses using topological summary statistics. The algorithm follows four key steps:

1. 3D shapes (represented as triangular meshes) are summarized by a collection of vectors (or curves) detailing their topology (e.g. Euler characteristics, persistence diagrams). 
2. A statistical model is used to classify the shapes based on their topological summaries. Here, we make use of a Gaussian process classification model with a probit link function.
3. After itting the model, an association measure is computed for each topological feature (e.g. centrality measures, posterior inclusion probabilities, p-values, etc).
4. Association measures are mapped back onto the original shapes via a reconstruction algorithm — thus, highlighting evidence of the physical (spatial) locations that best explain the variation between the two groups.

Through detailed simulations, we assess the power of our algorithm as a function of its inputs. Lastly, as an application of our pipeline, we conduct feature selection on a dataset consisting of mandibular molars from different genera of New World Monkeys and examine the physical properties of their teeth that summarize their phylogeny. 

## Package Details and Requirements

Code for implementing the SINATRA pipeline was written in R (version 3.5.3). As part of this procedure:

1. Inference for the Gaussian process classification (GPC) model was done using elliptical slice sampling (Murray, Prescott, and MacKay 2010) and carried out with the R package [FastGP](https://cran.r-project.org/web/packages/FastGP/index.html) (version 1.2).
2. Next association measures are computed for the Euler characteristic curves. While our pipeline is flexible and any feature selection algorithm can be implemented, we use the relative centrality criterion ([RATE](https://github.com/lorinanthony/RATE)), which is a variable selection measure for nonlinear and nonparametric statistical methods (Crawford et al. 2019; Ish-Horowicz et al. 2019). Alternative methods implemented include the [elastic net](https://cran.r-project.org/web/packages/elasticnet/elasticnet.pdf) (Zou and Hastie 2003) and Bayesian variable selection using variational inference ([varbvs](https://cran.r-project.org/web/packages/varbvs/index.html)) (Carbonetto, Zhou, and Stephens 2017).
3. Visualization of reconstructed regions outputted by SINATRA is done using the package [rgl](https://cran.r-project.org/web/packages/rgl/index.html) (version 0.100.19), and general utility functions for triangular meshes from the package [Rvcg](https://cran.r-project.org/web/packages/Rvcg/index.html) (version 0.18). 

Note that the package `rgl` requires X11 or XQuartz on macOS systems.

## R Package Download

To download the package, install [devtools](https://cran.r-project.org/web/packages/devtools/index.html) and run the command:

	devtools::install('software') 
	
Next, to load the package, use the command:

	library(sinatra)

Other common installation procedures may also apply.

## Code Usage

Other details of our implementation choices for the SINATRA algorithm are provided below.

### Topological Summary Statistics for 3D Shapes

In the first step of the SINATRA pipeline, we use a tool from differential geometry called the Euler characteristic (EC) transform (Turner, Mukherjee, and Boyer 2014; Crawford et al. 2017; Ghrist, Levanger, and Mai 2018) to represent 3D shapes as a collection of vector-valued topological summary statistics. To do so, after picking a set of directions on which to measure the ECs of each shape in our data, the algorithm runs the function `compute_standardized_ec_curve`.
	
If desired, resulting EC curve can then be transformed --- either smoothened or differentiated --- by using the function `update_ec_curve`.
	
For each shape in the dataset, EC curves are computed in every direction and then concatenated into a *p*-dimensional topological feature vector. For a study with *n*-shapes, we analayze an *n × p* design matrix, where the columns denote the Euler characteristic computed at a given filtration step and direction combination.

###  Gaussian Process Classification and Variable Selection

Recall we are interested in identifying physical features that best differentiate two classes of shapes. For this purpose, we use Gaussian process classification model to assess the relationship between topological summary statistics and the variance between class labels. To perform variable selection on these statistics, we use relative centrality measures: a criterion which evaluates how much information about the shape classification is lost when a particular topological feature is missing from the model.

Keeping ECs with centrality measures above a certain threshold then determines a selected set worth further investigation.

### Reconstruction of the Selected Sub-Image

After obtaining a select set of topological features, we map this information back onto the physical shape using the function `compute_selected_vertices_cones`. This generates the sub-image on a given shape in our dataset.

Alternatively, the function `reconstruct_vertices_on_shape` generates a heatmap on the shape which can be interpreted as visualizing the importance of each subset of Euler characteristics with respect to the class labels.
 
### Available Vignettes

Implementation of the code may be best understood by viewing the examples in `software/vignettes`. We provide tutorials for running the full SINATRA pipeline via the following cases:

- Simulations with perturbed spheres and shape caricaturizations;
- Simulations with shapes generated under the null hypothesis of relative centrality;
- How to generate power curves for simulated shapes;
- Real data analyses with primate molars.

Other code specific to analyses conducted in the paper can be found in the branch `SINATRA/paper_results`.

## Questions and Feedback

For questions or concerns, please contact [Bruce Wang](mailto:bruce.waang55@gmail.com), [Timothy Sudijono](mailto:timothy_sudijono@brown.edu), or [Lorin Crawford](mailto:lorin_crawford@brown.edu). We appreciate any feedback you may have with our repository and instructions.

## Relevant Citations

B. Wang*, T. Sudijono*, H. Kirveslahti*, T. Gao, D.M. Boyer, S. Mukherjee, and L. Crawford. SINATRA: a sub-image analysis pipeline for selecting features that differentiate classes of 3D shapes. _bioRxiv_. 701391.
