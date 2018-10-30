source('R/Feature_Selection_Functions.R')
source('R/SECT_Functions.R')
source('R/Simulation_Functions.R')
source('R/Tests and Simple Examples/GaussiansOnPlaneExample/Generate_fields.R')
library(plot3D)
library(Rcpp)
library(RcppArmadillo)
library(kernlab)
library(doParallel)
sourceCpp("R/BAKRGibbs.cpp")
source('R/RATE.R')
library(varbvs)


# specify the critical points for the flat plane with gaussians example.
points_1 = matrix(c(0.25,0.25,random_field_class_one(c(0.25,0.25)),
                    0.75,0.75,random_field_class_one(c(0.75,0.75)),
                    0.5,0.5,random_field_class_one(c(0.5,0.5)),
                    0.2,0.9,random_field_class_one(c(0.2,0.9))),nrow = 4, byrow = TRUE)
points_2 = matrix(c(0.25,0.25,random_field_class_two(c(0.25,0.25)),
                    0.75,0.75,random_field_class_two(c(0.75,0.75)),
                    0.5,0.5,random_field_class_two(c(0.5,0.5)),
                    0.8,0.1,random_field_class_two(c(0.8,0.1))),nrow = 4, byrow = TRUE)

#generate the desired directions
#directions = generate_equidistributed_points_hemisphere(25)
directions = matrix(c(0,0,1,
                      1/sqrt(3),1/sqrt(3),1/sqrt(3),
                      -1/sqrt(3),1/sqrt(3),1/sqrt(3),
                      1/sqrt(3),-1/sqrt(3),1/sqrt(3),
                      -1/sqrt(3),-1/sqrt(3),1/sqrt(3)),
                      ncol = 3, byrow = TRUE)

directions <- rbind(directions,-directions)

grid_length = 40
num_sim = 50

data = create_data_gaussian(num_sim = 50,grid_size = grid_length, dir = directions)

#### Visualize some GPs from each class ####

grf_mesh_posterior1 = visualize_random_field('GP_1',1)
grf_mesh_posterior2 = visualize_random_field('GP_2',2)



class_one_matrix <- generate_matrix_class_one(grid_length)

#Transform it into a simplicial complex
class_one_complex <- MatrixtoSimplicialComplex(class_one_matrix)

# Look at Ec curves
vertex_function <- class_one_complex$Vertices%*%c(0,0,0,1)
compute_discrete_ec_curve(class_one_complex,vertex_function,100,first_column_index = TRUE)



### Feature Selection and Reconstruction ###
# We use Rate. Fit a Gaussian process classifier to this data
ind <- sample(1:dim(data)[1],30)
gpc(data,ind)

#### Feature selection for indices, using RATE ####
want_indices=find_rate_variables(data,radius = 2,bandwidth = 0.02);


