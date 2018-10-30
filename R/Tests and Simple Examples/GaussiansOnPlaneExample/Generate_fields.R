library(mvtnorm)

## Create functions that make up the random field
random_field_class_one <- function(x){
  f1 <- function(x) dmvnorm(x,c(0.25,0.25),sigma = diag(2)*0.01)
  f2 <- function(x) dmvnorm(x,c(0.75,0.75),sigma = diag(2)*0.01)
  f3 <- function(x) dmvnorm(x,c(0.5,0.5),sigma = diag(2)*0.01)
  f4 <- function(x) dmvnorm(x,c(0.2,0.9),sigma = diag(2)*0.01)
    
  f1(x) + f2(x) + f3(x) - f4(x)
}

random_field_class_two <- function(x){
  f1 <- function(x) dmvnorm(x,c(0.25,0.25),sigma = diag(2)*0.01)
  f2 <- function(x) dmvnorm(x,c(0.75,0.75),sigma = diag(2)*0.01)
  f3 <- function(x) dmvnorm(x,c(0.5,0.5),sigma = diag(2)*0.01)
  f4 <- function(x) dmvnorm(x,c(0.8,0.1),sigma = diag(2)*0.01)
    
  f1(x) + f2(x) + f3(x) - f4(x) 
}

# Get a discretization of the mesh / graph of the function; add randomness within each class
generate_matrix_class_one <- function(grid_length){
  x <- seq(0,1,1/(grid_length-1))
  y <- seq(0,1,1/(grid_length-1))
  grid <- expand.grid(x,y)
  result <- matrix(random_field_class_one(grid),nrow = grid_length, byrow = TRUE)
  result <- result + matrix(runif(grid_length^2,min=-0.1,max=0.1),nrow = grid_length)
}

generate_matrix_class_two <- function(grid_length){
  x <- seq(0,1,1/(grid_length - 1))
  y <- seq(0,1,1/(grid_length - 1))
  grid <- expand.grid(x,y)
  result <- matrix(random_field_class_two(grid),nrow = grid_length, byrow = TRUE)
  result <- result + matrix(runif(grid_length^2,min=-0.1,max=0.1),nrow = grid_length)
}

create_data_gaussian <- function(num_sim, grid_size=20,curve_length = 25,dir){
  #Inputs:  num_sim: observations for each class
  #         grid_size (int) : number of training points on [0,1]
  #         length_scale (int) : parameter for RBF kernel
  #         dir (kxm) array: of directions
  #Ouputs: data: (2*num_sim x m array): of all the generated EC curves
  # create the GP object.
  # m here is the dimension of Euclidean space we are working in.
  data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]) )
  for (i in 1:num_sim){
    #new_data=generate_clustered_points_with_means(data_1,n=5,eta=0)
    m=generate_matrix_class_one(grid_size)
    complex <- MatrixtoSimplicialComplex(m)
    ec_curve <- matrix(NA,nrow = 1,ncol=0)
    for (j in 1:dim(dir)[1]){
      vertex_function <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,dir[j,1],dir[j,2],dir[j,3]))
      curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = TRUE)
      curve <- integrate_ec_curve(curve)
      # omit the length data, for now
      ec_curve <- c(ec_curve,curve[,2])
    }
    data <- rbind(data,c(1,ec_curve))
  }
  print('On class 2')
  #Second Class
  for (i in 1:num_sim){
    #new_data=generate_clustered_points_with_means(data_2,n=5,eta=0)
    m=generate_matrix_class_two(grid_size)
    complex <- MatrixtoSimplicialComplex(m)
    ec_curve <- matrix(NA,nrow = 1,ncol=0)
    for (j in 1:dim(dir)[1]){
      vertex_function <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,dir[j,1],dir[j,2],dir[j,3]))
      curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = TRUE)
      curve <- integrate_ec_curve(curve)
      # omit the length data, for now
      ec_curve <- c(ec_curve,curve[,2])
    }
    data <- rbind(data,c(-1,ec_curve))
  }
  return(data)
}

# Visualizes the random field in an OFF file, and outputs the computed mesh.
visualize_random_field <- function(filename,type){
  if(type == 1){
    posterior_matrix1=generate_matrix_class_one(grid_length)
  }else{
    posterior_matrix1=generate_matrix_class_two(grid_length)
  }
  posterior_matrix1=generate_matrix_class_two(grid_length)
  posterior_complex1 <- MatrixtoSimplicialComplexTriangular(posterior_matrix1,grid_length)
  try_vertex1=posterior_complex1$Vertices[,1:3]
  try_faces1=posterior_complex1$Faces
  ind1 = rep(1,dim(try_vertex1)[1])
  try_vertex1 = cbind(try_vertex1,ind1)
  
  grf_mesh_posterior1=tmesh3d(t(try_vertex1),indices=t(try_faces1))
  plot3d(grf_mesh_posterior1,col = alpha.col('grey',alpha=0.1))
  aspect3d(1,1,1)
  
  #Write to one OFF file
  vcgOffWrite(grf_mesh_posterior1,filename=filename)
  grf_mesh_posterior1
}

visualize_random_field_triangular <- function(filename){
  posterior_matrix1=generate_matrix_class_one(grid_length)
  posterior_complex1 <- MatrixtoSimplicialComplexTriangular(posterior_matrix1)
  
  vcgOffWrite(posterior1,filename=filename)
  grf_mesh_posterior1
  
}


