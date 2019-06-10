#library(Rcpp)
#library(mvtnorm)
#library(pdist)
#library(MASS)
#library(plyr)
#library(reshape2)
#library(FNN)

# cusps are the equidistributed points on the sphere - the 'grid' points.

#' Generate (S/D) EC curves for a collection of perturbed spheres
#'
#' @export
#' @import Rvcg
#' @import pdist
#' @import FNN
#'
#' @description Given a set of input parameters, \code{create_data_sphere_simulation} generates replicates of perturbed spheres.
#' Cusps are initialized by generating equidistributed points on the sphere. The cusps are created by taking the $k$ nearest vertices of these central points,
#' and depending on if the cusps are causal or shared, the cusps are indented or made protruding out.
#' These perturbed spheres are divided into different classes using the indented cusps.
#' These spheres then have the (S/D) EC curve computed over to create the design matrix $X$.
#'
#' @param num_sim (int) : The number of replicates of data.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param dir (matrix) : The matrix of directions we compute the EC curves over.
#' @param noise_points (int) : The number of shared points in each shared cusp.
#' @param causal_points (int) : The number of class specific causal points in each causal cusp.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param ec_type (string) : The type of EC we are computing. We currently support ECT, DECT and SECT.
#' @param subdivision (int) : The fineness of the mesh. We currently use subdivision = 3.
#' @param cusps (int) : The number of total cusps to be generated.
#' @param causal_regions_1 (vector) : The index of cusps to be denoted as the central points for the causal cusps of class 1.
#' @param causal_regions_2 (vector) : The index of cusps to be denoted as the central points for the causal cusps of class 2.
#' @param shared_regions (vector) : The index of cusps to be denoted as the central points for the shared cusps for both classes.
#'
#' @return data_list (list) : List of procedure specific data : The associated (S/D) EC matrix, indices of the shared & causal vertices,
#' and the metadata for locations of the cusps.
generate_data_sphere_simulation = function(nsim, curve_length, dir, noise_points = 5,causal_points = 5, ball = TRUE, ball_radius = 2,
                                           ec_type = "ECT",subdivision = 3, cusps, causal_regions_1 = c(1), causal_regions_2 = c(3),
                                           shared_regions = c(4)) {
  data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]) )
  regions =  generate_equidistributed_points(cusps,cusps)

  #Initiate the causal points
  sphere = vcgSphere(subdivision = subdivision)

  #closest_points_class2 = closest_points_class1
  print(paste('Causal Regions 1: '))
  print(causal_regions_1)
  print('Causal Regions 2: ')
  print(causal_regions_2)
  print('Shared Regions: ')
  print(shared_regions)

  complex_points = list()
  shared_points_list = list()
  total_shared_points = c()
  total_closest_points_class1 = c()
  total_closest_points_class2 = c()
  total_cusps_list = c()

  # dictionary; keys are the locations of the cusps, and the values are the vertices on the sphere closest to that cusp
  region_vertex_dictionary <- vector("list",dim(regions)[1])

  sphere_vertices <- asEuclidean(t(sphere$vb))

  #get distances between regions and vertices
  distances <- as.matrix(pdist(regions,sphere_vertices))

  for (i in 1:(dim(sphere_vertices))[1]){
    closest_region <- which.min(distances[,i])
    region_vertex_dictionary[[closest_region]] <- c(region_vertex_dictionary[[closest_region]],i)
  }

  vertex_region_dictionary <- apply(distances,2,FUN = which.min)


  ### Get the causal and shared regions on the sphere ###

  # iterate through class 1
  for (j in 1:length(causal_regions_1)){
    causal_dir1 = regions[causal_regions_1[j],]
    closest_points_class1 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir1,ncol = 3), k = causal_points)
    total_closest_points_class1 = c(total_closest_points_class1,closest_points_class1)
    total_cusps_list[[length(total_cusps_list) + 1]] = closest_points_class1

  }

  # iterate through class 2
  for (j in 1:length(causal_regions_2)){
    causal_dir2 = regions[causal_regions_2[j],]
    closest_points_class2 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir2,ncol = 3), k = causal_points)
    total_closest_points_class2 = c(total_closest_points_class2,closest_points_class2)
    total_cusps_list[[length(total_cusps_list) + 1]] = closest_points_class2
  }

  for (k in 1:length(shared_regions)){
    shared_dir = regions[shared_regions[k],]
    closest_points_shared = knnx.index(data = t(sphere$vb[-4,]),query = matrix(shared_dir,ncol = 3), k = noise_points)
    total_shared_points = c(total_shared_points,closest_points_shared)
  }

  ### Create the data ###
  for (i in 1:nsim){
    sphere1 = vcgSphere(subdivision = subdivision)
    sphere2 = vcgSphere(subdivision = subdivision)

    # Add noise to the sphere
    sphere1$vb[1:3,] = sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.035)
    sphere2$vb[1:3,] = sphere2$vb[1:3,]  * rnorm(dim(sphere2$vb)[2], mean = 1, sd = 0.035)

    # Descend the causal regions - Needs to be changed
    for (j in 1:length(causal_regions_1)){
      causal_dir1 = regions[causal_regions_1[j],]
      closest_points_class1 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir1,ncol = 3), k = causal_points)
      sphere1$vb[1:3,closest_points_class1] = sphere1$vb[1:3,closest_points_class1]  * 0.55 + rnorm(1, mean = 0, sd = 0.1)
    }

    for (j in 1:length(causal_regions_2)){
      causal_dir2 = regions[causal_regions_2[j],]
      closest_points_class2 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir2,ncol = 3), k = causal_points)
      sphere2$vb[1:3,closest_points_class2] = sphere2$vb[1:3,closest_points_class2]  * 0.55 + rnorm(1, mean = 0, sd = 0.1)
    }

    # Elevate the shared regions - Needs to be changed
    for (k in 1:length(shared_regions)){
      shared_dir = regions[shared_regions[k],]
      closest_points_shared = knnx.index(data = t(sphere$vb[-4,]),query = matrix(shared_dir,ncol = 3), k = noise_points)
      shared_points = sphere$vb[1:3,closest_points_shared]  * 1.35 + rnorm(1, mean = 0, sd = 0.1)
      sphere1$vb[1:3,closest_points_shared] = shared_points
      sphere2$vb[1:3,closest_points_shared] = shared_points

    }


    sphere_mesh1 = convert_off_file(sphere1)
    sphere_mesh2 = convert_off_file(sphere2)

    complex_points[[(2*i-1)]] = t(sphere1$vb[1:3,])
    complex_points[[2*i]] = t(sphere2$vb[1:3,])

    shared_points_list[[i]] = shared_points

    ec_curve_class1 <- matrix(NA,nrow = 1,ncol=0)
    ec_curve_class2 <- matrix(NA,nrow = 1,ncol=0)

    ### compute EC curves for both classes of curves
    for (j in 1:dim(dir)[1]){

      vertex_function_class_1 <- sphere_mesh1$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
      vertex_function_class_2 <- sphere_mesh2$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])

      if (ball == TRUE){

        curve1 <- compute_standardized_ec_curve(sphere_mesh1, vertex_function_class_1, curve_length-1, first_column_index = FALSE,ball_radius)
        curve2 <- compute_standardized_ec_curve(sphere_mesh2, vertex_function_class_2, curve_length-1, first_column_index = FALSE,ball_radius)

      }
      else{

        curve1 <- compute_discrete_ec_curve(sphere_mesh1, vertex_function_class_1, curve_length-1, first_column_index = FALSE)
        curve2 <- compute_discrete_ec_curve(sphere_mesh2, vertex_function_class_2, curve_length-1, first_column_index = FALSE)

      }

      # transform the ECT as desired
      curve1 <- update_ec_curve(curve1, ec_type)
      curve2 <- update_ec_curve(curve2, ec_type)

      # omit the length data, for now
      ec_curve_class1 <- c(ec_curve_class1,curve1[,2])
      ec_curve_class2 <- c(ec_curve_class2,curve2[,2])
    }

    data <- rbind(data,c(1,ec_curve_class1))
    data <- rbind(data,c(-1,ec_curve_class2))

  }

  #print(paste('Correlation of Central Direction: ',median(cor(t(data[,2:curve_length+1])))))
  #print(paste('Correlation of Matrix: ',median(cor(t(data[,-1])))))


  data_list=list(data = data, noise = total_shared_points, causal_points1 = total_closest_points_class1,
                 causal_points2 = total_closest_points_class2,complex_points = complex_points,
                 shared_points_list = shared_points_list,total_cusps_list = total_cusps_list,
                 region_vertex_dict = region_vertex_dictionary, vertex_region_dict = vertex_region_dictionary)

  return(data_list)
}
#' Generate (S/D) EC curves for a collection of interpolations on the plane
#'
#' @export
#' @import truncnorm
#' @description Given a set of input parameters, \code{create_data_normal_fixed} generates replicates of interpolations on the plane.
#' These interpolations are then converted to simplicial complexes, and have the (S/D) EC curve computed over them.
#'
#' @param num_sim (int) : The number of replicates of data.
#' @param dir (matrix) : The matrix of directions we compute the EC curves over.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param shared_points (int) : The number of shared points across both classes of shapes.
#' @param causal_points (int) : The number of class specific causal points in each class.
#' @param grid_size (int) : The fine-ness/granularity of the interpolated shapes.
#' @param func (function) : The function used to compute the interpolations.
#' @param eta (float) : The kernel shape parameter.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param ec_type (string) : The type of EC we are computing. We currently support ECT, DECT and SECT.
#'
#' @return data_list (list) : List of procedure specific data : The associated (S/D) EC matrix, coordinates of causal & shared points,
#' and the vertex coordinates of all the shapes.
#'
#'
create_data_normal_fixed = function(num_sim=25,dir,curve_length=10,shared_points=5,causal_points=5,
                                  grid_size=25,func=rbf_gauss,eta=5,ball = TRUE, ball_radius = 1,
                                  ec_type = "ECT"){

  data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]) )
  #Shared points
  n1=rtruncnorm(shared_points,a=0,b=1,sd=1)
  n2=rtruncnorm(shared_points,a=0,b=1,sd=1)
  #causal points
  x1=rtruncnorm(causal_points,a=-1,b=0,-0.5,sd=1)
  y1=rtruncnorm(causal_points,a=0,b=1,0.5,sd=1)
  x2=rtruncnorm(causal_points,a=0,b=1,mean=0.5,sd=1)
  y2=rtruncnorm(causal_points,a=-1,b=0,mean=-0.5,sd=1)
  noise_points=cbind(n1,n2)
  causal_points1=cbind(x1,y1)
  causal_points2=cbind(x2,y2)
  complex_points=list()
  for (i in 1:num_sim){
    total_complex=generate_normal_complex_fixed(grid_size=grid_size,noise_points=noise_points,causal_points1=causal_points1,
                                                causal_points2=causal_points2, func=func,eta=eta)
    if(inherits(total_complex,'try-error')){
      i=i-1
      next
    }
    else{
      complex_points[[(2*i-1)]] = total_complex[[6]]
      complex_points[[2*i]] = total_complex[[7]]
      complex=total_complex[[1]]
      ec_curve <- matrix(NA,nrow = 1,ncol=0)
      for (j in 1:dim(dir)[1]){
        vertex_function <- complex$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
        if (ball == TRUE){
          curve <- compute_standardized_ec_curve(complex, vertex_function, curve_length-1, first_column_index = FALSE,ball_radius)
        }
        else{
          curve = compute_discrete_ec_curve(complex, vertex_function , curve_length-1, first_column_index = FALSE)
        }

        ### transform the ECT as desired ###
        if (ec_type == "SECT"){
          curve <- integrate_ec_curve(curve)
        } else if(ec_type == "DECT"){
          curve <- differentiate_ec_curve(curve)
        } else{
          curve <- curve
        }
        # omit the length data, for now
        ec_curve <- c(ec_curve,curve[,2])
      }

      data <- rbind(data,c(1,ec_curve))
      complex=total_complex[[3]]
      ec_curve <- matrix(NA,nrow = 1,ncol=0)
      for (j in 1:dim(dir)[1]){
        vertex_function <- complex$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
        curve <- compute_standardized_ec_curve(complex, vertex_function, curve_length-1, first_column_index = FALSE,ball_radius)

        if (ec_type == "SECT"){
          curve <- integrate_ec_curve(curve)
        } else if(ec_type == "DECT"){
          curve <- differentiate_ec_curve(curve)
        } else {
          curve <- curve
        }

        # omit the length data, for now
        ec_curve <- c(ec_curve,curve[,2])
      }
      data <- rbind(data,c(-1,ec_curve))
    }
  }
  data_list=list(data=data,noise=noise_points,causal_points1=causal_points1,causal_points2=causal_points2, complex_points = complex_points)
  return(data_list)
}

#'Generate Kernel interpolated meshes.
#' @export
#'
#' @description Given a set of input data of shared points and causal points,
#'  the function generates two classes of interpolated 3d shapes with a specified kernel.
#'
#' @param grid_size (int) : The fineness of the grid for which the interpolation is on.
#' @param noise_points (matrix) : Matrix of points to be shared across both classees of shapes.
#' @param causal_points1 (matrix) : Matrix of points unique to class 1.
#' @param causal_points2 (matrix) : Matrix of points unique to class 2.
#' @param func (function)  : The function used for the kernel interpolation.
#' @param eta (float)  : The shape paramter for the kernel.
#'
#' @return complex (list) : List of data containing the mesh versions of the interpolated shapes, causal points, and the shared points.
generate_normal_complex_fixed=function(grid_size=25,noise_points,causal_points1,causal_points2,func=rbf_gauss,eta=5){
  noise=cbind(noise_points,rnorm(dim(noise_points)[1],1,0.25))
  samples1=cbind(causal_points1,rnorm(dim(causal_points1)[1],1,0.25))
  samples2=cbind(causal_points2,rnorm(dim(causal_points2)[1],1,0.25))

  real_samples1=rbind(samples1,noise)
  real_samples2=rbind(samples2,noise)

  predictions1=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=real_samples1,eta=eta)
  predictions2=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=real_samples2,eta=eta)

  complex1=matrix_to_simplicial_complex(predictions1,grid_length=grid_size)
  complex2=matrix_to_simplicial_complex(predictions2,grid_length=grid_size)
  complex=list(complex1=complex1,base_points1=samples1,complex2=complex2,base_points2=samples2,
               noise=noise, class1_points = real_samples1, class2_points = real_samples2)
}


### Code for putting gaussian points on random fields###
# generate these shapes and their ec curves
generate_data_gaussian_field <- function(nsim, curve_length, dir, shared_points = 2, causal_points = 1, ball_radius = 2,
                                         ec_type = "ECT", grid_size = 25){

  dir <- matrix(dir, ncol = 3)
  data <- matrix(NA,nrow=0,ncol = 1 + curve_length * dim(dir)[1] )


  #Shared points
  n1=runif(shared_points,-1,1)
  n2=runif(shared_points,-1,1)

  #causal points
  x1=runif(causal_points,-1,1)
  y1=runif(causal_points,-1,1)
  x2=runif(causal_points,-1,1)
  y2=runif(causal_points,-1,1)

  noise_points=cbind(n1,n2)
  causal_points1=cbind(x1,y1)
  causal_points2=cbind(x2,y2)
  complex_points=list()
  for (i in 1:nsim){
    total_complex=generate_gaussian_field(grid_size = grid_size, shared_points = noise_points, causal_points1 = causal_points1,
                                                causal_points2 = causal_points2, 0.01)
    if(inherits(total_complex,'try-error')){
      i=i-1
      next
    }
    else{
      complex_points[[(2*i-1)]] = total_complex[[6]]
      complex_points[[2*i]] = total_complex[[7]]
      complex=total_complex[[1]]
      ec_curve <- matrix(NA,nrow = 1,ncol=0)
      for (j in 1:dim(dir)[1]){
        vertex_function <- complex$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
        curve <- compute_standardized_ec_curve(complex, vertex_function, curve_length-1, first_column_index = FALSE,ball_radius)

        ### transform the ECT as desired ###
        if (ec_type == "SECT"){
          curve <- integrate_ec_curve(curve)
        } else if(ec_type == "DECT"){
          curve <- differentiate_ec_curve(curve)
        } else{
          curve <- curve
        }
        # omit the length data, for now
        ec_curve <- c(ec_curve,curve[,2])
      }

      data <- rbind(data,c(1,ec_curve))
      complex=total_complex[[3]]
      ec_curve <- matrix(NA,nrow = 1,ncol=0)
      for (j in 1:dim(dir)[1]){
        vertex_function <- complex$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
        curve <- compute_standardized_ec_curve(complex, vertex_function, curve_length-1, first_column_index = FALSE,ball_radius)

        if (ec_type == "SECT"){
          curve <- integrate_ec_curve(curve)
        } else if(ec_type == "DECT"){
          curve <- differentiate_ec_curve(curve)
        } else {
          curve <- curve
        }

        # omit the length data, for now
        ec_curve <- c(ec_curve,curve[,2])
      }
      data <- rbind(data,c(-1,ec_curve))
    }
  }
  data_list=list(data=data,noise=noise_points,causal_points1=causal_points1,causal_points2=causal_points2, complex_points = complex_points)
  return(data_list)

}

#' Place Gaussian Bumps on a Plane
#' @export
#'
#'
#' @description Place gaussian bumps on a plane at desired locations.
#'
#' @param shared_points (matrix) :rows are shared points between the classes.
#' @param causal_points1 (matrix) : should be a matrix whose rows are the causal points of class 1
#' @param causal_points2 (matrix) : should be a matrix whose rows are the causal points of class 2
#' @param grid_size (int): the fineness of the grid.
#'
#' @return complex (list) : Mesh form of the grid.
generate_gaussian_field <- function(grid_size = 25, shared_points, causal_points1, causal_points2,
                                    normal_std){

  class1_points <- rbind(causal_points1, shared_points)
  class2_points <- rbind(causal_points2, shared_points)

  class1_field <- generate_random_field_matrix(grid_length = grid_size, class1_points,
                                               normal_std)
  class2_field <- generate_random_field_matrix(grid_length = grid_size, class1_points,
                                               normal_std)

  complex1 <- matrix_to_simplicial_complex(class1_field, grid_length = grid_size)
  complex2 <- matrix_to_simplicial_complex(class2_field, grid_length = grid_size)

  complex <- list(complex1 = complex1, base_points1 = causal_points1,
                  complex2 = complex2, base_points2 = causal_points2,
                  shared_points = shared_points,
                  class1_points = class1_points,
                  class2_points = class2_points)
}


######################################################################################
######################################################################################
######################################################################################

#' Convert 2d function to simplicial complex
#' @export
#' @import Rvcg
#' @description Using Rvcg, we turn a matrix into a triangular simplicial complex.
#'
#' @param matrix (nx3 matrix) The matrix we want to convert into a simplicial complex
#' @param grid_length (float): The length of the grid for which the shape is defined on.
#'
#' @return complex (list) : List of metadata about our mesh; the vertex coordinates, edge positions, and face positions.
matrix_to_simplicial_complex <- function(matrix,grid_length){
  vertices <- matrix(NA,nrow = 0, ncol = 4)
  edges <- matrix(NA, nrow = 0, ncol = 2)
  faces <- matrix(NA,nrow = 0, ncol = 3)
  length_x <- dim(matrix)[1]
  length_y <- dim(matrix)[2]

  for(i in 1:(length_x - 1)){
    for(j in 1:(length_y - 1)){
      #The last three columns are the coordinates of the vertices
      #indices of the vertices of this pixel:
      vertex1 <- (i - 1)*(dim(matrix)[2]) + j
      vertex2 <- (i - 1)*(dim(matrix)[2]) + j+1
      vertex3 <- (i)*(dim(matrix)[2]) + j
      vertex4 <- (i)*(dim(matrix)[2]) + j+1

      # add the vertices to complex
      vertices <- rbind(vertices, c( vertex1, i-length_x/2, j-length_y/2, matrix[i,j]) )
      vertices <- rbind(vertices, c( vertex2, i-length_x/2, j+1-length_y/2, matrix[i,j+1]))
      vertices <- rbind(vertices, c( vertex3, i+1-length_x/2, j-length_y/2, matrix[i+1,j] ))
      vertices <- rbind(vertices, c( vertex4, i+1-length_x/2, j+1-length_y/2, matrix[i+1,j+1] ))

      # add the edges to complex
      edges <- rbind(edges, c(vertex1,vertex2))
      edges <- rbind(edges, c(vertex2,vertex4))
      edges <- rbind(edges, c(vertex1,vertex3))
      edges <- rbind(edges, c(vertex3,vertex4))
      #edges<-rbind(edges,c(vertex2,vertex3))
      # add the faces to complex
      #faces <- rbind(faces, c(vertex1,vertex2,vertex3,vertex4))
      faces <- rbind(faces, c(vertex1,vertex2,vertex3))
      faces <- rbind(faces, c(vertex2,vertex3,vertex4))
    }
  }

  vertices <- unique(vertices)
  vertices = vertices[order(vertices[,1]),]
  edges <- unique(edges)
  faces <- unique(faces)
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  complex$Vertices[,2]=complex$Vertices[,2]/(grid_length/2)
  complex$Vertices[,3]=complex$Vertices[,3]/(grid_length/2)
  complex$Vertices=complex$Vertices[,2:4]
  vertex=complex$Vertices
  faces=complex$Faces
  ind=rep(1,dim(vertex)[1])
  vertex=cbind(vertex,ind)
  mesh=tmesh3d(t(vertex),indices=t(faces))
  vertices=as.matrix(t(mesh$vb)[,1:3])
  faces=as.matrix(t(mesh$it))
  edges=Rvcg::vcgGetEdge(mesh)
  edges=as.matrix(edges[,1:2])
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  return(complex)
}
######################################################################################
######################################################################################
######################################################################################
### Radial Basis Function Generation ###

#' Euclidean Distance
#' @export
#'
#' @description Computes the Euclidean distance between two vectors
#'
#' @param x (vector) : The first vector
#' @param y (vector) : The second vector
#'
#' @return The Euclidean distance
difference=function(x,y){
  sumxy=sum((x-y)^2)
  return(sqrt(sumxy))
}
#' RBF Kernel
#' @export
#'
#' @description Computes the squared euclidean distance between two vectors
#'
#' @param x (vector) The first vector
#' @param y (vector) The second vector
#' @param eta (int) The varying weight parameter (not used here)
#'
#' @return The squared euclidean distance
rbf_gauss=function(x,y,eta=5){
  return(exp(-(difference(x,y))))
}

#' Inverse Quadratic Kernel
#' @export
#'
#' @description Computes the inverse quadratic distance between two vectors
#'
#' @param x (vector) The first vector
#' @param y (vector) The second vector
#' @param eta (float) The varying shape parameter.
#'
#' @return The inverse quadratic distance
inv_quad=function(x,y,eta=5){
  return(1/(sqrt(1+(difference(x,y)/eta)^2)))
}

#' Computes the RBF Kernel
#' @export
#'
#' @description Computes the inverse quadratic distance between two vectors
#'
#' @param x (vector) The first vector
#' @param y (vector) The second vector
#' @param eta (float) The varying shape parameter.
#'
#' @return rbf_final (list) : The rbf model values: rbf weights, input data, and interpolation function
compute_rbf_model=function(data,eta=5,func){
  x=data[,1:2]
  y=data[,3]
  samples=dim(x)[1]
  kernel_matrix=matrix(NA,ncol=samples,nrow=samples)
  for (i in 1:samples){
    kernel_matrix[i,]=apply(X=x,FUN=func,MARGIN=1,y=x[i,])
  }
  lambdas=solve(kernel_matrix)%*%y
  rbf_final=list(lambdas=lambdas,x=x,values=y,eta=eta,func=func)
}

#' RBF prediction
#' @export
#'
#' @description Computes the RBF kernel interpolation on a provided grid.
#'
#' @param rbf_model (list): The RBF kernel model containing metadata bout the model (kernels, rbf weights, interpolation function).
#' @param grid (nxn matrix): Grid of points to compute the RBF kernel on.
#'
#' @return predictions (matrix) : The rbf kernel interpolated values.
rbf_predict=function(rbf_model,grid){
  lambdas=rbf_model$lambdas
  x=rbf_model$x
  func=rbf_model$func
  num_preds=dim(grid)[1]
  num_centers=dim(x)[1]
  predictions=rep(0,num_preds)
  for (i in 1:num_centers){
    distances=apply(X = grid,FUN = func,MARGIN=1,y=x[i,])
    distances=distances*lambdas[i]
    predictions=predictions+distances
  }
  return(predictions)
}

#' Wrapper for RBF Kernel Interpolation
#' @export
#'
#' @description Provides Kernel Interpolated predictions.
#'
#' @param grid_size (int) The granularity of the grid.
#' @param func (function) The interpolation function used. RBF and Inverse quardrtic kernels are supported.
#' @param data (nx3 matrix) The data used to interpolate.
#' @param eta (float) Varying shape parameter of the kernel
#'
#' @return predictions (matrix) Matrix of Kernel predictions at each point of the grid.
rbf_on_grid=function(grid_size=25,func=rbf_gauss,data,eta=5){
  rbf_model=compute_rbf_model(data=data,eta=eta,func=func)
  grids=seq(-1,1,length=grid_size)
  grid=expand.grid(grids,grids)
  predictions=rbf_predict(rbf_model = rbf_model,grid=grid)
  predictions=matrix(predictions,nrow=grid_size)
  return(predictions)
}

# maps projective coordinates to euclidean coordinates
asEuclidean <- function(coordinates){
  return(coordinates[,1:3]/coordinates[,4])
}

######################################################################################
######################################################################################
######################################################################################
### Helper functions for Gaussian on plane functions ###

# Generates the values over a grid for a random field with
generate_random_field_matrix <- function(grid_length, points, normal_std){
  grid_ticks <- seq(-1, 1, length = grid_length)
  grid <- expand.grid(grid_ticks, grid_ticks)

  field_values <- function(x) sum_normals(x,points, normal_std)

  function_values <- matrix(apply(grid, 1, field_values), nrow = grid_length, byrow = TRUE)

  # add some noise to the grid values
  function_values + matrix(runif(grid_length^2, min=-0.1, max=0.1), nrow = grid_length)
}

# Generates a function that is the sum of normals centered at inputted points, with
# specified standard deviation
sum_normals <- function(test_point, points, normal_std){
  functions <- c()
  for(i in 1:(dim(points)[1])){
    functions[i] <- dmvnorm(test_point,points[i,],sigma = diag(2)*normal_std)
  }

  sum(functions)/10 #for renormalization
}

