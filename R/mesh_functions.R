#### Functions to help work with real meshes ####

#' Read in OFF files to list form
#' @export
#' @import Rvcg
#' @description  Processing OFF files in a format friendly to our EC curve computations.
#'
#' @param input_dir (string): The input directory of our off file.
#' @return complex (list): A list containing the Vertex coordinates, edge and face information.
process_off_file_v3=function(input_dir){
  off=vcgImport(input_dir,silent = TRUE)
  vertices=as.matrix(t(off$vb)[,1:3])
  faces=as.matrix(t(off$it))
  edges=vcgGetEdge(off)
  edges=as.matrix(edges[,1:2])
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  return(complex)
}

#### Processing OFF files from a Directory

#' Create EC matrix from meshes.
#' @export
#' @import stringr
#' @description  Create an EC curve matrix given an input directory and set of directions to compute the (S/D) EC curve over.
#' Each row corresponds to one mesh, and the columns correspond to the sub-level sets of the directions in the matrix.
#'
#'@param directory (string): :  The input directory of the meshes.
#'@param directions (nx3 matrix): :  The matrix of directions to compute the (S/D) EC curve over.
#'@param len (int): The number of sub-level sets to compute the (S/D) EC curve on in each direction.
#'@param ball_radius (float): The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.
#'@param ball (boolean): Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.
#'@param ec_type (string): What type of EC curve to compute. Currently we support ECT (Euler Characteristic Curve), SECT (Smooth Euler Characteristic Curve)
#' and DECT (Differentiated Euler Characteristic Curve). We use ECT in the papers.
#'
#' @return data (matrix): The matrix of (S/D) EC curves of the meshes in the input directory.
create_ec_matrix_mult_d=function(directory,directions,len,ball_radius = 1, ball = FALSE,ec_type = 'ECT'){
  curve_length = len
  file_names=list.files(path=directory,full.names = TRUE)
  file_names = file_names[stringr::str_detect(file_names,'off')]
  number_files=length(file_names)
  nvertices =  length(mesh_to_matrix(Rvcg::vcgImport(file_names[1])))
  if (ec_type == 'baseline'){
    data <- matrix(NA,nrow=number_files,ncol = nvertices)
  }
  else{
    data <- matrix(NA,nrow=number_files,ncol = curve_length*dim(directions)[1])
  }
  for (i in 1:number_files){
    print(paste('On File', i))
    off=process_off_file_v3(file_names[i])
    if (ec_type == 'baseline'){
      mesh = Rvcg::vcgImport(file_names[i])
      curve = mesh_to_matrix(mesh)
      data[i,] <- curve
    }
    else{
      curve_mult_d <- matrix(NA,nrow = 1,ncol=0)
      if (ball == FALSE){
        for (j in 1:dim(directions)[1]){
          vertex_function=off$Vertices%*%directions[j,]
          curve <- compute_discrete_ec_curve(off, vertex_function, len-1, first_column_index = FALSE)
          if (ec_type == 'ECT'){
            curve = curve
          }
          if (ec_type == 'SECT'){
            curve <- integrate_ec_curve(curve)
          }
          curve_mult_d <- cbind(curve_mult_d,t(curve[,2]))
        }
        data[i,] <- curve_mult_d
      }
      if (ball == TRUE){
        for (j in 1:dim(directions)[1]){
          vertex_function=off$Vertices%*%directions[j,]
          curve <- compute_standardized_ec_curve(off, vertex_function, len-1, first_column_index = FALSE,ball_radius = ball_radius)
          if (ec_type == 'ECT'){
            curve = curve
          }
          if (ec_type == 'SECT'){
            curve <- integrate_ec_curve(curve)
          }
          if (ec_type == 'DECT'){
            curve = differentiate_ec_curve(curve)
          }
          curve_mult_d <- cbind(curve_mult_d,t(curve[,2]))
          # omit the length data, for now
        }
        data[i,] <- curve_mult_d
      }
    }
  }
  print(dim(data))
  return(data)
}


#' Generate design matrix from meshes in two directories.
#'  @description Given two classes & directories, we generate the (S/D) EC curves and the associated class labels
#' @export
#'
#' @param directory1 (string): The first directory corresponding to class 1.
#' @param directory2 (string): The second directory correponding to class 2.
#'@param directions (nx3 matrix): The matrix of directions to compute the (S/D) EC curve over.
#'@param len (int): The number of sub-level sets to compute the (S/D) EC curve on in each direction.
#'@param ball (boolean): Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.
#'@param ball_radius (float): The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.
#'@param ec_type (string): What type of EC curve to compute. Currently we support ECT (Euler Characteristic Curve), SECT (Smooth Euler Characteristic Curve)
#' and DECT (Differentiated Euler Characteristic Curve). We use ECT in the papers.
#'
#' @return final_matrix (matrix): The matrix containing the EC curves of the meshes from the two classes and the associated class labels.
create_comparison_matrix_mult_d=function(directory1,directory2,directions,len, ball = FALSE, ball_radius = 1, ec_type = 'ECT'){
  matrix_1=create_ec_matrix_mult_d(directory1,directions,len,ball_radius = ball_radius, ball = ball, ec_type = ec_type)
  matrix_2=create_ec_matrix_mult_d(directory2,directions,len,ball_radius = ball_radius, ball = ball, ec_type = ec_type)
  matrix_1=cbind(rep(1,dim(matrix_1)[1]),matrix_1)
  matrix_2=cbind(rep(0,dim(matrix_2)[1]),matrix_2)
  final_matrix=rbind(matrix_1,matrix_2)
  return(final_matrix)
}

### find landmarks on the originl tooth
get_euclidean_fps_landmarks = function(mesh, num_landmarks){
  ### Gets landmarks using farthest point sampling
  landmarks = rdist::farthest_point_sampling(t(mesh$vb)[,1:3], metric = 'euclidean', num_landmarks)
  landmarks
}

create_landmark_comparison_matrix=function(directory1,directory2,base_shape_dir, num_landmarks){

  base_shape = Rvcg::vcgImport(base_shape_dir)
  landmark_indices = get_euclidean_fps_landmarks(base_shape, num_landmarks)
  #browser()
  matrix_1=create_landmark_summary(directory1, landmark_indices)
  matrix_2=create_landmark_summary(directory2, landmark_indices)

  matrix_1=cbind(rep(1,dim(matrix_1)[1]),matrix_1)
  matrix_2=cbind(rep(0,dim(matrix_2)[1]),matrix_2)
  final_matrix=rbind(matrix_1,matrix_2)
  return(comparison_matrix = final_matrix)
}

create_landmark_summary=function(directory, landmark_indices){
  file_names = list.files(path=directory,full.names = TRUE)
  file_names = file_names[stringr::str_detect(file_names,'off')]
  number_files=length(file_names)

  data <- matrix(NA,nrow=number_files,ncol = 3*length(landmark_indices))
  #browser()
  for (i in 1:number_files){
    mesh = Rvcg::vcgImport(file_names[i])
    curve = matrix(t(t(mesh$vb)[landmark_indices,1:3]), nrow = 1)
    # mfrow3d(1, 1, byrow = TRUE)
    # rgl::plot3d(t(mesh$vb)[landmark_indices,1:3])
    data[i,] <- curve
  }

  print(dim(data))
  return(data)
}


#### Summary Function for the Real Data ####

#' Creates Design Matrix and Run Rate
#'
#' @export
#' @description  Creates the EC comparison matrix for two classes of meshes, and runs variable selection on the two class matrix.
#' We currently fit a Gaussian Process Classifier and use RATE for variable selection.
#' @param dir1 (string): The first directory corresponding to class 1.
#' @param dir2 (string): The second directory correponding to class 2.
#'@param len (int): The number of sub-level sets to compute the (S/D) EC curve on in each direction.
#'@param ec_type (string): What type of EC curve to compute. Currently we support ECT (Euler Characteristic Curve), SECT (Smooth Euler Characteristic Curve)
#' and DECT (Differentiated Euler Characteristic Curve). We use ECT in the papers.
#'@param ball (boolean): Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.
#'@param ball_radius (float): The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.
#'
#'@return final_list (list): A list containing the design matrix of meshes from the two classes and the importances of each variable.
real_data_summary = function(shape_transformation, dir1,dir2,base_shape_dir,direction=dir,
                             len=len,ec_type = 'ECT',
                             ball = TRUE, ball_radius = 1.5,
                             mode = 'rate', reduce = max, alpha = 0.5, num_landmarks = 100){
  #Generate Matrix of EC curves and labels
  if (shape_transformation == 'ECT') {
    comparison_matrix = create_comparison_matrix_mult_d(dir1,dir2,direction,len,
                                                        ec_type = ec_type, ball = ball, ball_radius = ball_radius)
  } else if (shape_transformation == 'landmark') {
    comparison_matrix = create_landmark_comparison_matrix(dir1,dir2, base_shape_dir, num_landmarks)
  } else {
    print('shape transformation specified must be ECT or landmark')
  }

  #Feature Selection
  if (mode == 'rate'){
    want_indices=find_rate_variables_with_other_sampling_methods(comparison_matrix,bandwidth = 0.01, type = 'ESS')
    want_indices = want_indices[,2]
  }
  if (mode == 'elastic_net'){
    print(dir1)
    print(dir2)
    print(dim(comparison_matrix))
    rate_values = abs(find_elastic_variables(comparison_matrix,weights = TRUE, grouped_data = TRUE, alpha = alpha))
    print(dim(rate_values))
    rate_values[,1] = rep((1:(dim(rate_values)[1]/3)),each = 3)
    #browser()
    new_df = stats::aggregate(rate_values[,2],list(rate_values[,1]),reduce)
    want_indices = new_df$x

  }

  final_list=list(data = comparison_matrix, Rate2=want_indices, inference_result = rate_values)
  return(final_list)
}
