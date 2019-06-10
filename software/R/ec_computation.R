#### Computing the (S/D) ECT of a Simplicial Complex ####

#' Computes the EC curve on a ball
#'
#' @export
#' @description Computes the EC curve of a mesh on a bounding ball of specified radius $radius$. The filtration steps for which the
#' (S/D) EC curve are computed is done relative to the bounding ball. For filtrations steps that are computed relative to the shape, use
#' \code{compute_discrete_ec_curve}. For comparisons with multiple shapes, this is the recommended function.
#'
#' @param complex (list) : The list containing metadata about the mesh; use \code{process_off_filev3} to obtain this list from an off file, or
#' \code{convert_complex} to convert existing mesh3d objects to this list.
#' @param vertex_function (matrix): A matrix containing the projections of each vertex along a given direction. Computed as the dot product of the vertex onto the
#' direction on the unit sphere.
#' @param curve_length (int): The number of sub-level sets for which to compute the EC curve on in a given direction.
#' @param first_column_index (boolean): Specifying the vertex index is included in the vertex function matrix.
#' @param ball_radius (float): The radius of the bounding ball.
#'
#' @return ec_curve (nx2 matrix): A matrix containing the EC curve of the given simplicial complex with the first index as the projections for which the EC was computed on.
compute_standardized_ec_curve <- function(complex, vertex_function, curve_length, first_column_index = FALSE, ball_radius){

  V = complex$Vertices
  E = complex$Edges
  F = complex$Faces


  if ( dim(V)[1] != dim(vertex_function)[1] ) {
    print('The size of function should be same as the number of vertices')
    return(0)
  }

  # In the case where the vertex_indices are given in the first column. We pick out the birth times of each vertex in every
  # edge and face of the complex.
  if (first_column_index == TRUE){

    edge_birth_times = sapply(1:dim(E)[1],function(ind) max(vertex_function[ which(vertex_function[,1]%in%E[ind,]), 2 ] ))
    face_birth_times = sapply(1:dim(F)[1],function(ind) max(vertex_function[ which(vertex_function[,1]%in%F[ind,]), 2 ] ))
  }else{

    #Apply a function on indices on this vector? Count the number of edges that are present given
    # the value of the filtration (get max vertex value at each edge, face).
    edge_birth_times = sapply(1:dim(E)[1],function(ind) max(vertex_function[ E[ind,] ]))
    face_birth_times = sapply(1:dim(F)[1],function(ind) max(vertex_function[ F[ind,] ]))
  }

  if (first_column_index == FALSE){
    threshold = seq(from=-ball_radius,to=ball_radius,length =curve_length+1)
  } else{
    threshold = seq(from=min(vertex_function[,2])-buffer,to=max(vertex_function[,2])+buffer,length =curve_length+1)
  }


  ec_curve = matrix(0,curve_length+1,2)
  ec_curve[,1] <- threshold

  # count how many of each object is born given the threshold time.
  for ( i in 1:length(threshold) ) {
    v = length(which(vertex_function <= threshold[i]))
    e = length(which(edge_birth_times <= threshold[i]))
    f = length(which(face_birth_times <= threshold[i]))
    ec_curve[i,2] <- v - e + f;
  }
  return(ec_curve)
}

#' Computes the EC curve relative to the shape
#' @export
#'
#' @description Computes the EC curve of a mesh on relative to the shape. The filtration steps for which the
#' (S/D) EC curve are computed is done relative to the shape. For filtrations steps that are computed relative to a bounding ball, use
#' \code{compute_standardized_ec_curve}. For comparisons with multiple shapes, this is not the recommended function.
#'
#' @param complex (list) : The list containing metadata about the mesh; use \code{process_off_filev3} to obtain this list from an off file, or
#' \code{convert_complex} to convert existing mesh3d objects to this list.
#' @param vertex_function (matrix): A matrix containing the projections of each vertex along a given direction. Computed as the dot product of the vertex onto the
#' direction on the unit sphere.
#' @param curve_length (int): The number of sub-level sets for which to compute the EC curve on in a given direction.
#' @param first_column_index (boolean): Specifying the vertex index is included in the vertex function matrix.
#'
#' @return ec_curve (nx2 matrix): A matrix containing the EC curve of the given simplicial complex with the first index as the projections for which the EC was computed on.
compute_discrete_ec_curve <- function(complex, vertex_function, curve_length, first_column_index = FALSE){

  V = complex$Vertices
  E = complex$Edges
  F = complex$Faces
  buffer=1/(2*curve_length+5)

  if ( dim(V)[1] != dim(vertex_function)[1] ) {
    print('The size of function should be same as the number of vertices')
    return(0)
  }

  # In the case where the vertex_indices are given in the first column. We pick out the birth times of each vertex in every
  # edge and face of the complex.
  if (first_column_index == TRUE){

    edge_birth_times = sapply(1:dim(E)[1],function(ind) max(vertex_function[ which(vertex_function[,1]%in%E[ind,]), 2 ] ))
    face_birth_times = sapply(1:dim(F)[1],function(ind) max(vertex_function[ which(vertex_function[,1]%in%F[ind,]), 2 ] ))
  }else{

    #Apply a function on indices on this vector? Count the number of edges that are present given
    # the value of the filtration (get max vertex value at each edge, face).
    edge_birth_times = sapply(1:dim(E)[1],function(ind) max(vertex_function[ E[ind,] ]))
    face_birth_times = sapply(1:dim(F)[1],function(ind) max(vertex_function[ F[ind,] ]))
  }

  if (first_column_index == FALSE){
    threshold = seq(from=min(vertex_function)-buffer,to=max(vertex_function)+buffer,length =curve_length+1)
  } else{
    threshold = seq(from=min(vertex_function[,2])-buffer,to=max(vertex_function[,2])+buffer,length =curve_length+1)
  }


  ec_curve = matrix(0,curve_length+1,2)
  ec_curve[,1] <- threshold

  # count how many of each object is born given the threshold time.
  for ( i in 1:length(threshold) ) {
    v = length(which(vertex_function <= threshold[i]))
    e = length(which(edge_birth_times <= threshold[i]))
    f = length(which(face_birth_times <= threshold[i]))
    ec_curve[i,2] <- v - e + f;
  }
  return(ec_curve)
}

#' Differentiates an EC curve
#' @export
#'
#' @description Computes the numerical derivative of an Euler Characteristic Curve curve using successive differences.
#' @param ec_curve (matrix) : The Euler Characteristic Curve from \code{compute_standardized_ec_curve} or \code{compute_discrete_ec_curve}.
#'
#' @return ec_curve (matrix): The differentiated Euler Characteristic Curve (DECT)
differentiate_ec_curve <- function(ec_curve){
  differences <- diff(ec_curve[,2],lag = 1)
  ec_curve[,2] <- c(ec_curve[1,2], differences)

  ec_curve
}

#' Integrates an EC curve
#' @export
#'
#' @description Computes the numerical integral of an Euler Characteristic Curve curve by computinig the area under the EC curve.
#' @param ec_curve (matrix) : The Euler Characteristic Curve from \code{compute_standardized_ec_curve} or \code{compute_discrete_ec_curve}.
#'
#' @return ec_curve (matrix): The Integrated Euler Characteristic Curve (SECT)
integrate_ec_curve <- function(ec_curve){
  length <- length(ec_curve[,2])
  ec_curve[,2] <- ec_curve[,2] - mean(ec_curve[,2])
  ec_curve[,2] <- cumsum(ec_curve[,2])* (ec_curve[length,1] - ec_curve[1,1])/length

  ec_curve
}

#' Updates the EC curve
#' @export
#'
#' @description Transforms the Euler Characteristic Curve into our desired form.
#' @param ec_curve (matrix) : The Euler Characteristic Curve from \code{compute_standardized_ec_curve} or \code{compute_discrete_ec_curve}.
#' @param ec_type (string): The type of EC curve we want. Currently the SECT, DECT and ECT are supported.
#'
#' @return ec_curve (matrix): The (S/D) Euler Characteristic Curve.
update_ec_curve = function(curve,ec_type){
  if (ec_type == "SECT"){
    curve <- integrate_ec_curve(curve)
  } else if(ec_type == "DECT"){
    curve <- differentiate_ec_curve(curve)
  } else{
    curve <- curve
  }
  curve
}
