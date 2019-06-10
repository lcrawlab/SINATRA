#### Shape Reconstruction ####

#' Reconstruct Vertices
#'
#' @export
#' @description Given a set of variable importances, threshold, the function reconstructs the vertices above the threshold
#' by intersecting the sub-level sets of the directions in a cone with variable importance greater than the threshold.
#' This function loops over each cone, and computes the reconstructed vertices using intersections of sub-level sets. The vertices
#' from each cone are then unioned in the end for the final set of reconstructed vertices.
#'
#'@param dir (nx3 matrix) The matrix of directions that were used compute the (S/D) EC curve over.
#'@param complex (list) : The list containing metadata of the Vertices, Edges, and Faces of the mesh (use process_off_file_v3 to obtain).
#'@param rate_vals (vector) : Vector of variable importances for each sub-level set across each direction
#'@param len (int) : The number of sub-level sets to compute the (S/D) EC curve on in each direction.
#'@param threshold (float) : The threshold for determining which sub-level sets are used for the reconstruction.
#'@param cone_size (int) : The number of directions in each cone.
#'@param ball_radius (float) : The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.
#'@param ball (boolean) : Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.
#'@param ec_type (string) : What type of EC curve to compute. Currently we support ECT (Euler Characteristic Curve), SECT (Smooth Euler Characteristic Curve)
#'@param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include.
#'Setting Radius = 1 is recommened.
#'
#'@return total_selected_vertices (vector) : Vector of the vertex indices that are reconstructed.

compute_selected_vertices_cones = function(dir, complex, rate_vals, len, threshold=-1, cone_size, ball = TRUE, ball_radius,radius = 1){
  if (threshold==-1){
    threshold=1/length(rate_vals)
  }
  if ((dim(dir)[1] %% cone_size) != 0){
    print('Number of Cones not a multiple of directions')
    return(0)
  }
  coned_vertices=list()
  for (j in 1:(dim(dir)[1] / cone_size)){
    cone_dirs=matrix(dir[((j-1)*(cone_size)+1):(j*cone_size),],ncol = 3)
    cone_rate_vals=rate_vals[(j-1)*(cone_size*len)+1:(j*cone_size*len)]
    coned_vertices[[j]]=summarize_vertices(dir = cone_dirs, complex, rate_vals = cone_rate_vals, len,
                                           reduction_operation = intersect, threshold, cone_size, ball = ball, ball_radius, radius = radius)
  }
  total_selected_vertices=Reduce(union,coned_vertices)
  return(total_selected_vertices)
}
#' Reconstruct Faces
#' @description Given a set of variable importances, threshold, the function reconstructs the faces above the threshold
#' by intersecting the sub-level sets of the directions in a cone with variable importance greater than the threshold.
#'
#' @export
#'@param dir (nx3 matrix) The matrix of directions that were used to compute the (S/D) EC curve over.
#'@param complex (list) : The list containing metadata of the Vertices, Edges, and Faces of the mesh (use process_off_file_v3 to obtain).
#'@param rate_vals (vector) : Vector of variable importances for each sub-level set across each direction
#'@param len (int) : The number of sub-level sets to compute the (S/D) EC curve on in each direction.
#'@param threshold (float) : The threshold for determining which sub-level sets are used for the reconstruction.
#'@param cone_size (int) : The number of directions in each cone.
#'@param ball_radius (float) : The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.
#'@param ball (boolean) : Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.
#'@param ec_type (string) : What type of EC curve to compute. Currently we support ECT (Euler Characteristic Curve), SECT (Smooth Euler Characteristic Curve)
#'@param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include.
#'Setting Radius = 1 is recommened.
#'
#'@return total_selected_faces (vector) : Vector of the face indices that are reconstructed.

compute_selected_faces_cones = function(dir, complex, rate_vals, len, threshold=-1, cone_size, ball = TRUE, ball_radius,radius = 0){
  if (threshold==-1){
    threshold=1/length(rate_vals)
  }
  if ((dim(dir)[1] %% cone_size) != 0){
    print('Number of Cones not a multiple of directions')
    return(0)
  }
  coned_vertices=list()
  for (j in 1:(dim(dir)[1] / cone_size)){
    cone_dirs=matrix(dir[((j-1)*(cone_size)+1):(j*cone_size),],ncol = 3)
    cone_rate_vals=rate_vals[(j-1)*(cone_size*len)+1:(j*cone_size*len)]
    coned_vertices[[j]]=summarize_vertices(dir = cone_dirs, complex, rate_vals = cone_rate_vals, len,
                                           reduction_operation = intersect, threshold, cone_size, ball = ball, ball_radius, radius = radius)
  }
  total_selected_vertices=Reduce(union,coned_vertices)
  reconstructed_faces = apply(X = complex$Faces,MARGIN = 1,function(x) any(x %in% total_selected_vertices))
  reconstructed_faces = which(reconstructed_faces == TRUE)
  return(reconstructed_faces)
}

#' Reconstruct Vertices in a Cone
#'
#' @description Given a cone of directions and a threshold, We take the sub-level sets that are above the given threshold, and for each direction in the cone,
#' we compute the set of vertices that are associated with the selected sub-level sets. We then intersect the set of
#' vertices that are associated with the sub-level sets across each direction in the cone; only taking the vertices that are reconstructed in each direction of the cone.
#'
#' @export
#'@param dir (nx3 matrix) : The matrix of directions in the cone that were used to compute the (S/D) EC curve over.
#'@param complex (list) : The list containing metadata of the Vertices, Edges, and Faces of the mesh (use process_off_file_v3 to obtain).
#'@param rate_vals (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#'@param len (int) : The number of sub-level sets to compute the (S/D) EC curve on in each direction.
#'@param reduction_operation (function) The function to summarize the sets of vertices. We default to intersect, and recommend it.
#'@param threshold (float) : The threshold for determining which sub-level sets are used for the reconstruction.
#'@param cone_size (int) : The number of directions in each cone.
#'@param ball_radius (float) : The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.
#'@param ball (boolean) : Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.
#'@param ec_type (string) : What type of EC curve to compute. Currently we support ECT (Euler Characteristic Curve), SECT (Smooth Euler Characteristic Curve)
#'@param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include.
#'Setting Radius = 1 is recommened.
#'
#'@return final_selected_vertices (vector) : Vector of the vertex indices that are reconstructed for a given cone.

summarize_vertices=function(dir,complex,rate_vals,len,reduction_operation=intersect,threshold,cone_size, ball = TRUE, ball_radius = 1, radius = 0){
  picked_indices=which(rate_vals>=threshold)
  indices=c()
  for (j in 0:radius){
      indices=c(indices,picked_indices+j)
      indices=c(indices,picked_indices-j)
  }
  selected_vertices=list()

  # Count how many projections are selected for
  for(i in 1:dim(dir)[1]){

    vtx_projection <- complex$Vertices[,1:3]%*%dir[i,]
    if (ball == TRUE){
      buckets <- seq(-ball_radius,ball_radius,length.out = len+1)
    }
    else{
      buckets <- seq(min(vtx_projection),max(vtx_projection),length.out = len+1)
    }

    # map vertex projection to the feature index
    projection_bucket <- cut(vtx_projection, buckets, labels = FALSE)

    # update index to reflect rate values
    projection_bucket <- projection_bucket + (i - 1)*len

    selected_vertices[[i]] <- which(projection_bucket %in% indices)

  }
  final_selected_vertices <- Reduce(reduction_operation,selected_vertices)

  return(final_selected_vertices)
}


