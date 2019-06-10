
#### Generate ROC curves ####

#' Generates ROC curve averaged over multiple runs.
#'
#' @description Generates ROC curve averaged over multiple runs. We specify in the function what shapes to simulate,
#' paramters for the about the EC comptutation, as well as assessment scheme.
#' @export
#'
#' @param runs (int): Number of runs to average the curves over
#' @param num_sim (int) : The number of replicates of data.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param grid_size (int) : The fine-ness/granularity of the interpolated shapes.
#' @param distance_to_causal_point (float) : For interpolated shapes, the distance from a vertex to the causal points to be considered a "causal vertex"
#' @param causal_points (int) : The number of causal points in each causal cusp, or number of causal points used for interpolations.
#' @param shared_points (int) : The number of shared points in the shared cusps, or the number of shared points in the interpolations.
#' @param num_cones (int): The number of cones to compute the (S/D) EC curves for the generated shapes over.
#' @param eta (float) : The kernel shape parameter.
#' @param truncated (int) : The number of "cuts" to compute TPR/FPR for the ROC curve over. Used to speed up ROC computations.
#' @param two_curves (boolean) : Whether or not to compute ROC curves using class specific causal points, or the set of all causal points.
#' Setting two_curves = TRUE will provide two curves, for each class.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param type (string) : The assessment scheme. We currently support 'vertex' (finding causal vertices), 'feature' (finding causal sub-level sets),
#' 'cusp' (finding causal cusps for spheres).
#' @param min_points (int) : Used when type = 'feature'. The mininum number of causal vertices for a sub-level set to be associated with to be considered a causal 'feature'.
#' @param directions_per_cone (int): The number of directions we want generated within each cone.
#' @param cap_radius (float): The radius of the cones we generate (determines the size of each cone).
#' @param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include (during reconstruction).
#' @param mode (string) : The data generation scheme. We currently support 'sphere', 'gaussian_grid", or rbf interpolations (default).
#' @param subdivision (int) : The fineness of the sphere meshes (if mode == 'sphere'). We currently use subdivision = 3.
#' @param num_causal_region (int) : The number of causal cusps (for when mode == 'sphere').
#' @param num_shared_region (int) : The number of shared cusps (for when mode == 'sphere').
#' @param ec_type (string) : The type of EC we are computing. We currently support ECT, DECT and SECT.
generate_averaged_ROC_with_coned_directions  <- function(runs = 5, nsim = 50, curve_length = 10, grid_size = 25, distance_to_causal_point = 0.1, causal_points = 10,
                                                         shared_points = 3, num_cones = 5, eta = 0.1, truncated = FALSE, two_curves = FALSE,
                                                         ball_radius = 2, ball = TRUE, type = 'vertex',min_points = 2,directions_per_cone = 5, cap_radius = 0.15,
                                                         radius = 0, mode = 'sphere',
                                                         subdivision = 3, num_causal_region = 5, num_shared_region = 5,
                                                         ec_type = 'ECT'){
  if (type == 'vertex'){
    roc_curves = list()
    roc_curves2 = list()
    roc_curves3 = list()
    i = 1
    while (i<runs+1){
      roc_curve = try(generate_ROC_with_coned_directions(nsim = nsim, curve_length = curve_length, grid_size = grid_size, distance_to_causal_point = distance_to_causal_point,
                                                         causal_points = causal_points,shared_points = shared_points,num_cones = num_cones,
                                                         eta = eta,truncated = truncated, two_curves = two_curves,
                                                         ball = ball, ball_radius = ball_radius,type = type, min_points = min_points,num_causal_region = num_causal_region,
                                                         num_shared_region = num_shared_region,cap_radius = cap_radius,radius = radius,
                                                         directions_per_cone = directions_per_cone, mode = mode,ec_type = ec_type))
      if (inherits(roc_curve,'try-error')){
        next
      }
      else{
        roc_curves[[i]] = roc_curve
        if (two_curves == TRUE){
          roc_curves2[[i]] = as.matrix(roc_curve[which(roc_curve[,3]==1),])[,1:2]
          roc_curves3[[i]] = as.matrix(roc_curve[which(roc_curve[,3]==2),])[,1:2]
        }
        else{
          roc_curves2[[i]] = as.matrix(roc_curve[which(roc_curve[,3]==0),])[,1:2]
        }
        i = i+1
      }
    }
    total_roc = matrix(0, nrow = dim(roc_curves[[1]]), ncol = dim(roc_curves[[1]])[2])
    for (j in 1:runs){
      total_roc = total_roc + roc_curves[[j]]
    }
    total_roc = total_roc/runs
    return(total_roc)
  }
  if (type == 'feature' | type == 'cone' | type == 'cusp'){
    roc_curves = list()
    i = 1
    while (i<runs+1){
      roc_curve = try(generate_ROC_with_coned_directions(nsim = nsim, curve_length = curve_length, grid_size = grid_size, distance_to_causal_point = distance_to_causal_point,
                                                         causal_points = causal_points,shared_points = shared_points,num_cones = num_cones,
                                                         eta = eta,truncated = truncated, two_curves = two_curves,
                                                         ball = ball, ball_radius = ball_radius,type = type, min_points = min_points,num_causal_region = num_causal_region,
                                                         num_shared_region = num_shared_region,cap_radius = cap_radius,radius = radius,
                                                         directions_per_cone = directions_per_cone, mode = mode,ec_type = ec_type))
      if (inherits(roc_curve,'try-error')){
        next
      }
      else{
        roc_curves[[i]] = roc_curve
      }
      i = i+1
    }
    total_roc = matrix(0, nrow = dim(roc_curves[[1]]), ncol = dim(roc_curves[[1]])[2])
    for (j in 1:runs){
      total_roc = total_roc + roc_curves[[j]]
    }
    total_roc = total_roc/runs
    return(total_roc)
  }
}


#' This function generates an ROC curve using the cone reconstruction idea.
#'
#' @description The set of directions that we choose are grouped into cones, the centers of which are random. To select the vertices that are the output
#' of the reconstruction process, we consider evidence restricted only to the cones. For each cone of directions, take all the vertices whose
#' projections onto each of these directions in the cone are selected by the GPC/RATE procedure. Now take the union of such vertices, over all
#' the cones in our set of directions.
#'
#' For the actual ROC idea: we start with the RATE Values for each feature, and vary the threshold 1/p at which to consider the RATE values significant.
#' Consider a threshold t, and the set of all features whose rate values are above this threshold. For this collection of features, do the cone reconstruction
#' process outlined above to obtain a collection of vertices, which we consider to be 'positive'. Those that aren't selected are the 'negative' ones. We regard
#' a vertex to be a True Positive if it is within some small distance of a causal point, and conversely with a False Positive. True Negative and False Negative
#' vertices are defined similarly.
#' @export
#'
#' @param nsim (int) : The number of replicates of data.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param grid_size (int) : The fine-ness/granularity of the interpolated shapes.
#' @param distance_to_causal_point (float) : For interpolated shapes, the distance from a vertex to the causal points to be considered a "causal vertex"
#' @param causal_points (int) : The number of causal points in each causal cusp, or number of causal points used for interpolations.
#' @param shared_points (int) : The number of shared points in the shared cusps, or the number of shared points in the interpolations.
#' @param num_cones (int): The number of cones to compute the (S/D) EC curves for the generated shapes over.
#' @param eta (float) : The kernel shape parameter.
#' @param truncated (int) : The number of "cuts" to compute TPR/FPR for the ROC curve over. Used to speed up ROC computations.
#' @param two_curves (boolean) : Whether or not to compute ROC curves using class specific causal points, or the set of all causal points.
#' Setting two_curves = TRUE will provide two curves, for each class.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param type (string) : The assessment scheme. We currently support 'vertex' (finding causal vertices), 'feature' (finding causal sub-level sets),
#' 'cusp' (finding causal cusps for spheres).
#' @param min_points (int) : Used when type = 'feature'. The mininum number of causal vertices for a sub-level set to be associated with to be considered a causal 'feature'.
#' @param directions_per_cone (int): The number of directions we want generated within each cone.
#' @param cap_radius (float): The radius of the cones we generate (determines the size of each cone).
#' @param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include (during reconstruction).
#' @param ec_type (string) : The type of EC we are computing. We currently support ECT, DECT and SECT.
#' @param mode (string) : The data generation scheme. We currently support 'sphere', 'gaussian_grid", or interpolations (default).
#' @param subdivision (int) : The fineness of the sphere meshes (if mode == 'sphere'). We currently use subdivision = 3.
#' @param num_causal_region (int) : The number of causal cusps (for when mode == 'sphere').
#' @param num_shared_region (int) : The number of shared cusps (for when mode == 'sphere').
#'
#' @return roc_curve (matrix): The ROC curve for both classes of shapes

generate_ROC_with_coned_directions <- function(nsim = 10, curve_length = 25, grid_size = 25, distance_to_causal_point = 0.1,
                                               causal_points = 10,shared_points = 3, num_cones = 5, eta = 0.1,
                                               truncated = 300, two_curves = TRUE, ball = TRUE, ball_radius = 2.5, type = 'vertex',
                                               min_points = 3,directions_per_cone = 4, cap_radius = 0.15, radius = 0,ec_type = 'ECT',
                                               mode = 'sphere',
                                               subdivision = 3,num_causal_region = 5, num_shared_region = 5){
  print("generating directions")


  print("generating data")
  # generate data

  if (mode == 'sphere'){
    #cusps = 2*num_causal_region + num_shared_region + 1
    cusps = 2*num_causal_region + num_shared_region + 1
    causal_dirs = generate_equidistributed_points(cusps,cusps)
    causal_regions_1 = sample(1:cusps,num_causal_region)
    causal_regions_2 = sample((1:cusps)[-causal_regions_1],num_causal_region)
    shared_regions = sample((1:cusps)[-c(causal_regions_1,causal_regions_2)],num_shared_region)
    directions <- generate_equidistributed_cones(num_cones,cap_radius,directions_per_cone)
    data <- generate_data_sphere_simulation(nsim = nsim,dir = directions, curve_length = curve_length,noise_points = shared_points,
                                           causal_points = causal_points, ball_radius = ball_radius, subdivision = subdivision,
                                           cusps = cusps, causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                           shared_regions = shared_regions, ec_type = ec_type)
    directions <- directions
    ec_curve_data <- data$data

  }
  else if (mode == 'gaussian_grid'){
    print(mode)
    directions <- generate_equidistributed_cones(num_cones,cap_radius,directions_per_cone)
    data <- generate_data_gaussian_field(nsim = nsim, dir = directions, curve_length = curve_length,shared_points = shared_points,
                                     causal_points = causal_points,grid_size = grid_size,ball_radius = ball_radius,
                                     ec_type = ec_type)
    directions <- directions
    ec_curve_data <- data$data
  }
  else{
    print(mode)
    directions <- generate_equidistributed_cones(num_cones,cap_radius,directions_per_cone)
    data <- create_data_normal_fixed(num_sim = nsim, dir = directions, curve_length = curve_length,shared_points = shared_points,
                                     causal_points = causal_points,grid_size = grid_size,eta = eta,ball_radius = ball_radius,
                                     ec_type = ec_type)
    directions <- directions
    ec_curve_data <- data$data
  }

  num_cones <- dim(directions)[1]/directions_per_cone

  print("getting rate values")
  rate_values <- find_rate_variables_with_other_sampling_methods(ec_curve_data, bandwidth = 0.01, type = 'ESS')[,2]
  #Indices for Two Classes
  index1 = seq(1,nsim,2)
  complex_points1 = data$complex_points[index1]

  index2 = seq(2,nsim,2)
  complex_points2 = data$complex_points[index2]
  #Compute ROC using training data
  if (type == 'vertex'){
    if (two_curves == TRUE){
      roc_curve1 =  compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                             curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                             eta = eta, directions_per_cone = directions_per_cone, directions = directions, class = 1,truncated = truncated,
                                             ball_radius = ball_radius, radius = radius, mode = mode,subdivision = subdivision)
      roc_curve1 = cbind(roc_curve1, rep(1,dim(roc_curve1)[1]))
      roc_curve1 = cbind(roc_curve1,(1:dim(roc_curve1)[1]))

      roc_curve2 =  compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                             curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                             eta = eta, directions_per_cone = directions_per_cone, directions = directions,class = 2,truncated = truncated,
                                             ball_radius = ball_radius, radius = radius, mode = mode, subdivision = subdivision)
      roc_curve2 = cbind(roc_curve2, rep(2,dim(roc_curve2)[1]))
      roc_curve2 = cbind(roc_curve2,(1:dim(roc_curve2)[1]))

      roc_curve = rbind(roc_curve1,roc_curve2)
    }
    else{
      roc_curve = compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                           curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                           eta = eta, directions_per_cone = directions_per_cone, directions = directions,  class = 0, truncated = truncated,
                                           ball = ball, ball_radius = ball_radius, radius = radius, mode = mode, subdivision = subdivision)
      roc_curve = cbind(roc_curve,(1:dim(roc_curve)[1]))
    }
    return(roc_curve)
  }
  if (type == 'feature'){
    roc_curve = compute_roc_curve_features(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                           curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                           eta = eta, directions_per_cone = directions_per_cone, class = 0, truncated = truncated,
                                           ball = ball,ball_radius = ball_radius,
                                           dir = directions, min_points = min_points,mode = mode, subdivision = subdivision)
    roc_curve = cbind(roc_curve,(1:dim(roc_curve)[1]))
    return(roc_curve)
  }
  if (type == 'cone'){
    print('cone')
    roc_curve = compute_roc_curve_cones(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                        curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                        eta = eta, directions_per_cone = directions_per_cone, class = 0, truncated = truncated,
                                        ball = ball, ball_radius = ball_radius,
                                        dir = directions,  min_points = min_points, radius = radius,mode = mode, subdivision = subdivision)
    roc_curve = cbind(roc_curve,(1:dim(roc_curve)[1]))
    return(roc_curve)
  }
  if (type == 'cusp'){
    print('cusp')
    roc_curve = compute_roc_curve_modified_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                                  curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                                  eta = eta, directions_per_cone = directions_per_cone, class = 0, truncated = truncated,
                                                  ball = ball, ball_radius = ball_radius,
                                                  dir = directions, radius = radius,mode = mode, subdivision = subdivision)
    roc_curve = cbind(roc_curve,(1:dim(roc_curve)[1]))
    return(roc_curve)
  }
}

#'Computes the ROC curve by assessing the overlap of reconstructed vertices and causal vertices.
#'
#' @description We compute the ROC curve by assessing the overlap of reconstructed vertices and causal vertices.
#'We do this for every complex in the data set then average the ROC curves.
#'
#' @export
#'
#'
#' @param data (list) : Metadata about the simulated shapes (vertex coordinates, etc.)
#' @param class_1_causal_points : Vertex indices of causal points for class 1.
#' @param class_2_causal_points : Vertex indices of causal points for class 2.
#' @param distance_to_causal_point (float) : For interpolated shapes, the distance from a vertex to the causal points to be considered a "causal vertex"
#' @param rate_values (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#' @param grid_size (int) : The fine-ness/granularity of the interpolated shapes.
#' @param eta (float) : The kernel shape parameter.
#' @param directions_per_cone (int): The number of directions we want generated within each cone.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param directions (nx3 matrix):  The matrix of directions for which the (S/D) EC curve were computed over.
#' @param truncated (int) : The number of "cuts" to compute TPR/FPR for the ROC curve over. Used to speed up ROC computations.
#' @param class (int) : The class of the group of shapes we compute the ROC curve against.
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include (during reconstruction).
#' @param mode (string) : The data generation scheme. We currently support 'sphere', 'gaussian_grid", or interpolations (default).
#' @param subdivision (int) : The fineness of the sphere meshes (if mode == 'sphere'). We currently use subdivision = 3.
#'
#' @return total_rate_roc (matrix) : The ROC curve.

compute_roc_curve_vertex = function(data,class_1_causal_points,class_2_causal_points,distance_to_causal_point = 0.1,
                                    rate_values,grid_size,eta = 0.1,directions_per_cone, curve_length,directions, truncated = -1,class = 0,
                                    ball_radius = ball_radius, ball = TRUE, radius = 0, mode = 'sphere', subdivision = 3){
  print('Computing ROC curve...')
  #Initializing the number of vertices
  num_vertices = grid_size^2
  data_points = data$complex_points
  remove = c()
  counter = 0

  #Initializing the aggregate ROC curve frame
  if (truncated == -1){
    total_rate_roc = matrix(0, nrow = length(rate_values),ncol = 2)
  }
  else{
    total_rate_roc = matrix(0, nrow = min(truncated,length(rate_values)),ncol = 2)
  }

  roc_list = list()

  for (j in 1:length(data_points)){
    if (class == 1 && mod(j,2) != 1){
      next
    }
    if (class == 2 && mod(j,2) != 0){
      next
    }
    if (mode == 'grid'){
      class_1_causal_points = data$causal_points1
      class_2_causal_points = data$causal_points2
      predictions=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=data_points[[j]],eta=eta)
      complex=matrix_to_simplicial_complex(predictions,grid_length=grid_size)

      #Starting to Compute the ROC curve for a given complex
      class_1_true_vertices = c()
      class_2_true_vertices = c()

      for (j in 1:num_vertices){
        #computes the 2D euclidean distance on the grid between the points
        dist1=apply(X = class_1_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:2])
        dist2=apply(X = class_2_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:2])

        if (min(dist1)< distance_to_causal_point) class_1_true_vertices=c(class_1_true_vertices,j)
        if (min(dist2)< distance_to_causal_point) class_2_true_vertices=c(class_2_true_vertices,j)
      }
    }

    if (mode == 'gaussian_grid'){
      class_1_causal_points = data$causal_points1
      class_2_causal_points = data$causal_points2
      predictions=generate_random_field_matrix(grid_length=grid_size, points=data_points[[j]], 0.01)
      complex=matrix_to_simplicial_complex(predictions,grid_length=grid_size)

      #Starting to Compute the ROC curve for a given complex
      class_1_true_vertices = c()
      class_2_true_vertices = c()

      for (j in 1:num_vertices){
        #computes the 2D euclidean distance on the grid between the points
        dist1=apply(X = class_1_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:2])
        dist2=apply(X = class_2_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:2])

        if (min(dist1)< distance_to_causal_point) class_1_true_vertices=c(class_1_true_vertices,j)
        if (min(dist2)< distance_to_causal_point) class_2_true_vertices=c(class_2_true_vertices,j)
      }
    }

    if (mode == 'sphere'){
      class_1_true_vertices = data$causal_points1
      class_2_true_vertices = data$causal_points2

      shared_vertices = data$noise
      shared_points = data$shared_points_list

      sphere1 = vcgSphere(subdivision = subdivision)
      #sphere2 = vcgSphere(subdivision = subdivision)

      #sphere1$vb[1:3,class_1_true_vertices] = t(data_points[[j]])
      num_vertices = dim(sphere1$vb)[2]
      #sphere2$vb[1:3,class_2_true_vertices] = t(data_points[[j+1]])
      #sphere1$vb[1:3,shared_vertices] = shared_points[[(j+1)/2]]
      #sphere2$vb[1:3,shared_vertices] = shared_points[[(j+1)/2]]
      sphere1$vb[1:3,] = t(data_points[[j]])

      complex = convert_off_file(sphere1)

      #class_1_causal_points = data$causal_points1
      #class_2_causal_points = data$causal_points2
      #class_1_true_vertices = c()
      #class_2_true_vertices = c()
      #
      #for (j in 1:num_vertices){
      #  #computes the 2D euclidean distance on the grid between the points
      #  dist1=apply(X = t(sphere1$vb[1:3,data$causal_points1]),MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:3])
      #  dist2=apply(X = t(sphere1$vb[1:3,data$causal_points2]),MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:3])
      #
      #  if (min(dist1)< distance_to_causal_point) class_1_true_vertices=c(class_1_true_vertices,j)
      #  if (min(dist2)< distance_to_causal_point) class_2_true_vertices=c(class_2_true_vertices,j)
      #}
    }

    combined_true_vertices = union(class_1_true_vertices,class_2_true_vertices)

    class_1_false_vertices = setdiff(1:num_vertices, class_1_true_vertices)
    class_2_false_vertices = setdiff(1:num_vertices, class_2_true_vertices)

    combined_false_vertices = setdiff(1:num_vertices, combined_true_vertices)
    true_vertices = combined_true_vertices
    false_vertices = combined_false_vertices

    rate_ROC <- matrix(0,nrow = 0,ncol = 2)


    if (length(true_vertices) == 0 || length(false_vertices) == 0){
      remove = c(remove,j)
      next
    }

    counter = counter + 1

    # build the ROC by varying the ROC; we bucket the rate values into quantiles and select the thresholds that way; should make length.out = 1000, or higher
    # can also recover the case where we add rate values one at a time by taking length.out to be the number of rate values.
    if (truncated == -1){
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){

        #sink("/dev/null")
        rate_positive_vertices <- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                 len = curve_length, threshold = threshold,
                                                                 cone_size = directions_per_cone, ball_radius = ball_radius, radius = radius)
        print(length(rate_positive_vertices))
        #sink()

        rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

        # calculate the TPR, FPR
        # To do the case where we consider true vertices to be close to any causal point, replace class_1, class_2 true vertices with the combined_true vertices
        # Otherwise, use class_1_true_Vertices, class_2_true_vertices
        if (class == 0)
        {
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        combined_true_vertices,combined_false_vertices))
          true_vertices = combined_true_vertices
          false_vertices = combined_false_vertices
        }
        else if(class == 1){
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        class_1_true_vertices,class_1_false_vertices))

          true_vertices = class_1_true_vertices
          false_vertices = class_1_false_vertices
        }
        else{
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        class_2_true_vertices,class_2_false_vertices))
          true_vertices = class_2_true_vertices
          false_vertices = class_2_false_vertices
        }
        if (isTRUE((all.equal(c(1,1),previous_tpr_fpr)))){
          rate_ROC2 = matrix(1,ncol = 2, nrow = dim(total_rate_roc) - dim(rate_ROC)[1])
          rate_ROC = rbind(rate_ROC,rate_ROC2)
          break
        }
      }
    }
    else{
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated))){


        rate_positive_vertices<- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                 len = curve_length, threshold = threshold,
                                                                 cone_size = directions_per_cone, ball_radius = ball_radius, radius = radius)


        rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

        # calculate the TPR, FPR
        # To do the case where we consider true vertices to be close to any causal point, replace class_1, class_2 true vertices with the combined_true vertices
        # Otherwise, use class_1_true_Vertices, class_2_true_vertices
        if (class == 0)
        {
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        combined_true_vertices,combined_false_vertices))
          true_vertices = combined_true_vertices
          false_vertices = combined_false_vertices
        }
        else if(class == 1){
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        class_1_true_vertices,class_1_false_vertices))
          true_vertices = class_1_true_vertices
          false_vertices = class_1_false_vertices
        }
        else{
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        class_2_true_vertices,class_2_false_vertices))
          true_vertices = class_2_true_vertices
          false_vertices = class_2_false_vertices
        }
        previous_tpr_fpr = calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                             true_vertices,false_vertices)
        if (isTRUE((all.equal(c(1,1),previous_tpr_fpr)))){
          rate_ROC2 = matrix(1,ncol = 2, nrow = dim(total_rate_roc) - dim(rate_ROC)[1])
          rate_ROC = rbind(rate_ROC,rate_ROC2)
          break
        }
      }
    }

    total_rate_roc = total_rate_roc + rate_ROC


  }
  total_rate_roc = (total_rate_roc / counter)



  print(total_rate_roc)

  return(total_rate_roc)
}

#'Computes the ROC curve by assessing the overlap of reconstructed cusps and causal cusps.
#'
#' @export
#' @description We compute the ROC curve by assessing the overlap of reconstructed cusps and causal cusps. Reconstructing a causal cusp here means
#' identifying one vertex that is in the causal cusp.
#'We do this for every complex in the data set then average the ROC curves.
#'
#' @param data (list) : Metadata about the simulated shapes (vertex coordinates, etc.)
#' @param class_1_causal_points : Vertex indices of causal points for class 1.
#' @param class_2_causal_points : Vertex indices of causal points for class 2.
#' @param distance_to_causal_point (float) : For interpolated shapes, the distance from a vertex to the causal points to be considered a "causal vertex"
#' @param rate_values (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#' @param grid_size (int) : The fine-ness/granularity of the interpolated shapes.
#' @param eta (float) : The kernel shape parameter.
#' @param directions_per_cone (int): The number of directions we want generated within each cone.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param directions (nx3 matrix):  The matrix of directions for which the (S/D) EC curve were computed over.
#' @param truncated (int) : The number of "cuts" to compute TPR/FPR for the ROC curve over. Used to speed up ROC computations.
#' @param class (int) : The class of the group of shapes we compute the ROC curve against.
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include (during reconstruction).
#' @param mode (string) : The data generation scheme. We currently support 'sphere', 'gaussian_grid", or interpolations (default).
#' @param subdivision (int) : The fineness of the sphere meshes (if mode == 'sphere'). We currently use subdivision = 3.
#'
#' @return total_rate_roc (matrix) : The ROC curve.
compute_roc_curve_modified_vertex = function(data,class_1_causal_points,class_2_causal_points,distance_to_causal_point = 0.1,
                                             rate_values,grid_size,eta = 0.1,directions_per_cone, curve_length,directions, truncated = 0,class = 0,
                                             ball_radius = ball_radius, ball = TRUE, radius = 0, mode = 'grid', subdivision = 3){
  print('Computing ROC curve...')
  #Initializing the number of vertices
  num_vertices = grid_size^2
  data_points = data$complex_points
  remove = c()

  #Initializing the aggregate ROC curve frame
  if (truncated == 0){
    total_rate_roc = matrix(0, nrow = length(rate_values),ncol = 2)
  }
  else{
    total_rate_roc = matrix(0, nrow = min(truncated,length(rate_values)),ncol = 2)
  }
  roc_list = list()
  for (j in 1:length(data_points)){
    if (mode == 'grid'){
      class_1_causal_points = data$causal_points1
      class_2_causal_points = data$causal_points2
      predictions=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=data_points[[i]],eta=eta)
      complex=matrix_to_simplicial_complex(predictions,grid_length=grid_size)

      #Starting to Compute the ROC curve for a given complex
      class_1_true_vertices = c()
      class_2_true_vertices = c()

      for (j in 1:num_vertices){
        #computes the 2D euclidean distance on the grid between the points
        dist1=apply(X = class_1_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:2])
        dist2=apply(X = class_2_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:2])

        if (min(dist1)< distance_to_causal_point) class_1_true_vertices=c(class_1_true_vertices,j)
        if (min(dist2)< distance_to_causal_point) class_2_true_vertices=c(class_2_true_vertices,j)
      }
    }

    if (mode == 'sphere'){
      class_1_true_vertices = data$causal_points1
      class_2_true_vertices = data$causal_points2
      total_cusps = data$total_cusps_list

      shared_vertices = data$noise
      shared_points = data$shared_points_list

      sphere1 = vcgSphere(subdivision = subdivision)
      sphere2 = vcgSphere(subdivision = subdivision)

      #sphere1$vb[1:3,class_1_true_vertices] = t(data_points[[j]])
      num_vertices = dim(sphere1$vb)[2]
      #sphere2$vb[1:3,class_2_true_vertices] = t(data_points[[j+1]])
      #sphere1$vb[1:3,shared_vertices] = shared_points[[(j+1)/2]]
      #sphere2$vb[1:3,shared_vertices] = shared_points[[(j+1)/2]]
      sphere1$vb[1:3,] = t(data_points[[i]])

      complex = convert_off_file(sphere1)

      #class_1_causal_points = data$causal_points1
      #class_2_causal_points = data$causal_points2
      #class_1_true_vertices = c()
      #class_2_true_vertices = c()
      #
      #for (j in 1:num_vertices){
      #  #computes the 2D euclidean distance on the grid between the points
      #  dist1=apply(X = t(sphere1$vb[1:3,data$causal_points1]),MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:3])
      #  dist2=apply(X = t(sphere1$vb[1:3,data$causal_points2]),MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:3])
      #
      #  if (min(dist1)< distance_to_causal_point) class_1_true_vertices=c(class_1_true_vertices,j)
      #  if (min(dist2)< distance_to_causal_point) class_2_true_vertices=c(class_2_true_vertices,j)
      #}
    }
    combined_true_vertices = union(class_1_true_vertices,class_2_true_vertices)
    class_1_false_vertices = setdiff(1:num_vertices, class_1_true_vertices)
    class_2_false_vertices = setdiff(1:num_vertices, class_2_true_vertices)
    combined_false_vertices = setdiff(1:num_vertices, combined_true_vertices)
    true_vertices = combined_true_vertices
    false_vertices = combined_false_vertices
    #if(mod(i,2) == 1){
    #  true_vertices = class_1_true_vertices
    #  false_vertices = class_1_false_vertices
    #}
    #else{
    #  true_vertices = class_2_true_vertices
    #  false_vertices = class_2_false_vertices
    #}
    if (length(true_vertices) == 0 || length(false_vertices) == 0){
      remove = c(remove,i)
      next
    }
    # build the ROC by varying the ROC; we bucket the rate values into quantiles and select the thresholds that way; should make length.out = 1000, or higher
    # can also recover the case where we add rate values one at a time by taking length.out to be the number of rate values.
    if (truncated == 0){
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){


        rate_positive_vertices<- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                 len = curve_length, threshold = threshold,
                                                                 cone_size = directions_per_cone, ball_radius = ball_radius,ball = ball, radius = radius)

        rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

        # calculate the TPR, FPR
        # To do the case where we consider true vertices to be close to any causal point, replace class_1, class_2 true vertices with the combined_true vertices
        # Otherwise, use class_1_true_Vertices, class_2_true_vertices
        rate_positive_vertices <- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                  len = curve_length, threshold = threshold,
                                                                  cone_size = directions_per_cone, ball = ball, ball_radius = ball_radius, radius = radius)

        rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

        TPR_FPR <- calculate_TPR_FPR_cusps(rate_positive_vertices,rate_negative_vertices,
                                           total_cusps,false_vertices)
        rate_ROC <- rbind(rate_ROC, TPR_FPR)
      }
    }
    else{
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated))){


        rate_positive_vertices<- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                 len = curve_length, threshold = threshold,
                                                                 cone_size = directions_per_cone, ball_radius = ball_radius,ball = ball, radius = radius)

        rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

        rate_positive_vertices <- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                  len = curve_length, threshold = threshold,
                                                                  cone_size = directions_per_cone, ball = ball, ball_radius = ball_radius, radius = radius)

        rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

        TPR_FPR <- calculate_TPR_FPR_cusps(rate_positive_vertices,rate_negative_vertices,
                                           total_cusps,false_vertices)
        rate_ROC <- rbind(rate_ROC, TPR_FPR)
      }
    }

    total_rate_roc = total_rate_roc + rate_ROC

  }
  total_rate_roc = (total_rate_roc / length(data_points))



  print(total_rate_roc)

  return(total_rate_roc)
}

#'Computes the ROC curve by assessing the overlap of selected features and causal features.
#'
#' @export
#' @description We compute the ROC curve by assessing the overlap of selected features and causal features. A causal feature here
#' means to be associated to more than the (min_points) number of causal vertices.
#'We do this for every complex in the data set then average the ROC curves.
#'
#' @param data (list) : Metadata about the simulated shapes (vertex coordinates, etc.)
#' @param class_1_causal_points : Vertex indices of causal points for class 1.
#' @param class_2_causal_points : Vertex indices of causal points for class 2.
#' @param distance_to_causal_point (float) : For interpolated shapes, the distance from a vertex to the causal points to be considered a "causal vertex"
#' @param rate_values (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#' @param grid_size (int) : The fine-ness/granularity of the interpolated shapes.
#' @param eta (float) : The kernel shape parameter.
#' @param directions_per_cone (int): The number of directions we want generated within each cone.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param truncated (int) : The number of "cuts" to compute TPR/FPR for the ROC curve over. Used to speed up ROC computations.
#' @param class (int) : The class of the group of shapes we compute the ROC curve against.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param dir (nx3 matrix):  The matrix of directions for which the (S/D) EC curve were computed over.
#' @param min_points (int) : Used when type = 'feature'. The mininum number of causal vertices for a sub-level set to be associated with to be considered a causal 'feature'.
#' @param mode (string) : The data generation scheme. We currently support 'sphere', 'gaussian_grid", or interpolations (default).
#' @param subdivision (int) : The fineness of the sphere meshes (if mode == 'sphere'). We currently use subdivision = 3.
#'
#' @return total_rate_roc (matrix) : The ROC curve.
compute_roc_curve_features <- function(data,class_1_causal_points,class_2_causal_points,distance_to_causal_point = 0.1,
                                       rate_values,grid_size,eta = 0.1,directions_per_cone, curve_length, truncated = -1,
                                       class = 0, ball = TRUE, ball_radius = ball_radius, dir, min_points = 2,
                                       mode = 'grid', subdivision = 3){
  #Assign index to Rate Values
  print('Computing ROC curve...')
  data_points = data$complex_points

  #Initializing the number of vertices
  num_vertices = grid_size^2
  #Initializing the aggregate ROC curve frame
  if (truncated == 0){
    total_rate_roc = matrix(0, nrow = length(rate_values)+1,ncol = 2)
  }
  else{
    total_rate_roc = matrix(0, nrow = truncated+1,ncol = 2)
  }
  roc_list = list()
  num_tests = 0
  for (j in seq(1,length(data_points),2)){

    #Interpolating based on the causal and shared points in R^3 for each shape
    if (mode == 'grid'){
      num_vertices = grid_size^2
      class_1_causal_points = data$causal_points1
      class_2_causal_points = data$causal_points2
      predictions1=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=data_points[[j]],eta=eta)
      predictions2=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=data_points[[j+1]],eta=eta)

      #Produce the Simplicial Complex
      complex1=MatrixtoSimplicialComplexTriangular(predictions1,grid_length=grid_size)
      complex2=MatrixtoSimplicialComplexTriangular(predictions2,grid_length=grid_size)

      #Starting to Compute the ROC curve for a given complex
      class_1_true_vertices = c()
      class_2_true_vertices = c()

      for (i in 1:num_vertices){
        #computes the 2D euclidean distance on the grid between the points
        dist1=apply(X = class_1_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex1$Vertices[i,1:2])
        dist2=apply(X = class_2_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex2$Vertices[i,1:2])

        if (min(dist1)< distance_to_causal_point) class_1_true_vertices=c(class_1_true_vertices,i)
        if (min(dist2)< distance_to_causal_point) class_2_true_vertices=c(class_2_true_vertices,i)
      }
    }
    if (mode == 'sphere'){
      class_1_true_vertices = data$causal_points1
      class_2_true_vertices = data$causal_points2

      shared_vertices = data$noise
      shared_points = data$shared_points_list

      sphere1 = vcgSphere(subdivision = subdivision)
      sphere2 = vcgSphere(subdivision = subdivision)

      num_vertices = dim(sphere1$vb)[2]

      sphere1$vb[1:3,] = t(data_points[[j]])
      sphere2$vb[1:3,] = t(data_points[[j+1]])

      complex1 = convert_off_file(sphere1)
      complex2 = convert_off_file(sphere2)

      class_1_causal_points = data$causal_points1
      class_2_causal_points = data$causal_points2
    }
    combined_true_vertices = union(class_1_true_vertices,class_2_true_vertices)

    #class_1_false_vertices = setdiff(1:num_vertices, class_1_true_vertices)
    #class_2_false_vertices = setdiff(1:num_vertices, class_2_true_vertices)
    combined_false_vertices = setdiff(1:num_vertices, combined_true_vertices)

    rate_ROC <- matrix(0,nrow = 1,ncol = 2)
    feature_with_index  = cbind(1:length(rate_values),rate_values)
    #Order the Rate Values
    ordered_rate_values = feature_with_index[order(feature_with_index[,2], decreasing = TRUE),]
    vertex_list1 = feature_vertex_association(dir = dir, complex = complex1, len = curve_length,ball_radius = ball_radius, ball = ball)
    vertex_list2 = feature_vertex_association(dir = dir, complex = complex2, len = curve_length,ball_radius = ball_radius, ball = ball)


    #Check if features are positive or not.
    positive_features = c()
    for (i in 1:dim(ordered_rate_values)[1]){
      index=ordered_rate_values[i,1]
      found_vertices1 = intersect(vertex_list1[[index]],class_1_true_vertices)
      found_vertices2 = intersect(vertex_list2[[index]],class_2_true_vertices)
      if ((length(found_vertices1)+length(found_vertices2))>=min_points){
        positive_features = c(positive_features, index)
      }
    }
    empty_list = c()
    negative_features <- setdiff(1:length(rate_values),positive_features)
    if (length(positive_features) == 0 || length(negative_features) == 0){
      next
    }
    if (truncated == 0){
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values)))){
        selected_features = matrix(ordered_rate_values[(which(ordered_rate_values[,2]>=threshold)),],ncol=2)[,1]

        #rate_positive_features = c(empty_list,intersect(positive_features, selected_features))
        #rate_negative_features = c(empty_list,intersect(negative_features, selected_features))
        rate_positive_features = selected_features
        rate_negative_features = setdiff(1:length(rate_values), selected_features)
        rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_features,rate_negative_features,
                                                      positive_features,negative_features))
      }
    }
    else{
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated))){

        selected_features = matrix(ordered_rate_values[(which(ordered_rate_values[,2]>=threshold)),],ncol=2)[,1]

        #rate_positive_features = c()
        #rate_negative_features = c()
        #rate_positive_features = c(rate_positive_features,intersect(positive_features, selected_features))
        rate_positive_features = selected_features
        rate_negative_features = setdiff(1:length(rate_values), selected_features)
        #rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_features,rate_negative_features,
        #                                              positive_features,negative_features))
        rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_features,rate_negative_features,
                                                      positive_features,negative_features))
      }
    }
    # print(rate_ROC)
    # roc_list[[j]] = rate_ROC
    if (NaN %in% rate_ROC){
      next
    }
    else{
      num_tests = num_tests + 1
      total_rate_roc = total_rate_roc + rate_ROC
    }
  }
  total_rate_roc = (total_rate_roc / num_tests)
  print(total_rate_roc)
  if (NaN %in% total_rate_roc){
    stop('NAs in ROC, moving onto next iteration')
  }

  return(total_rate_roc)
}

#'Computes the ROC curve by assessing the overlap of selected bands and causal cone.
#'
#' @export
#'
#' @description We compute the ROC curve by assessing the overlap of selected features and causal features. A causal cone here
#' means to have every direction in a given cone that are associated with a causal vertex to be selected.
#' We do this for every complex in the data set then average the ROC curves.
#'
#' @param data (list) : Metadata about the simulated shapes (vertex coordinates, etc.)
#' @param class_1_causal_points : Vertex indices of causal points for class 1.
#' @param class_2_causal_points : Vertex indices of causal points for class 2.
#' @param distance_to_causal_point (float) : For interpolated shapes, the distance from a vertex to the causal points to be considered a "causal vertex"
#' @param rate_values (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#' @param grid_size (int) : The fine-ness/granularity of the interpolated shapes.
#' @param eta (float) : The kernel shape parameter.
#' @param directions_per_cone (int): The number of directions we want generated within each cone.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param truncated (int) : The number of "cuts" to compute TPR/FPR for the ROC curve over. Used to speed up ROC computations.
#' @param class (int) : The class of the group of shapes we compute the ROC curve against.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param dir (nx3 matrix):  The matrix of directions for which the (S/D) EC curve were computed over.
#' @param min_points (int) : Used when type = 'feature'. The mininum number of causal vertices for a sub-level set to be associated with to be considered a causal 'feature'.
#' @param mode (string) : The data generation scheme. We currently support 'sphere', 'gaussian_grid", or interpolations (default).
#' @param subdivision (int) : The fineness of the sphere meshes (if mode == 'sphere'). We currently use subdivision = 3.
#'
#' @return ROC_dataframe (dataframe) : The ROC curve in Dataframe form.

compute_roc_curve_cones <- function(data,class_1_causal_points,class_2_causal_points,distance_to_causal_point = 0.1,
                                    rate_values,grid_size,eta = 0.1,directions_per_cone, curve_length, truncated = 0,
                                    class = 0, ball = ball, ball_radius = ball_radius, dir, min_points = 2, radius = 2,
                                    mode = 'grid', subdivision = 3){
  #Assign index to Rate Values
  print('Computing ROC curve...')
  #Initializing the number of vertices
  num_vertices = grid_size^2
  data_points = data$complex_points
  #Initializing the aggregate ROC curve frame
  if (truncated == 0){
    total_rate_roc = matrix(0, nrow = length(rate_values)+1,ncol = 2)
  }
  else{
    total_rate_roc = matrix(0, nrow = truncated+1,ncol = 2)
  }
  num_tests = 0
  remove = c()
  for (j in seq(1,length(data_points),2)){
    if (mode == 'grid'){
      num_vertices = grid_size^2
      class_1_causal_points = data$causal_points1
      class_2_causal_points = data$causal_points2
      predictions1=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=data_points[[j]],eta=eta)
      predictions2=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=data_points[[j+1]],eta=eta)

      #Produce the Simplicial Complex
      complex1=MatrixtoSimplicialComplexTriangular(predictions1,grid_length=grid_size)
      complex2=MatrixtoSimplicialComplexTriangular(predictions2,grid_length=grid_size)

      #Starting to Compute the ROC curve for a given complex
      class_1_true_vertices = c()
      class_2_true_vertices = c()

      for (i in 1:num_vertices){
        #computes the 2D euclidean distance on the grid between the points
        dist1=apply(X = class_1_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex1$Vertices[i,1:2])
        dist2=apply(X = class_2_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex2$Vertices[i,1:2])

        if (min(dist1)< distance_to_causal_point) class_1_true_vertices=c(class_1_true_vertices,i)
        if (min(dist2)< distance_to_causal_point) class_2_true_vertices=c(class_2_true_vertices,i)
      }
    }
    if (mode == 'sphere'){
      class_1_true_vertices = data$causal_points1
      class_2_true_vertices = data$causal_points2

      shared_vertices = data$noise
      shared_points = data$shared_points_list

      sphere1 = vcgSphere(subdivision = subdivision)
      sphere2 = vcgSphere(subdivision = subdivision)

      #sphere1$vb[1:3,class_1_true_vertices] = t(data_points[[j]])
      num_vertices = dim(sphere1$vb)[2]
      #sphere2$vb[1:3,class_2_true_vertices] = t(data_points[[j+1]])
      #sphere1$vb[1:3,shared_vertices] = shared_points[[(j+1)/2]]
      #sphere2$vb[1:3,shared_vertices] = shared_points[[(j+1)/2]]
      sphere1$vb[1:3,] = t(data_points[[j]])
      sphere2$vb[1:3,] = t(data_points[[j+1]])

      complex1 = convert_off_file(sphere1)
      complex2 = convert_off_file(sphere2)

      class_1_causal_points = data$causal_points1
      class_2_causal_points = data$causal_points2
      # class_1_true_vertices = c()
      # class_2_true_vertices = c()
      #
      # for (i in 1:num_vertices){
      #   #computes the 2D euclidean distance on the grid between the points
      #   dist1=apply(X = t(sphere1$vb[1:3,data$causal_points1]),MARGIN = 1,FUN = difference,y=complex1$Vertices[i,1:3])
      #   dist2=apply(X = t(sphere1$vb[1:3,data$causal_points2]),MARGIN = 1,FUN = difference,y=complex2$Vertices[i,1:3])
      #
      #   if (min(dist1)< distance_to_causal_point) {
      #     class_1_true_vertices=c(class_1_true_vertices,i)
      #   }
      #   if (min(dist2)< distance_to_causal_point){
      #     class_2_true_vertices=c(class_2_true_vertices,i)
      #   }
      # }
    }

    combined_true_vertices = union(class_1_true_vertices,class_2_true_vertices)
    combined_false_vertices = setdiff(1:num_vertices, combined_true_vertices)

    rate_ROC <- matrix(0,nrow = 1,ncol = 2)

    feature_with_index  = cbind(1:length(rate_values),rate_values)
    ordered_rate_values = feature_with_index[order(feature_with_index[,2], decreasing = TRUE),]
    vertex_list1 = feature_vertex_association(dir = dir, complex = complex1, len = curve_length, ball_radius = ball_radius, ball = ball)
    vertex_list2 = feature_vertex_association(dir = dir, complex = complex2, len = curve_length, ball_radius = ball_radius, ball = ball)

    total_reconstructions = curve_length * (dim(dir)[1]/directions_per_cone)
    positive_reconstructions = c()
    negative_reconstructions = c()
    reconstruction_vertices1 = list()
    reconstruction_vertices2 = list()
    reconstruction_features  = list()
    num_cones = dim(dir)[1]/directions_per_cone
    for (i in 1:num_cones){
      if (num_cones == 1){
        cone_dirs = dir[1:directions_per_cone,]
      }
      else{
        cone_dirs=dir[((i-1)*(directions_per_cone)+1):(i*directions_per_cone),]
      }
      #Offset for the cones
      constant = (i-1) * curve_length * directions_per_cone
      for (l in 1:curve_length){
        reconstruction_vector = rep(0,curve_length*directions_per_cone)
        reconstruction_number = (i-1)*curve_length + l
        reconstruction_index  = constant + l
        reconstruction_features[[reconstruction_number]] = list()
        reconstruction_vector[reconstruction_index] = 1
        for (k in 1:(directions_per_cone)){
          radius_index = (reconstruction_index + (k-1)*curve_length  - radius):(reconstruction_index + (k-1)*curve_length + radius)
          #Finding the indices of a certain direction along the curve
          direction_min_index = (i-1) * curve_length * directions_per_cone + (k-1) * curve_length + 1
          direction_max_index = (i-1) * curve_length * directions_per_cone + (k) * curve_length
          direction_vector = direction_min_index:direction_max_index
          #Intersecting the Radius of points with the index
          radius_points = intersect(radius_index,direction_vector)
          radius_points = radius_points - constant
          reconstruction_vector[radius_points] = 1
          reconstruction_features[[reconstruction_number]][[k]] = radius_points
        }
        reconstruction_vertices1[[reconstruction_number]] = summarize_vertices(dir = cone_dirs,complex = complex1,rate_vals = reconstruction_vector,len = curve_length,
                                                                               threshold = 0.5,cone_size = directions_per_cone,ball_radius = ball_radius,ball = ball)
        reconstruction_vertices2[[reconstruction_number]] = summarize_vertices(dir = cone_dirs,complex = complex2,rate_vals = reconstruction_vector,len = curve_length,
                                                                               threshold = 0.5,cone_size = directions_per_cone,ball_radius = ball_radius,ball = ball)
      }
    }
    reconstruction_vertices = c()
    for (i in 1:total_reconstructions){
      found_vertices1 = intersect(reconstruction_vertices1[[i]],combined_true_vertices)
      found_vertices2 = intersect(reconstruction_vertices2[[i]],combined_true_vertices)
      if ((length(found_vertices1)+length(found_vertices2))>=min_points){
        positive_reconstructions = c(positive_reconstructions, i)
        reconstruction_vertices = c(reconstruction_vertices,found_vertices1,found_vertices2)
      }
      else{
        negative_reconstructions = c(negative_reconstructions, i)
      }
    }
    empty_list = c()
    if (length(positive_reconstructions) == 0 || length(negative_reconstructions) == 0){
      remove = c(remove,((j+1)/2))
      next
    }
    if (truncated == 0){
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values)))){
        selected_features = matrix(ordered_rate_values[(which(ordered_rate_values[,2]>=threshold)),],ncol=2)[,1]

        rate_positive_reconstructions = c()

        for (l in 1:total_reconstructions){
          if (check_reconstruction(selected_features =selected_features, reconstruction_features_index = reconstruction_features[[l]])){
            rate_positive_reconstructions = c(rate_positive_reconstructions,l)
          }
        }
        rate_negative_reconstructions = setdiff(1:total_reconstructions, rate_positive_reconstructions)
        rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_reconstructions,rate_negative_reconstructions,
                                                      positive_reconstructions,negative_reconstructions))
      }
      if (NaN %in% rate_ROC){
        next
      }
      else{
        num_tests = num_tests + 1
        total_rate_roc = total_rate_roc + rate_ROC
      }
    }
    else{
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated))){

        selected_features = matrix(ordered_rate_values[(which(ordered_rate_values[,2]>=threshold)),],ncol=2)[,1]

        rate_positive_reconstructions = c()

        for (l in 1:total_reconstructions){
          if (check_reconstruction(selected_features =selected_features, reconstruction_features_index = reconstruction_features[[l]])){
            rate_positive_reconstructions = c(rate_positive_reconstructions,l)
          }
        }
        rate_negative_reconstructions = setdiff(1:total_reconstructions, rate_positive_reconstructions)
        rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_reconstructions,rate_negative_reconstructions,
                                                      positive_reconstructions,negative_reconstructions))
      }
    }
    if (NaN %in% rate_ROC){
      next
    }
    else{
      num_tests = num_tests + 1
      total_rate_roc = total_rate_roc + rate_ROC
    }
  }
  total_rate_roc = (total_rate_roc / num_tests)
  if (num_tests == 0){
    stop('Bad Simulation!')
  }
  if (NaN %in% total_rate_roc){
    stop('NAs in ROC, moving onto next iteration')
  }

  print(total_rate_roc)
  # plot the ROC curve, should be a list of points
  # R staircase function might come in handy

  # append labels to each class' ROC for coloring purposes
  rate_ROC_labeled <- cbind(total_rate_roc, rep(class,dim(total_rate_roc)[1] ) )
  ROC_dataframe <- data.frame(rate_ROC_labeled)

  return(ROC_dataframe)
}

#'Computes the averaged ROC curve of caricatured teeth.
#'
#' @export
#' @description We compute the ROC curve by assessing the overlap of reconstructed vertices and causal vertices.
#'We do this for every tooth in the directory then average the ROC curves. The user must pass in the locations of the directories for the two meshes,
#'and also the vector of transition probabilities for the caricaturization procedure.
#'
#' @param data_dir1 (string) : Location of the first class of meshes.
#' @param data_dir2 (string) : Location of the second class of meshes.
#' @param gamma (float) : The value at which if a vertex has transition probability greater than $\gamma$, then the vertex is considered causal.
#' @param class_1_probs (vector/matrix) : The vector/matrix of transition probibilities to propogate the weights for class 1 shapes.
#' @param class_2_probs (vector/matrix) : The vector/matrix of transition probibilities to propogate the weights for class 2 shapes.
#' @param rate_values (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#' @param directions_per_cone (int): The number of directions we want generated within each cone.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param directions (nx3 matrix):  The matrix of directions for which the (S/D) EC curve were computed over.
#' @param truncated (int) : The number of "cuts" to compute TPR/FPR for the ROC curve over. Used to speed up ROC computations.
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include (during reconstruction).
#' @param two_curves (boolean) : Whether or not to compute ROC curves using class specific causal points, or the set of all causal points.
#'
#' @return roc_curve (matrix) : The roc curves for both classes of shapes.
compute_roc_curve_teeth = function(data_dir1, data_dir2, gamma, class_1_probs, class_2_probs,
                                   rate_values, directions_per_cone, curve_length,directions, truncated = 0,
                                   ball_radius = ball_radius, ball = TRUE , radius = 0,two_curves = FALSE){
  if (two_curves == FALSE){
    roc_curve_label1 = 1
    roc_curve_label2 = 2
  }
  else{
    roc_curve_label1 = 0
    roc_curve_label2 = 0
  }
  roc_curve1 =  compute_roc_curve_teeth_vertex(data_dir = data_dir1, gamma = gamma, class_1_probs = class_1_probs, class_2_probs = class_2_probs,
                                               curve_length = curve_length,  rate_values = rate_values,
                                               directions_per_cone = directions_per_cone, directions = directions, class = roc_curve_label1,truncated = truncated,
                                               ball_radius = ball_radius,radius = radius)
  roc_curve1 = cbind(roc_curve1, rep(1,dim(roc_curve1)[1]))

  roc_curve2 =  compute_roc_curve_teeth_vertex(data_dir = data_dir2, gamma = gamma, class_1_probs = class_1_probs, class_2_probs = class_2_probs,
                                               curve_length = curve_length,  rate_values = rate_values,
                                               directions_per_cone = directions_per_cone, directions = directions, class = roc_curve_label2,truncated = truncated,
                                               ball_radius = ball_radius, radius = radius)

  roc_curve2 = cbind(roc_curve2, rep(2,dim(roc_curve2)[1]))

  roc_curve = rbind(roc_curve1,roc_curve2)
  if (two_curves == FALSE){
    roc_curve = (roc_curve1 + roc_curve2)/2
    roc_curve[,3] = 0
  }
  return(roc_curve)}

#'Computes the ROC curve by assessing the overlap of reconstructed vertices and causal vertices on the teeth.
#' @export
#' @import Rvcg
#'
#' @description We compute the ROC curve by assessing the overlap of reconstructed vertices and causal vertices.
#'We do this for every tooth in the directory then average the ROC curves. The user must pass in the locations of the directories for the two meshes,
#'and also the vector of transition probabilities for the caricaturization procedure.
#'
#' @param data_dir (string) : Location of the meshes.
#' @param gamma (float) : The value at which if a vertex has transition probability greater than $\gamma$, then the vertex is considered causal.
#' @param class_1_probs (vector/matrix) : The vector/matrix of transition probibilities to propogate the weights for class 1 shapes.
#' @param class_2_probs (vector/matrix) : The vector/matrix of transition probibilities to propogate the weights for class 2 shapes.
#' @param rate_values (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#' @param directions_per_cone (int): The number of directions we want generated within each cone.
#' @param curve_length (int) : Number of sub-level sets in each EC computation.
#' @param directions (nx3 matrix):  The matrix of directions for which the (S/D) EC curve were computed over.
#' @param truncated (int) : The number of "cuts" to compute TPR/FPR for the ROC curve over. Used to speed up ROC computations.
#' @param class (int) : The class of the roc curve. If class = 0, we let the set of causal vertices be the union of causal vertices from class 1 and 2.
#' @param ball_radius (float) : The radius of the bounding ball used if we compute the balled EC curve.
#' @param ball (boolean) : Denotes whether or not to compute the EC curves over a ball for uniform measurements
#' @param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include (during reconstruction).
#'
#' @return total_rate_roc (matrix): The ROC curve.

compute_roc_curve_teeth_vertex = function(data_dir,gamma,class_1_probs,class_2_probs,
                                          rate_values,directions_per_cone, curve_length,directions, truncated = 0,class = 0,
                                          ball_radius = ball_radius, ball = TRUE , radius = 0){
  print('Computing ROC curve...')
  #Initializing the number of vertices
  remove = c()
  counter = 0

  #Initializing the aggregate ROC curve frame
  if (truncated == 0){
    total_rate_roc = matrix(0, nrow = length(rate_values),ncol = 2)
  }
  else{
    total_rate_roc = matrix(0, nrow = min(truncated,length(rate_values)),ncol = 2)
  }
  roc_list = list()
  files = list.files(data_dir,full.names = TRUE)
  for (j in 1:length(files)){
    print(j)
    mesh = vcgImport(files[j])
    complex = convert_off_file(mesh)
    num_vertices = dim(complex$Vertices)[1]
    class_1_true_vertices = which(class_1_probs > gamma)
    class_2_true_vertices = which(class_2_probs > gamma)

    combined_true_vertices = union(class_1_true_vertices,class_2_true_vertices)

    class_1_false_vertices = setdiff(1:num_vertices, class_1_true_vertices)
    class_2_false_vertices = setdiff(1:num_vertices, class_2_true_vertices)

    combined_false_vertices = setdiff(1:num_vertices, combined_true_vertices)
    true_vertices = combined_true_vertices
    false_vertices = combined_false_vertices

    rate_ROC <- matrix(0,nrow = 0,ncol = 2)


    if (length(true_vertices) == 0 || length(false_vertices) == 0){
      remove = c(remove,j)
      next
    }
    counter = counter + 1
    # build the ROC by varying the ROC; we bucket the rate values into quantiles and select the thresholds that way; should make length.out = 1000, or higher
    # can also recover the case where we add rate values one at a time by taking length.out to be the number of rate values.
    if (truncated == 0){
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){

        #sink("/dev/null")
        rate_positive_vertices<- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                 len = curve_length, threshold = threshold,
                                                                 cone_size = directions_per_cone, ball_radius = ball_radius,ball = ball, radius = radius)
        #sink()

        rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

        # calculate the TPR, FPR
        # To do the case where we consider true vertices to be close to any causal point, replace class_1, class_2 true vertices with the combined_true vertices
        # Otherwise, use class_1_true_Vertices, class_2_true_vertices
        if (class == 0)
        {
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        combined_true_vertices,combined_false_vertices))
          true_vertices = combined_true_vertices
          false_vertices = combined_false_vertices
        }
        else if(class == 1){
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        class_1_true_vertices,class_1_false_vertices))
          true_vertices = class_1_true_vertices
          false_vertices = class_1_false_vertices
        }
        else{
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        class_2_true_vertices,class_2_false_vertices))
          true_vertices = class_2_true_vertices
          false_vertices = class_2_false_vertices
        }
        previous_tpr_fpr = calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                             true_vertices,false_vertices)
        if (isTRUE((all.equal(c(1,1),previous_tpr_fpr)))){
          rate_ROC2 = matrix(1,ncol = 2, nrow = dim(total_rate_roc) - dim(rate_ROC)[1])
          rate_ROC = rbind(rate_ROC,rate_ROC2)
          break
        }
      }
    }
    else{
      for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated))){



        rate_positive_vertices<- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                 len = curve_length, threshold = threshold,
                                                                 cone_size = directions_per_cone, ball_radius = ball_radius,ball = ball, radius = radius)

        rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

        # calculate the TPR, FPR
        # To do the case where we consider true vertices to be close to any causal point, replace class_1, class_2 true vertices with the combined_true vertices
        # Otherwise, use class_1_true_Vertices, class_2_true_vertices
        if (class == 0)
        {
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        combined_true_vertices,combined_false_vertices))
          true_vertices = combined_true_vertices
          false_vertices = combined_false_vertices
        }
        else if(class == 1){
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        class_1_true_vertices,class_1_false_vertices))
          true_vertices = class_1_true_vertices
          false_vertices = class_1_false_vertices
        }
        else{
          rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                        class_2_true_vertices,class_2_false_vertices))
          true_vertices = class_2_true_vertices
          false_vertices = class_2_false_vertices
        }
        previous_tpr_fpr = calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                             true_vertices,false_vertices)
        if (isTRUE((all.equal(c(1,1),previous_tpr_fpr)))){
          rate_ROC2 = matrix(1,ncol = 2, nrow = dim(total_rate_roc) - dim(rate_ROC)[1])
          rate_ROC = rbind(rate_ROC,rate_ROC2)
          break
        }
      }
    }

    total_rate_roc = total_rate_roc + rate_ROC
    power = try(TPR_at_specified_FPR_metric(0.1,rate_ROC))
    print(power)
  }
  total_rate_roc = (total_rate_roc / counter)


  print(total_rate_roc)

  return(total_rate_roc)
}

######################################################################################
######################################################################################
######################################################################################
### Helper Functions ###

#' Returns a list that displays the relation of sub-level sets to vertices.
#'
#' @description Computes the vertex-sub-level set association.
#'
#' @export
#'@param dir (nx3 matrix):  The matrix of directions for which the (S/D) EC curve were computed over.
#'@param complex (list) : The list containing metadata of the Vertices, Edges, and Faces of the mesh (use process_off_file_v3 to obtain).
#'@param rate_vals (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#'@param ball_radius (float) : The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.
#'@param ball (boolean) : Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.
#'
#'@return vertex_list (list) : The association of vertices to sub-level sets, where the keys are the sub-level set indices, and the values are the vertex indices.
feature_vertex_association=function(dir,complex,len,ball_radius = 0, ball = FALSE){
  num_dir=dim(dir)[1]
  vertex_list=list()
  for(i in 1:num_dir){
    if (ball == FALSE){
      vertex_index=((i-1)*len+1):(i*len)
      projections <- complex$Vertices[,1:3]%*%dir[i,]
      buckets <- seq(min(projections),max(projections),length.out = len)
      step_length <- (max(projections) - min(projections))/len
      projection_buckets <- apply((projections - min(projections))/step_length,1, function(float) as.integer(float)) + len*(i-1)
      projection_buckets=projection_buckets+1-(len*(i-1))
      for (j in 1:len){
        feature_index=(i-1)*len+j
        associated_vertices=which(projection_buckets == j)
        vertex_list[[feature_index]]=associated_vertices
      }
    }
    else{
      vertex_index=((i-1)*len+1):(i*len)
      projections <- complex$Vertices[,1:3]%*%dir[i,]
      buckets <- seq(-ball_radius,ball_radius,length.out = len+1)

      #bucket these projections into curve_length number of groups; could have also solved this with the cut function
      step_length <- (2*ball_radius)/(len+1)
      projection_buckets <- apply((projections - min(projections))/step_length,1, function(float) as.integer(float)) + (len)*(i-1)
      projection_buckets=projection_buckets+1
      projection_buckets=projection_buckets+1-(len*(i-1))
      for (j in 1:len){
        feature_index=(i-1)*len+j
        associated_vertices=which(projection_buckets == j)
        vertex_list[[feature_index]]=associated_vertices
      }
    }
  }
  return(vertex_list)
}
