
######################################################################################
######################################################################################
######################################################################################

### Metric Curve Code ###

generate_metric_curve <- function(nsim = 20, curve_length = 20, grid_size = 30, distance_to_causal_point = 0.1, causal_points = 8,
                                  shared_points = 5, num_cones = 15, directions_per_cone = 4, eta = 0.1, ball_radius = 2.5, type){


  print("generating directions")
  # generate directions, length num_directions
  initial_cones <- 150
  directions <- generate_equidistributed_cones(initial_cones,0.1,directions_per_cone)
  num_cones = dim(directions)[1]/(directions_per_cone)


  print("generating data")
  # generate data
  data <- create_data_normal_fixed(num_sim = nsim, dir = directions, curve_length = curve_length,shared_points = shared_points,
                                   causal_points = causal_points,grid_size = grid_size,eta = eta,ball_radius = ball_radius)

  # Prune directions
  print("pruning data")
  temp <- prune_directions_to_desired_number(data = data$data[,-1],directions, initial_cones,curve_length,directions_per_cone,num_cones)
  directions <- temp[[1]]
  ec_curve_data <- temp[[2]]
  num_cones <- dim(directions)[1]/directions_per_cone

  metric_curve <- matrix(0, nrow = num_cones, ncol = 2)

  #want to vary directions, add on a new cone at each step
  for (i in 1:num_cones){
    print(sprintf("Analyzing Direction %i", i))

    # Generate the Rate using the data; take only data corresponding to desired directions
    print("computing rate values")
    rate_values <- find_rate_variables_with_other_sampling_methods(ec_curve_data[,1:(i*curve_length*directions_per_cone)],bandwidth = 0.01,type = 'ESS')[,2]

    print("computing metrics")
    if (type == 'vertex'){
      metrics <- compute_metrics_vertex(data_points = data$complex_points, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                        curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                        eta = eta, directions_per_cone = directions_per_cone, directions = directions[1:(i*directions_per_cone),], ball_radius = ball_radius,
                                        ball = ball)

    }
    if (type == 'feature'){
      metrics <- compute_metrics_feature(data_points = data$complex_points, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                         curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                         eta = eta, directions_per_cone = directions_per_cone, dir = directions[1:(i*directions_per_cone),], ball = ball,ball_radius = ball_radius,
                                         min_points = min_points)
    }
    if (type == 'cone'){
      metrics <- compute_metrics_cone(data_points = data$complex_points, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                      curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                      eta = eta, directions_per_cone = directions_per_cone, dir = directions[1:(i*directions_per_cone),],
                                      ball = ball, ball_radius = ball_radius,  min_points = min_points, radius = radius)
    }

    metric_curve[i,] <- metrics
  }

  return(metric_curve)
}


######################################################################################
######################################################################################
######################################################################################

# Computing Metrics
compute_metrics_vertex <- function(data_points,class_1_causal_points,class_2_causal_points,distance_to_causal_point = 0.1,
                                   rate_values,grid_size,eta = 0.1,directions_per_cone, curve_length,directions, ball_radius, ball = ball){

  num_vertices = grid_size^2
  #Initializing the aggregate ROC curve frame
  total_metric = matrix(0, nrow = length(data_points),ncol = 2)



  # go down the list of complexes?
  for (i in 1:length(data_points)){

    #Interpolating based on the causal and shared points in R^3 for each shape
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
    combined_true_vertices = union(class_1_true_vertices,class_2_true_vertices)

    class_1_false_vertices = setdiff(1:num_vertices, class_1_true_vertices)
    class_2_false_vertices = setdiff(1:num_vertices, class_2_true_vertices)
    combined_false_vertices = setdiff(1:num_vertices, combined_true_vertices)

    rate_ROC <- matrix(0,nrow = 1,ncol = 2)
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = floor(length(rate_values)/5)) ) ){

      rate_positive_vertices <- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                len = curve_length, threshold = threshold,
                                                                cone_size = directions_per_cone, ball_radius = ball_radius)

      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

      TPR_FPR <- calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                   class_1_true_vertices,class_1_false_vertices)
      rate_ROC <- rbind(rate_ROC, TPR_FPR)

      # if FPR > 0.05, break
      if(TPR_FPR[1] > 0.05) break();
    }

    metrics <- c(TPR_at_specified_FPR_metric(0.05,rate_ROC),
                 size_of_intersection_metric(combined_true_vertices, rate_positive_vertices))

    total_metric[i,] <- metrics
  }

  averaged_metrics <- colMeans(total_metric)

  return(averaged_metrics)
}

compute_metrics_feature <- function(data_points,class_1_causal_points,class_2_causal_points,distance_to_causal_point = 0.1,
                                    rate_values,grid_size,eta = 0.1,directions_per_cone, curve_length,
                                    ball = TRUE, ball_radius = ball_radius, dir, min_points = 2){
  num_vertices = grid_size^2
  #Initializing the aggregate ROC curve frame
  total_metric = matrix(0, nrow = length(data_points),ncol = 2)
  for (j in seq(1,length(data_points),2)){
    #Interpolating based on the causal and shared points in R^3 for each shape
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
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values)))){
      selected_features = matrix(ordered_rate_values[(which(ordered_rate_values[,2]>=threshold)),],ncol=2)[,1]

      #rate_positive_features = c(empty_list,intersect(positive_features, selected_features))
      #rate_negative_features = c(empty_list,intersect(negative_features, selected_features))
      rate_positive_features = selected_features
      rate_negative_features = setdiff(1:length(rate_values), selected_features)
      TPR_FPR <- calculate_TPR_FPR(rate_positive_features,rate_negative_features,
                                   positive_features,negative_features)
      rate_ROC <- rbind(rate_ROC, TPR_FPR)

      # if FPR > 0.05, break
      if(TPR_FPR[1] > 0.05) break();
    }
    if (NaN %in% rate_ROC){
      total_metric = total_metric[-((j+1)/2),]
      next
    }
    else{
      metrics <- c(TPR_at_specified_FPR_metric(0.05,rate_ROC))
      total_metric[((j+1)/2),] <- c(metrics,0)
    }
  }

  averaged_metrics <- colMeans(total_metric)

  return(averaged_metrics)
}

compute_metrics_cone <- function(data_points,class_1_causal_points,class_2_causal_points,distance_to_causal_point = 0.1,
                                 rate_values,grid_size,eta = 0.1,directions_per_cone, curve_length,dir,
                                 ball = ball, ball_radius = ball_radius, min_points = 2, radius = 2){
  num_vertices = grid_size^2
  total_metric = matrix(0, nrow = length(data_points)/2,ncol = 2)
  remove = c()
  for (j in seq(1,length(data_points),2)){
    #Interpolating based on the causal and shared points in R^3 for each shape
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
    combined_true_vertices = union(class_1_true_vertices,class_2_true_vertices)
    combined_false_vertices = setdiff(1:num_vertices, combined_true_vertices)

    rate_ROC <- matrix(0,nrow = 1,ncol = 2)

    feature_with_index  = cbind(1:length(rate_values),rate_values)
    ordered_rate_values = feature_with_index[order(feature_with_index[,2], decreasing = TRUE),]
    vertex_list1 = feature_vertex_association(dir = dir, complex = complex1, len = curve_length, ball_radius = ball_radius, ball = ball)
    vertex_list2 = feature_vertex_association(dir = dir, complex = complex2, len = curve_length, ball_radius = ball_radius, ball = ball)

    total_reconstructions = curve_length * (dim(dir)[1]/directions_per_cone)
    positive_reconstructions = c()
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
        reconstruction_vector = rep(0,curve_length*dim(dir)[1])
        reconstruction_number = i*l
        reconstruction_index  = constant + l -1
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
          reconstruction_vector[radius_points] = 1
          reconstruction_features[[reconstruction_number]][[k]] = radius_points
        }
        reconstruction_vertices1[[reconstruction_number]] = summarize_vertices(dir = cone_dirs,complex = complex1,rate_vals = reconstruction_vector,len = curve_length,
                                                                               threshold = 0.5,cone_size = directions_per_cone,ball_radius = ball_radius,ball = ball)
        reconstruction_vertices2[[reconstruction_number]] = summarize_vertices(dir = cone_dirs,complex = complex2,rate_vals = reconstruction_vector,len = curve_length,
                                                                               threshold = 0.5,cone_size = directions_per_cone,ball_radius = ball_radius,ball = ball)
      }
    }
    for (i in 1:total_reconstructions){
      found_vertices1 = intersect(reconstruction_vertices1[[i]],class_1_true_vertices)
      found_vertices2 = intersect(reconstruction_vertices2[[i]],class_2_true_vertices)
      if ((length(found_vertices1)+length(found_vertices2))>=min_points){
        positive_reconstructions = c(positive_reconstructions, i)
      }
    }
    empty_list = c()
    negative_reconstructions <- setdiff(1:total_reconstructions,positive_reconstructions)
    if (length(positive_reconstructions) == 0 || length(negative_reconstructions) == 0){
      remove = c(remove,((j+1)/2))
      next
    }
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values)))){
      selected_features = matrix(ordered_rate_values[(which(ordered_rate_values[,2]>=threshold)),],ncol=2)[,1]

      rate_positive_reconstructions = c()

      for (l in 1:total_reconstructions){
        if (check_reconstruction(selected_features =selected_features, reconstruction_features_index = reconstruction_features[[l]])){
          rate_positive_reconstructions = c(rate_positive_reconstructions,l)
        }
      }
      rate_negative_reconstructions = setdiff(1:total_reconstructions, rate_positive_reconstructions)
      TPR_FPR <- calculate_TPR_FPR(rate_positive_reconstructions,rate_negative_reconstructions,
                                   positive_reconstructions,negative_reconstructions)
      rate_ROC <- rbind(rate_ROC, TPR_FPR)
      if(TPR_FPR[1] > 0.05) break();
    }
    if (NaN %in% rate_ROC){
      remove = c(remove,((j+1)/2))
      next
    }
    else{
      metrics <- c(TPR_at_specified_FPR_metric(0.05,rate_ROC))
      total_metric[((j+1)/2),] <- c(metrics,0)
    }
  }
  if (length(remove) == dim(total_metric)[1]){
    print('no good directions')
    averaged_metrics <- colMeans(total_metric)
  }
  if (length(remove) == 0){
    averaged_metrics <- colMeans(total_metric)
  }
  else{
    total_metric = matrix(total_metric[-remove,],ncol = 2)
    averaged_metrics <- colMeans(total_metric)
  }


  return(averaged_metrics)
}

######################################################################################
######################################################################################
######################################################################################

### Helper Functions ###
#' TPR/FPR
#' @export
calculate_TPR_FPR <- function(positive_vertices, negative_vertices, true_vertices, false_vertices){
  TP <- length(intersect(positive_vertices, true_vertices))
  FP <- length(positive_vertices) - TP
  TN <- length(intersect(negative_vertices, false_vertices))
  FN <- length(negative_vertices) - TN

  TPR = TP/(TP + FN)
  FPR = FP/(FP + TN)

  return(c(FPR,TPR))
}

#Input an ROC curve: a matrix of size n x 2, where n is the threshold size.
# Find the (FPR,TPR) pair such that FPR < specified FPR.
TPR_at_specified_FPR_metric <- function(FPR, ROC){
  sorted_ROC <- ROC[order(ROC[,1]),]
  desired_index <- which(sorted_ROC[,1] >= FPR)[1]
  return(sorted_ROC[desired_index,2])
}

# This is extremely similar to a TPR; intersection over total union.
# Inputs are lists of vertex indices, 1 to n. Take causal points to be the union of class1,class2 causal points;
# selected points should be the output of the rate & cone recontruction idea.
size_of_intersection_metric <- function(causal_points,selected_points){
  int <- intersect(causal_points,selected_points)
  union <- union(causal_points,selected_points)
  return(length(int)/length(union))
}
