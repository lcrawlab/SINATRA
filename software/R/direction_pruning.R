######## Direction Pruning Code ########

# assume that the cones include their central direction. Outputs the inter-class correlations and intra-class correlations
compute_total_correlations <- function(data,num_cones, curve_length, dir_per_cone){
  median_correlations <- 1:num_cones
  for (i in 0:(num_cones-1)){
    # test the correlations of only the central directions
    central_dir_index <- i*curve_length*dir_per_cone
    ec_curves <- data[,(central_dir_index+1):(central_dir_index + curve_length)]
    median_correlations[i+1] <- median(cor(t(ec_curves)))
  }
  return(median_correlations)
}

# input which cone to get- the function then computes the correlation matrix, both the interclass and intraclass correlation matrices
# works with central cone direction but can be easily modified.
compute_cone_class_correlations <- function(data, cone, curve_length, dir_per_cone){
  central_dir_index <- cone*curve_length*dir_per_cone
  num_data <- dim(data)[1]
  class_1_ec_curves <- data[seq(1,num_data,2),(central_dir_index+1):(central_dir_index + curve_length)]
  class_2_ec_curves <- data[seq(2,num_data,2),(central_dir_index+1):(central_dir_index + curve_length)]


  class_1_correlations <- cor(t(class_1_ec_curves))
  class_2_correlations <- cor(t(class_2_ec_curves))
  interclass_correlations <- cor(t(class_1_ec_curves),t(class_2_ec_curves))

  return(list(class1 = class_1_correlations, class2 = class_2_correlations, inter = interclass_correlations))
}

prune_directions_to_desired_number <- function(data, directions, num_cones, curve_length, dir_per_cone,desired_number){
  cors <- compute_total_correlations(data, num_cones, curve_length,dir_per_cone)

  # get the desired number of best cones by correlation, by removing the ones with high correlation
  idxs <- order(cors)[(desired_number+1):num_cones]
  idxs <- (idxs - 1)*dir_per_cone + 1

  return(update_data_and_directions(idxs,data,directions,curve_length,dir_per_cone))
}

# prune directions + data with correlations greater than 0.98, say.
# We might want to keep only a certain number; pick the k smallest ones?
prune_directions <- function(data, directions, num_cones, curve_length, dir_per_cone){
  cors <- compute_total_correlations(data, num_cones, curve_length,dir_per_cone)

  # get the central directions to remove, and then the associated directions in the cone
  idxs <- which(cors > 0.98)
  idxs <- (idxs - 1)*dir_per_cone + 1

  return(update_data_and_directions(idxs,data,directions,curve_length,dir_per_cone))
}

# prune directions that repeat already observed information, along with low variance directions - use interclass/intraclass metrics
# rank the directions with the most intraclass variance, least interclass variance
prune_low_var_and_collinear_directions <- function(data, directions, num_cones,curve_length, dir_per_cone){
  # rank the directions by least variance in a class, and most variance between classes. We can summarize the
  # associated correlation matrices by some statistic?



  # add variables that aren't collinear. Work with the central directions
}

# Helper function for updating data / directions, given the central cone indices to prune
update_data_and_directions <- function(idxs, data, directions, curve_length,dir_per_cone){
  # get the associated directions by shifting the idxs of the central directions
  direction_idxs <- idxs
  for (i in 1:(dir_per_cone-1)){
    direction_idxs <- c(direction_idxs, idxs + i)
  }
  # remove the directions associated with each of these indices
  pruned_directions <- directions[-direction_idxs,]

  # Remove the data corresponding to these cones as well
  temp_idxs <- (direction_idxs - 1)*curve_length+1
  data_idxs <- temp_idxs
  for (i in 1:(curve_length-1)){
    data_idxs <- c(data_idxs, temp_idxs + i)
  }
  pruned_data <- data[,-data_idxs]

  return(list(pruned_directions,pruned_data))
}



find_directions_with_power <- function(runs = 1, nsim = 50, curve_length = 10, grid_size = 25, distance_to_causal_point = 0.1,
                                       causal_points = 10,shared_points = 3, num_directions = 10, eta = 0.1,
                                       truncated = FALSE, two_curves = TRUE, ball = TRUE, ball_radius = 1, type = 'feature',
                                       min_points = 3,directions_per_cone = 5, cap_radius = 0.15, radius = 1){
  # Generate the random cones
  set.seed(1230)
  total_directions = generate_equidistributed_points(num_directions)
  dir_powers = cbind(total_directions,rep(0,dim(total_directions)[1]))
  dir_powers = cbind(dir_powers,rep(0,dim(total_directions)[1]))
  dir_powers = cbind(dir_powers,rep(0,dim(total_directions)[1]))
  dir_powers = cbind(dir_powers,rep(-1,dim(total_directions)[1]))

  print(dim(total_directions))
  for (i in 1:(dim(total_directions)[1])){
    print(paste("Onto Direction", i))
    total_cor = c()
    class1_cor = c()
    class2_cor = c()
    directions <- rodriq(total_directions[i,],cap_radius,directions_per_cone)

    # generate data
    num_vertices <- grid_size^2
    for (j in 1:runs){
      data <- create_data_normal_fixed(num_sim = nsim,dir = directions,curve_length = curve_length,shared_points = shared_points,
                                       causal_points = causal_points,grid_size = grid_size,eta = eta,ball = ball, ball_radius = ball_radius)

      #Checking if RATE runs on the direction, for just computing correlations, we can just let rate_values be any value.
      #rate_values <- try(find_rate_variables_with_other_sampling_methods(data$data, bandwidth = 0.01,type = 'ESS')[,2])
      rate_values <- 3
      if (inherits(rate_values,'try-error')){
        dir_powers[i,4] = -1
        dir_powers[i,5] = -1
        dir_powers[i,6] = -1
        break
      }
      else{
        cors =  median(cor(t(data$data[,-1])))
        #Indices for Two Clases
        index1 = seq(1,2*nsim,2)
        complex_data1 = data$data[index1,-1]

        index2 = seq(2,2*nsim,2)
        complex_data2 = data$data[index2,-1]

        class1_cor = c(class1_cor, median(cor(t(complex_data1))))
        class2_cor = c(class2_cor, median(cor(t(complex_data2))))
        total_cor  = c(total_cor, cors)
        next
        # If we want to assess accuracy too, remove the next.
        roc_curve = try(compute_roc_curve_cones(data = data$complex_points, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                                curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                                eta = eta, directions_per_cone = directions_per_cone, class = 0, truncated = truncated,
                                                ball = ball, ball_radius = ball_radius,
                                                dir = directions,  min_points = min_points, radius = radius))
        if (inherits(roc_curve,'try-error')){
          next
        }
        else{
          power = try(TPR_at_specified_FPR_metric(0.1,roc_curve))
          if (inherits(power,'try-error')){
            next
          }
          else{
            dir_powers[i,7] = power
          }
        }
      }
    }
    dir_powers[i,4] = median(class1_cor)
    dir_powers[i,5] = median(class2_cor)
    dir_powers[i,6] = median(total_cor)
  }
  return(dir_powers)
}
