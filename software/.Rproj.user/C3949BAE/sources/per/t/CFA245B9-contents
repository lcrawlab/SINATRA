#' Try to understand how parameters of SINATRA affect reconstruction
#'
#' We work with a flat plane and a Gaussian bump on it
#'
#' We also test the reconstruction with a sphere and Gaussian bump on it.
#'
#' Implement Grid Search for this? - visualize the reconstructed regions side by side.

library(sinatra)

desired_num_cones <- 35
cap_radius <- 0.1
directions_per_cone <- 5


### Generate directions ###
dir <- generate_equidistributed_cones(desired_num_cones,cap_radius,directions_per_cone)

### Generate Data ###

nsim <- 50
curve_length <- 100
ball_radius <- 1.5
subdivision <- 3
ec_type <- 'ECT'

data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]))

for (i in 1:nsim){

  sphere1 = vcgSphere(subdivision = subdivision)
  sphere2 = vcgSphere(subdivision = subdivision)

  # Add noise to the sphere
  sphere1$vb[1:3,] = sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.02)
  sphere2$vb[1:3,] = sphere2$vb[1:3,]  * rnorm(dim(sphere2$vb)[2], mean = 1, sd = 0.02)

  sphere_mesh1 = convert_off_file(sphere1)
  sphere_mesh2 = convert_off_file(sphere2)

  ec_curve_class1 <- matrix(NA,nrow = 1,ncol=0)
  ec_curve_class2 <- matrix(NA,nrow = 1,ncol=0)

  ### compute EC curves for both classes of curves
  for (j in 1:dim(dir)[1]){

    vertex_function_class_1 <- sphere_mesh1$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
    vertex_function_class_2 <- sphere_mesh2$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])

    curve1 <- compute_standardized_ec_curve(sphere_mesh1, vertex_function_class_1, curve_length-1, first_column_index = FALSE,ball_radius)
    curve2 <- compute_standardized_ec_curve(sphere_mesh2, vertex_function_class_2, curve_length-1, first_column_index = FALSE,ball_radius)

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


### Run the model + select features with RATE
# how does bandwidth impact reconstruction?
rate_values <- find_rate_variables_with_other_sampling_methods(data,radius = 0, bandwidth = 0.1,
                                                               weights = TRUE, type = 'ESS')[,2]

plot(rate_values)
### Plot it back onto shape, and make rotating plot
sphere1 <- vcgSphere(subdivision = subdivision)
sphere1$vb[1:3,] <- sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.02)
complex <- convert_off_file(sphere1)

# reconstruct birth times of vertices
vert_matrix <- reconstruct_vertices_on_shape(dir, complex, rate_values, curve_length, cuts = length(rate_values),
                                             directions_per_cone, ball_radius, TRUE)

# define heatmap colors
color1='blue'
color2='lightgreen'
color3='orangered'
color3 = 'red'
col_pal=c(color1,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

# plot, using absolute birth times
# vert_heat1 <- colfunc(cuts)[vert_matrix1[,1]] #absolute
vert_heat1 = colfunc(1 + max(vert_matrix[,1]) - min(vert_matrix[,1]))[1 + vert_matrix[,1] - min(vert_matrix[,1])] # relative
plot3d(sphere1, col = vert_heat1, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')
rglwidget()
