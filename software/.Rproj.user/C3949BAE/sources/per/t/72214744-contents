#' Purpose of simulation" understand impact of different kernels in the power of SINATRA procedure.
#'
#' We are currently using the squared exponential / Gaussian kernel, but there are other ideas:
#' some suggestions say that using polynomial kernels might offer advantages for high dimensional classification problems
#' We could also use a more geometrically inclined kernel like in 'Gaussian Process Landmarking' - at least one that reflects the
#' EC-cone feature structure of our data.
#'
#' Perhaps other covariance functions like the Matern, or OU covariance functions can help.
#'
#' We should also try custom ones.
#'
#'
#' Simulation: measure performance via ROC curves on the easy medium difficult simulations?
#'
#' Produce plots of performance as a function of parameters

library(pracma)
library(sinatra)
library(ggplot2)
library(doParallel)

no_cores <- detectCores() - 10
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

n.simulations <- 2

simulation_results <- foreach(i=1:n.simulations, .combine = 'rbind') %:%
  foreach(j=c(0.1,1,10), .combine = 'rbind') %dopar% {
    set.seed(5*i+j)
    kernel_function <- function(X) polynomial_kernel(X, j, 2)

    res <- tryCatch( run_simulation_with_kernel(kernel_function, nsim = 10, curve_length = 25, grid_size = 25, distance_to_causal_point = 0.1,
                                                        causal_points = 10, shared_points = 10, num_cones = 40, eta = 0.1,
                                                        truncated = 500, ball = TRUE, ball_radius = 1.5,
                                                        min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 0,ec_type = 'ECT',
                                                        subdivision = 3, num_causal_region = 3, num_shared_region = 6),
                     error = function(x) {
                       return(matrix(nrow = 0,ncol = 3))
                     }
    )
    ### Label the results for each trial and directions ###
    rdf <- cbind(res, rep(j, dim(res)[1]) )
    rdf <- data.frame(cbind(rdf,rep(i,dim(rdf)[1])))
    rdf <- plyr::rename(rdf,c("X1" = "FPR","X2" = "TPR","X3" = "Class","X4" = "Index","X5" = "kernel_param","X6" = "Trial"))
  }


stopCluster(cl)

rdfmeans <- aggregate(simulation_results[c("FPR","TPR")],
                      by = list("Kernel.Param" = simulation_results$kernel_param,
                                "Index" = simulation_results$Index,
                                "Class" = simulation_results$Class), mean)

rdfmeans$Kernel.Param <- as.factor(rdfmeans$Kernel.Param)

### Plot results ###
ROC_curve_plt <- ggplot(data <- rdfmeans[rdfmeans$Class == 1,],aes(x = FPR, y = TPR, color = Kernel.Param)) +
  geom_line(stat = "identity") +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Kernel Parameter") +
  ggtitle(sprintf("Medium Simulation: Polynomial Kernel Deg 2")) +
  #geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)



########################################################
########################################################

run_simulation_with_kernel <- function(kernel_function, nsim = 10, curve_length = 25, grid_size = 25, distance_to_causal_point = 0.1,
                                       causal_points = 10, shared_points = 3, num_cones = 50, eta = 0.1,
                                       truncated = 300, ball = TRUE, ball_radius = 2.5,
                                       min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 0,ec_type = 'ECT',
                                       subdivision = 3, num_causal_region = 3, num_shared_region = 6){
  print("generating directions")
  print("generating data")
  # generate data

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


  num_cones <- dim(directions)[1]/directions_per_cone

  print("getting rate values")
  rate_values <- find_rate_variables_kernel(ec_curve_data, kernel_function, type = 'ESS')[,2]

  #Indices for Two Classes
  index1 = seq(1,nsim,2)
  complex_points1 = data$complex_points[index1]

  index2 = seq(2,nsim,2)
  complex_points2 = data$complex_points[index2]

  #Compute ROC using training data
  roc_curve1 =  compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                         curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                         eta = eta, directions_per_cone = directions_per_cone, directions = directions, class = 1,truncated = truncated,
                                         ball_radius = ball_radius, radius = radius, mode = 'sphere',subdivision = subdivision)
  roc_curve1 = cbind(roc_curve1, rep(1,dim(roc_curve1)[1]))
  roc_curve1 = cbind(roc_curve1,(1:dim(roc_curve1)[1]))

  roc_curve2 =  compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                         curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                         eta = eta, directions_per_cone = directions_per_cone, directions = directions,class = 2,truncated = truncated,
                                         ball_radius = ball_radius, radius = radius, mode = 'sphere', subdivision = subdivision)
  roc_curve2 = cbind(roc_curve2, rep(2,dim(roc_curve2)[1]))
  roc_curve2 = cbind(roc_curve2,(1:dim(roc_curve2)[1]))

  roc_curve = rbind(roc_curve1,roc_curve2)

  return(roc_curve)

}
########################################################
########################################################
########################################################

find_rate_variables_kernel <- function(gp_data, kernel, type = 'ESS'){
  n <- dim(gp_data)[1]
  X <- gp_data[,-1]
  gp_data[,1]=ifelse(gp_data[,1]>0,1,-1)
  #RATE
  f <- rep(0,n)
  Kn <- kernel(t(X))
  diag(Kn)=1

  # Get the samples of the latent posterior. The code for each of these methods can be found in 'GPC_Approximate_Inference'
  if ( type == 'Laplace' ){
    params = LaplaceApproximation(Kn,gp_data[,1])
    mu <- params[[1]]
    sigma <- params[[2]]
    fhat.samples = rmvnorm(1e4,mu,sigma)
  } else if ( type == 'EP' ){
    params = ExpectationPropagation(Kn,gp_data[,1])
    mu <- params[[1]]
    sigma <- params[[2]]
    fhat.samples = rmvnorm(1e4,mu,sigma)
  } else if ( type == 'ESS' ){
    fhat.samples = sinatra::Elliptical_Slice_Sampling(Kn,gp_data[,1],1e5,probit = TRUE)
  } else {
    stop(" Input one of 'Laplace','EP','ESS' as methods for Gaussian Process Inference ")
  }

  #is.positive.semi.definite(q_hat)
  # use RATE:
  cores= detectCores()
  res = RATE(X=X,f.draws=fhat.samples,prop.var = 1,snp.nms = colnames(X),cores = cores)
  rates = res$RATE
  return(cbind(1:length(rates),rates))
  }

########################################################
########################################################
########################################################
### Kernels ###

polynomial_kernel <- function(X, bias, exponent){
  n <- dim(X)[2]
  p <- dim(X)[1]

  K <- matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      K[i,j] <- (pracma::dot(X[,i],X[,j]) + bias)^exponent
    }
  }

  return(K)
}

OU_kernel <- function(X,length_scale){
  n <- dim(X)[2]
  p <- dim(X)[1]

  K <- matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        K[i,j] <- exp(-sqrt(sum((X[,i] - X[,j])^2))/length_scale)
      }
    }
  }

  return(K)
}

squared_exp_kernel <- function(X, bandwidth){
  n <- dim(X)[2]
  p <- dim(X)[1]

  K <- matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
        K[i,j] <- exp(-sum((X[,i] - X[,j])^2)/(2*p*bandwidth^2))
    }
  }
  return(K)
}

GaussKernel <- function(X,bandwidth){
  n <- dim(X)[2]
  p <- dim(X)[1]

  K <- matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        K[i,j] <- exp(-bandwidth*sum((X[,i] - X[,j])^2)/p)
      }
    }
  }
  return(K + t(K))
}

### custom kernels


