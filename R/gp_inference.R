####### RATE Code ######

#' Gaussian Kernel Implementation
#'
#' @description \code{GaussKernel} computes the covariance matrix given by the Gaussian kernel below.
#'
#' @param X (p x n matrix): the design matrix where columns are observations
#' @param bandwidth (float): free parameter for the Gaussian kernel
#'
GaussKernel <- function(X,bandwidth){
  n <- dim(X)[2]
  p <- dim(X)[1]

  K <- matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
        if(j > i){
          K[i,j] <- exp(-bandwidth*sum((X[,i] - X[,j])^2)/p)
        }
    }
  }
  K = K + t(K)
  diag(K) = 1
  return(K)
}


#'
#' Computes the gradient of the Gauss kernel with respect to the bandwidth
gradGaussKernel <- function(X,bandwidth){
  n <- dim(X)[2]
  p <- dim(X)[1]

  K <- matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      if(j > i){
        K[i,j] <- -exp(-bandwidth*sum((X[,i] - X[,j])^2)/p)*sum((X[,i] - X[,j])^2)/p
      }
    }
  }

  grad_cov <- K + t(K)
  diag(grad_cov) <- 0

  return(grad_cov)
}

#' Gaussian Kernel Implementation
#'
#' @description \code{UnstandardizedGaussKernel} computes the covariance matrix given by the Gaussian kernel below, without
#' setting the diagonal to zero.
#'
#' @param X (p x n matrix): the design matrix where columns are observations
#' @param bandwidth (float): free parameter for the Gaussian kernel
#'
SqExpKernel <- function(X,bandwidth){
  n <- dim(X)[2]
  p <- dim(X)[1]

  K <- matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      K[i,j] <- exp(-sum((X[,i] - X[,j])^2)/(2*bandwidth^2))
    }
  }

  return(K)
}

#' Gaussian Kernel Implementation
#'
#' @description \code{gradUnstandardizedGaussKernel} computes the gradient covariance matrix given by the Gaussian kernel below, without
#' setting the diagonal to zero.
#'
#' @param X (p x n matrix): the design matrix where columns are observations
#' @param bandwidth (float): free parameter for the Gaussian kernel
#'
gradSqExpKernel <- function(X,bandwidth){
  n <- dim(X)[2]
  p <- dim(X)[1]

  K <- matrix(0,nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      K[i,j] <- exp(-sum((X[,i] - X[,j])^2)/(2*bandwidth^2))*(sum((X[,i] - X[,j])^2)/(bandwidth^3))
    }
  }

  return(K)
}



#' Derive Rate Values
#'
#' @export
#'
#' @import mvtnorm
#' @import parallel
#'
#' @description \code{find_rate_variables_with_other_sampling_methods} returns a vector of variable importances using RATE.
#' We fit a GPC classifier to the data, draw samples from the latent posterior using one of Lapace Approximation, Elliptical Slice Samping,
#' or Expectation Propogation. After posterior inference, we fit RATE to derive association measure for each sub-level set of the design matrix.
#'
#' @param gp_data (matrix) : The design matrix of (S/D) EC curves and the associated class labels
#' @param bandwidth (float) : The bandwidth of the Gaussian Kernel used to fit the GPC.
#' @param type (string) : The sampling method used. We currently support Laplace's method, Elliptical Slice Sampling, and Expectation Propogation.
#'
#' @return rate_values (nx2 matrix) : The derived variable importance values and the row number denoting which sub-level set it corresponds to.
find_rate_variables_with_other_sampling_methods <- function(gp_data,bandwidth = 0.01,type = 'Laplace'){
  n <- dim(gp_data)[1]
  X <- gp_data[,-1]
  gp_data[,1]=ifelse(gp_data[,1]>0,1,-1)
  h <- bandwidth #median(dist(X)) #need a better choice of this; how does bandwidth affect kernel choice?
  #RATE
  f <- rep(0,n)
  Kn <- GaussKernel(t(X),1/(2*h^2))

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
    fhat.samples = Elliptical_Slice_Sampling(Kn,gp_data[,1],1e5,probit = TRUE)
  } else {
    stop(" Input one of 'Laplace','EP','ESS' as methods for Gaussian Process Inference ")
  }

  #is.positive.semi.definite(q_hat)
  # use RATE:
  print('Running RATE')
  cores= detectCores()
  #  cores = 8
  res = RATE(X=X,f.draws=fhat.samples,prop.var = 1,snp.nms = colnames(X),cores = cores)
  rates = res$RATE
  return(cbind(1:length(rates),rates))
}


###### Gaussian Approximation to Posterior of Latent Values ######

##############
#' Use Expectation Propogation to Approximate mean & covariance
#' @description \code{ExpectationPropogation} Approximates the latent posterior with a Gaussian distributions; it does so by moment matching.
#' Pseudocode taken from Rasmussen and Williams, Chapter 3. This function outputs the mean and covariance of the
#' approximated posterior. To actually generate samples from the latent posterior, generate samples from a multivariate
#' normal with the parameters returned by this function.
#'
#' @export
#' @import mvtnorm
#'
#' @param K (matrix): the covariance matrix for the GP model
#' @param class_labels (vector): +/- 1 values indicating the class labels of the data points
#'
#' @return params (list): list of the posterior mean and variances.
#'
ExpectationPropagation <- function(K, class_labels){
  n <- length(class_labels)
  nu_tilde <- matrix(0,nrow = n,ncol = 1)
  tau_tilde <- matrix(0,nrow = n,ncol = 1)
  sigma <- K
  mu <- matrix(0,nrow = n,ncol = 1)

  #initialize vectors
  nu_minus<- matrix(0,nrow = n,ncol = 1)
  tau_minus <- matrix(0,nrow = n,ncol = 1)
  mu_minus<- matrix(0,nrow = n,ncol = 1)
  sigma_minus <- matrix(0,nrow = n,ncol = 1)
  mu_hat <- matrix(0,nrow = n,ncol = 1)
  sigma_hat <- matrix(0,nrow = n,ncol = 1)


  # change this stopping condition
  for(j in 1:100){
    for (i in 1:n){
      tau_minus[i] <- sigma[i,i]^-2 - tau_tilde[i]
      nu_minus[i] <- (sigma[i,i]^-2)*mu[i] - nu_tilde[i]

      mu_minus[i] <- nu_minus[i]/tau_minus[i]
      sigma_minus[i] <- tau_minus[i]^-0.5

      # Compute Marginal Moments
      z_i <- class_labels[i]*mu_minus[i]/(sqrt(1+sigma_minus[i]^2))
      mu_hat[i] <- mu_minus[i] + (class_labels[i] * sigma_minus[i]^2 * stats::dnorm(z_i) )/(stats::pnorm(z_i)*sqrt(1+sigma_minus[i]^2))
      sigma_hat[i] <- sqrt( sigma_minus[i]^2 - (sigma_minus[i]^4*stats::dnorm(z_i))/((1+sigma_minus[i]^2)*stats::pnorm(z_i))*(z_i + stats::dnorm(z_i)/stats::pnorm(z_i)) )

      # Update Site Parameters
      delta_tau_tilde <- sigma_hat[i]^-2 - tau_minus[i] - tau_tilde[i]
      tau_tilde[i] <- tau_tilde[i] + delta_tau_tilde
      nu_tilde[i] <- (sigma_hat[i]^-2)*mu_hat[i] - nu_minus[i]

      # Update Sigma, mu - the parameters of the posterior
      sigma <- sigma - (( delta_tau_tilde^-1 + sigma[i,i])^-1)*sigma[,i]%*%t(sigma[,i])
      mu <- sigma%*%nu_tilde
    }

    #Recompute posterior parameters
    S_tilde <- diag(as.vector(tau_tilde))
    B <- diag(n) + sqrt(S_tilde)%*%K%*%sqrt(S_tilde)
    L <- chol(B)
    V <- solve(t(L),sqrt(S_tilde)%*%K)
    sigma <- K - t(V)%*%V
    mu <- sigma%*%nu_tilde
  }
  # Compute Log Marginal Likelihood
  B <- diag(n) + sqrt(S_tilde) %*% K %*% sqrt(S_tilde)
  T <- diag(as.vector(tau_minus))
  term1 <- 0.5*sum(log(1 + tau_tilde/tau_minus)) - sum(log(diag(L)))
  term2 <- 0.5*sum(log(stats::pnorm(class_labels*mu_minus)/sqrt(1 + sigma_minus^2)))
  term3 <- 0.5*t(nu_tilde)%*%(K - K %*% sqrt(S_tilde) %*% solve(B) %*% sqrt(S_tilde) %*% K - (T + solve(S_tilde)) ) %*% (nu_tilde)
  term4 <- 0.5* t(mu_minus)%*% T %*% solve(S_tilde + T) %*% (S_tilde %*% mu_minus - 2*nu_tilde)

  log_Z <- term1 + term2 + term3 + term4


  params <- list(mu,sigma, nu_tilde, tau_tilde, log_Z)
  return(params)
}

##############
#' Use Laplace Approximation to Approximate mean & covariance
#'
#' @export
#'
#' @description \code{LaplaceApproximation}
#' Approximates the latent posterior with a Gaussian distribution; it does so by finding the mode of the posterior, and
#' using the Hessian (second order Taylor expansion) as an approximation of the covariance. Newton Raphson is used to find
#' the mode of the posterior.
#'
#' Pseudocode taken from Rasmussen and Williams, Chapter 3. This function outputs the mean and covariance of the
#' approximated posterior. To actually generate samples from the latent posterior, generate samples from a multivariate
#' normal with the parameters returned by this function.
#'
#'
#'  @param Kn (matrix): the covariance matrix for the GP model
#'  @param class_labels (vector): +/- 1 values indicating the class labels of the data points
#'
#' @return params (list): list of the posterior mean and variances.


LaplaceApproximation <- function(Kn, class_labels){
  f <- rep(0,length(class_labels))
  #Change stopping conditions for convergence
  # do Newton IRLS procedure to use Laplace / Gaussian approximation.
  for(k in 1:1000){
    W <- diag(as.vector(sigmoid(f)*(1-sigmoid(f))))
    B <- diag(x = 1,length(class_labels)) + sqrt(W) %*% Kn %*% sqrt(W)
    #Kinda show
    L <- chol(B)
    b <- W%*%f + (class_labels+1)/2 - sigmoid(f)
    #Kinda slow
    #L has components that
    #Try to not use all the solves
    a <- b - solve(sqrt(W)%*%t(L),solve(L,sqrt(W)%*%Kn%*%b))
    f <- Kn%*%a
  }
  v = solve(L, sqrt(W)%*%Kn)

  mean <- f
  covariance <- Kn - t(v)%*%v
  params <- list(mean, covariance)

}

####### MCMC Inference for Posterior of Latent Values #######

#'  Draw samples from posterior using Elliptical Slice Sampling
#'
#'  @export
#'
#'  @import FastGP
#'
#'   @description  \code{Elliptical_Slice_Sampling} Based on Iain Murray's paper 'Elliptical Slice Sampling'. Implemented using the FastGP package. The function returns
#' the desired number of mcmc samples
#'  @param K (matrix): the covariance matrix of the GP model
#'  @param class_labels (vector): the class labels for each data point, +/- 1.
#'  @param num_mcmc_samples (int): the number of desired mcmc samples to be returned
#'  @param probit (boolean): set TRUE if the link function in the model is probit; otherwise the function uses the logistic link.
#'
#'  @return samples (vector): Vector of Samples obtained from ESS sampling.
Elliptical_Slice_Sampling <- function(K,class_labels,num_mcmc_samples, probit = TRUE){
  if(probit){
    samples <- FastGP::ess(probit_log_likelihood, class_labels,K, num_mcmc_samples,
                           burn_in = 1000, N = length(class_labels), TRUE)
    return(samples)

  } else {
    samples <- FastGP::ess(logistic_log_likelihood, class_labels,K, num_mcmc_samples,
                           burn_in = 1000, N = length(class_labels), TRUE)
    return(samples)
  }

}

############################################################
############################################################
############################################################
##### Model Selection Code #####

#' Compute Marginal Likelihood Gradient
#'
#' @export
#'
#' @description \code{compute_marginal_likelihood_gradient} Computes Marginal likeilhood Gradient for the GPC Model, when using Expectation Prop.
#' Used for as model evidence in order to compute optimal parameters/ hyperparameters for the model. Implemented so far
#' only for the GaussKernel kernel function. Might be expanded in the future.
#'
#' @param X (matrix) : The design matrix of (S/D) EC curves
#' @param kernel_param (float) : Bandwidth parameter for the Gaussian Kernel
#' @param class_labels (vector) : vector of +1,-1 of class labels.
#'
#' @return grad_log_Z (float): the gradient log likelihood given the parameters.
compute_marginal_likelihood_gradient <- function(X, kernel_param, class_labels){
  n <- length(class_labels)
  K <- SqExpKernel(t(X), kernel_param)
  K <- K + 0.001*diag(length(class_labels))
  params <- ExpectationPropagation(K, class_labels)
  nu_tilde <- params[[3]]
  tau_tilde <- params[[4]]

  # compute gradient
  S_tilde = diag(as.vector(tau_tilde))
  L <- chol(diag(n) + sqrt(S_tilde) %*% K %*% sqrt(S_tilde))
  b <- nu_tilde - sqrt(S_tilde) %*% solve(L) %*% solve(t(L)) %*% sqrt(S_tilde) %*% K %*% nu_tilde
  R <- b%*%t(b) - sqrt(S_tilde) %*% solve(t(L)) %*% solve(L) %*% sqrt(S_tilde)
  C <- gradSqExpKernel(t(X), kernel_param)

  grad_log_Z <- 0.5*sum(diag(R %*% C))

  return(grad_log_Z)
}

#' Compute Marginal Likelihood
#'
#' @export
#'
#' @description \code{compute_marginal_likelihood} Computes Marginal likeilhood for the GPC Model, when using Expectation Prop.
#' Used for as model evidence in order to compute optimal parameters/ hyperparameters for the model. Implemented so far
#' only for the GaussKernel kernel function. Might be expanded in the future.
#'
#' @param X (matrix) : The design matrix of (S/D) EC curves
#' @param kernel_param (float) : Bandwidth parameter for the Gaussian Kernel
#' @param class_labels (vector) : vector of +1,-1 of class labels.
#'
#' @return log_Z (float): the EP log likelihood given the parameters.
compute_marginal_likelihood <- function(X, kernel_param, class_labels){
  K <- SqExpKernel(t(X), kernel_param)
  K <- K + 0.001*diag(length(class_labels))
  params <- ExpectationPropagation(K, class_labels)
  log_Z <- params[[5]]

  return(log_Z)
}

#' Pick best kernel parameters

#' @export
#'
#' @description \code{optimize_params_EP} Uses nonlinear conjugate gradient
#' optimizer to pick out the best value of the kernel params
#'
#' @param X (matrix) : The design matrix of (S/D) EC curves
#' @param kernel_param (float) : Bandwidth parameter for the Gaussian Kernel
#' @param class_labels (vector) : vector of +1,-1 of class labels.
#'
#' @return Rcgmin output
optimize_params_EP <- function(X, class_labels, initial_params){
  fn <- function(param) compute_marginal_likelihood(X,param,class_labels)
  gr <- function(param) compute_marginal_likelihood_gradient(X,param,class_labels)

  return(Rcgmin::Rcgmin(initial_params, fn, gr))
}


############################################################
############################################################
############################################################
##### Helper Functions #####

#'Sigmoid Transformation
#'
#'
#' @export
#'
#'@description \code{sigmoid} Applies the sigmoid transformation
#'
#'@param x (vector): vector of values to apply the sigmoid transformation
#'
#'@return  sigmoid_x (vector): vector of values post sigmoid transformation
sigmoid <- function(x){
  return(1/(1+exp(-x)))
}

#' Probit Log Likelihood
#' @export
probit_log_likelihood <- function(latent_variables, class_labels){
  return(sum(log(stats::pnorm(latent_variables*class_labels))))
}

#' Probit Log Likelihood
#' @export
#' @description When the link function is logistic
logistic_log_likelihood <- function(latent_variables, class_labels){
  return(-sum(log(1+exp(latent_variables*class_labels))))
}

#### Other Feature Selection Functions ####

#### Bayesian Variable Selection ####

#' Conducts the Bayesian Feature Selection
#'
#' @description \code{find_bayesian_variables} returns a vector of the selected features/indices of the curve
#'
#' The function uses the package varbvs to conduct feature selection. See Carbonetto and Stephens 2012.
#' Once the model is fitted, we set the features who's pip (posterior inclusion probability) is greater than a user specified cutoff (typically 0.1,0.5)
#' These features are returned along with a window of neighboring features if the user wishes (typically recommended).
#'
#' @param data (nxm matrix) A nxm matrix containing the covariates and responses.
#' @param param (float (0-1)) The cutoff for PIP for including features.
#' @param radius (positive integer ). An integer parameter specifying if, and how many neighboring features should be considered for feature selection.
#' This is done to capture critical points that may be 'close' to the selected feature.
#' @return The output is a vector of indices/features to be selected.
#'
find_bayesian_variables=function(data,param=0.5,radius=0,weights=FALSE){
  #Convert the response to 0 or 1, for binary classification in the varbvs model
  data[,1]=ifelse(data[,1]>0,1,0)
  #Fit the model
  fit=varbvs(X = data[,-1],Z = NULL, y= data[,1],family='binomial',verbose = FALSE)
  #Find the pips for each feature.
  bayesian_probs=fit$pip
  #Only keep the features that are above the threshold (pip).
  if (weights==TRUE){
    return(cbind(1:length(bayesian_probs),abs(bayesian_probs)))
  }
  want_indices_bayesian=which(param<bayesian_probs)
  #If the radius is not zero, we include the 'n' neighboring features as well.
  real_indices=c()
  for (i in 1:(length(want_indices_bayesian))){
    for (j in 0:radius){
      real_indices=c(real_indices,want_indices_bayesian[i]+j)
      real_indices=c(real_indices,want_indices_bayesian[i]-j)
    }
  }
  real_indices=unique(real_indices)
  return(real_indices)
}

#### Elastic Net ####
#' Conducts the Lasso Feature Selection
#'
#' @export
#' @description \code{find_elastic_variables} returns a vector of the selected features/indices of the curve
#'
#' The function uses Elastic Net to conduct feature selection. A cross validated elastic net model is fit on the binary response data,
#' and the feature importances are extracted.
#' Similar to lasso, the selected features are those that are non-zero.
#' These features are returned along with a window of neighboring features if the user wishes (typically recommended).
#'
#' @param data (nxm matrix) A nxm matrix containing the covariates and responses.
#' @param radius (positive integer ). An integer parameter specifying if, and how many neighboring features should be considered for feature selection.
#' This is done to capture critical points that may be 'close' to the selected feature.
#' @return The output is a vector of indices/features to be selected.
#'
find_elastic_variables=function(data,radius=0,weights=FALSE, grouped_data = TRUE, alpha = 0.5){
  #Fit the model
  #browser()

  # need intercept = True
  if(grouped_data == FALSE) {
    regression_model=glmnet::cv.glmnet(data[,-1], data[,1], alpha = alpha,intercept = TRUE,family='binomial') # need groups added here
    #Extract the coefficients
    tmp_coeffs = as.matrix(stats::coef(regression_model, s = "lambda.min"))
    regression_coeff= as.vector(tmp_coeffs)
  } else {
    groups <- rep(1:(dim(data)[2]/3),each = 3)
    model <- gglasso::cv.gglasso(x = data[,-1], y = 2*data[,1] - 1, group = groups,
                                 intercept = TRUE, loss = "logit", nlambda = 500, pred.loss = "loss",nfolds = 10)
    # cv with misclassification error doesn't work?

    #browser()
    # plot(lasso)
    plot(model)
    regression_coeff <- as.matrix(stats::coef(model, s = 'lambda.min'))

  }

  #Don't include the intercept
  regression_coeff=regression_coeff[2:length(regression_coeff)]

  if (weights==TRUE){
    return(cbind(1:length(regression_coeff),abs(regression_coeff)))
  }
  real_indices=c()
  #Only keep the non-zero indices.
  want_indices_elastic=which(0<abs(regression_coeff))
  #If the radius is not zero, we include the 'n' neighboring features as well.
  for (i in 1:(length(want_indices_elastic))){
    for (j in 0:radius){
      real_indices=c(real_indices,want_indices_elastic[i]+j)
      real_indices=c(real_indices,want_indices_elastic[i]-j)
    }
  }
  real_indices=unique(real_indices)
  return(real_indices)
}


#### Lasso ####
#' Conducts the Lasso Feature Selection
#'
#' @description \code{find_lasso_variables} returns a vector of the selected features/indices of the curve
#'
#' The function uses Lasso to conduct feature selection. A cross validated lasso model is fit on the binary response data, and the feature importances are extracted
#' As lasso sets many of the features to 0, the selected features are those that are non-zero.
#' These features are returned along with a window of neighboring features if the user wishes (typically recommended).
#'
#' @param data (nxm matrix) A nxm matrix containing the covariates and responses.
#' @param radius (positive integer ). An integer parameter specifying if, and how many neighboring features should be considered for feature selection.
#' This is done to capture critical points that may be 'close' to the selected feature.
#' @return The output is a vector of indices/features to be selected.
#'
find_lasso_variables=function(data,radius=0,weights=FALSE){
  #Transforming the data to -1, and 1 if it isn't already for logistic regression purposes.
  data[,1]=ifelse(data[,1]>0,1,0)
  #Initialize lasso model
  regression_model=glmnet::cv.glmnet(data[,-1], data[,1], alpha = 0.95,intercept = FALSE,family='binomial')
  #Extract the coefficients of the variables.
  tmp_coeffs = as.matrix(stats::coef(regression_model,s=regression_model$lambda.1se))
  regression_coeff=as.vector(tmp_coeffs)
  #Don't include the intercept.
  regression_coeff=regression_coeff[2:length(regression_coeff)]
  if (weights==TRUE){
    return(cbind(1:length(regression_coeff),abs(regression_coeff)))
  }
  real_indices=c()
  #Include the non-zero features only to be selected.
  want_indices_lasso=which(0<abs(regression_coeff))
  #If the radius is not zero, we include the 'n' neighboring features as well.
  for (i in 1:(length(want_indices_lasso))){
    for (j in 0:radius){
      real_indices=c(real_indices,want_indices_lasso[i]+j)
      real_indices=c(real_indices,want_indices_lasso[i]-j)
    }
  }
  real_indices=unique(real_indices)
  return(real_indices)
}
