# Implement Expectation Propagation for Gaussian Process Classification
library(Matrix)
library(FastGP)

###### Gaussian Approximation to Posterior of Latent Values ######
##############
# Approximates the latent posterior with a Gaussian distributions; it does so by moment matching.
# Pseudocode taken from Rasmussen and Williams, Chapter 3. This function outputs the mean and covariance of the
# approximated posterior. To actually generate samples from the latent posterior, generate samples from a multivariate
# normal with the parameters returned by this function.
#
# Inputs:
#  K : the covariance matrix for the GP model
#  class_labels: +/- 1 values indicating the class labels of the data points

ExpectationPropagation <- function(K, class_labels){
  n = length(class_labels)
  nu_tilde = matrix(0,nrow = n,ncol = 1)
  tau_tilde = matrix(0,nrow = n,ncol = 1)
  sigma = K
  mu = matrix(0,nrow = n,ncol = 1)
  
  # change this stopping condition
  for(j in 1:500){
    for (i in 1:n){
      tau_minus_i = sigma[i,i]^-2 - tau_tilde[i]
      nu_minus_i = (sigma[i,i]^-2)*mu[i] - nu_tilde[i]
      
      mu_minus_i = nu_minus_i/tau_minus_i
      sigma_minus_i = tau_minus_i^-0.5
      
      # Compute Marginal Moments
      z_i = class_labels[i]*mu_minus_i/(sqrt(1+sigma_minus_i^2))
      mu_hat_i = mu_minus_i + (class_labels[i] * sigma_minus_i^2 * dnorm(z_i) )/(pnorm(z_i)*sqrt(1+sigma_minus_i^2))
      sigma_hat_i = sqrt( sigma_minus_i^2 - (sigma_minus_i^4*dnorm(z_i))/((1+sigma_minus_i^2)*pnorm(z_i))*(z_i + dnorm(z_i)/pnorm(z_i)) )
      
      # Update Site Parameters
      delta_tau_tilde = sigma_hat_i^-2 - tau_minus_i - tau_tilde[i]
      tau_tilde[i] = tau_tilde[i] + delta_tau_tilde
      nu_tilde[i] = (sigma_hat_i^-2)*mu_hat_i - nu_minus_i
      
      # Update Sigma, mu - the parameters of the posterior
      sigma = sigma - (( delta_tau_tilde^-1 + sigma[i,i])^-1)*sigma[,i]%*%t(sigma[,i])
      mu = sigma%*%nu_tilde
    }
    #Recompute posterior parameters
    S_tilde = diag(as.vector(tau_tilde))
    L = chol(diag(n) + sqrt(S_tilde)%*%K%*%sqrt(S_tilde))
    V = solve(t(L),sqrt(S_tilde)%*%K)
    sigma = K - t(V)%*%V
    mu = sigma%*%nu_tilde
  }
  params <- list(mu,sigma)
  return(params)
}
##############
# Approximates the latent posterior with a Gaussian distribution; it does so by finding the mode of the posterior, and
# using the Hessian (second order Taylor expansion) as an approximation of the covariance. Newton Raphson is used to find
# the mode of the posterior.
#
# Pseudocode taken from Rasmussen and Williams, Chapter 3. This function outputs the mean and covariance of the
# approximated posterior. To actually generate samples from the latent posterior, generate samples from a multivariate
# normal with the parameters returned by this function.
#
# Inputs:
#  Kn : the covariance matrix for the GP model
#  class_labels: +/- 1 values indicating the class labels of the data points

LaplaceApproximation <- function(Kn, class_labels){
  f <- rep(0,n)
  #Change stopping conditions for convergence
  # do Newton IRLS procedure to use Laplace / Gaussian approximation.
  for(k in 1:1000){
    W <- diag(as.vector(sigmoid(f)*(1-sigmoid(f))))
    B <- diag(x = 1,n) + sqrt(W) %*% Kn %*% sqrt(W)
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

# Based on Iain Murray's paper 'Elliptical Slice Sampling'. Implemented using the FastGP package. The function returns
# the desired number of mcmc samples
# Inputs:
#  K: the covariance matrix of the GP model
#  class_labels: the class labels for each data point, +/- 1.
#  num_mcmc_samples: the number of desired mcmc samples to be returned
#  probit: set TRUE if the link function in the model is probit; otherwise the function uses the logistic link.
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

ess_personal_implementation <- function(n){
  
}

# Based on Radford Neal's paper 
Neal_Metropolis_Hastings <- function(n){
  
}

# Based on Titsias and Lawrence
MCMC_Control_Variables <- function(n){
  
}

##### Helper Functions ######
sigmoid <- function(x){
  return(1/(1+exp(-x)))
}

# When the link function is probit
probit_log_likelihood <- function(latent_variables, class_labels){
  return(sum(log(pnorm(latent_variables*class_labels))))
}

# When the link function is logistic 
logistic_log_likelihood <- function(latent_variables, class_labels){
  return(-sum(log(1+exp(latent_variables*class_labels))))
}


