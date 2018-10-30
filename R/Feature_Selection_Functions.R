
#Feature Selection Functions
library(Rcpp)
library(kernlab)
library(doParallel)
sourceCpp("sinatra/R/BAKRGibbs.cpp")
source('sinatra/R/RATE.R')
library(truncnorm)
library(varbvs)
library(matrixcalc)
library(glmnet)
library(svd)
library(mvtnorm)
### Rate ###

#' Conducts the Rate Feature Selection
#'
#' \code{find_rate_variables} returns a vector of the selected features/indices of the curve
#'
#' The function uses RATE (Crawford 2018) to conduct feature selection. The covariance matrix is estimated via Lapace approximation.
#' Samples are then drawn from the posterior with the approximated covariance matrix, and fed into the RATE function. See Crawford et.al for RATE details.
#' Features are deemed significant if they have weight or 'RATE' 1/p, where p is the number of total variables. These features are returned along with a window
#' of neighboring features if the user wishes (typically recommended).
#'
#' Another use of this function is to return the weights for each feature, for informative purposes. This is given in the parameter 'weights'.
#'
#' @param gp_data (nxm matrix) A nxm matrix containing the covariates and responses.
#' @param radius (positive integer ). An integer parameter specifying if, and how many neighboring features should be considered for feature selection.
#' This is done to capture critical points that may be 'close' to the selected feature.
#' @param bandwidth: (float (0-\infinity))
#' @param weights: (TRUE/FALSE) Returns the RATE weights for each feature if TRUE, otherwise it returns the 'selected' weights.
#' @return The output is a vector of indices/features to be selected.
#'
find_rate_variables=function(gp_data,radius=0,bandwidth = 0.01,weights=FALSE){
  n <- dim(gp_data)[1]
  X <- gp_data[,-1]
  gp_data[,1]=ifelse(gp_data[,1]>0,1,-1)
  h <- bandwidth #median(dist(X)) #need a better choice of this; how does bandwidth affect kernel choice?
  #RATE
  f <- rep(0,n)
  Kn <- GaussKernel(t(X),1/(2*h^2))
  diag(Kn)=1
  #Change stopping conditions for convergence
  # do Newton IRLS procedure to use Laplace / Gaussian approximation.
  for(k in 1:600){
    W <- diag(as.vector(sigmoid(f)*(1-sigmoid(f))))
    B <- diag(x = 1,n) + sqrt(W) %*% Kn %*% sqrt(W)
    #Kinda show
    L <- chol(B)
    b <- W%*%f + (gp_data[,1]+1)/2 - sigmoid(f)
    #Kinda slow
    #L has components that
    #Try to not use all the solves
    a <- b - solve(sqrt(W)%*%t(L),solve(L,sqrt(W)%*%Kn%*%b))
    f <- Kn%*%a
  }
  v = solve(L, sqrt(W)%*%Kn)
  # generate samples from approximate posterior
  #fhat samples may not be right dimensions
  fhat.samples = rmvnorm(1e4,f, Kn - t(v)%*% v)
  #is.positive.semi.definite(q_hat)
  # use RATE:
  cores = cores=detectCores()
  res = RATEv2(X=X,f.draws=fhat.samples,prop.var = 1,snp.nms = colnames(X),cores = cores)

  #The weights of RATE
  rates=res$RATE
  #Clear the elements from memory
  fhat.samples=0
  B=0
  W=0
  L=0
  f=0
  Kn=0
  if (weights==TRUE){
    return(cbind(1:length(rates),rates))
  }
  #Pick out the features/indices with weight 1/p
  want=rates>1/length(rates)
  numeric_want=as.numeric(want)
  want_indices=which(1==numeric_want)
  #If the radius is not zero, we include the 'n' neighboring features as well.
  real_indices=c()
  for (i in 1:(length(want_indices))){
    for (j in 0:radius){
      real_indices=c(real_indices,want_indices[i]+j)
      real_indices=c(real_indices,want_indices[i]-j)
    }
  }
  real_indices=unique(real_indices)
  return(real_indices)
}

#### Lasso ####
#' Conducts the Lasso Feature Selection
#'
#' \code{find_lasso_variables} returns a vector of the selected features/indices of the curve
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
find_lasso_variables=function(data,radius=0){
  #Transforming the data to -1, and 1 if it isn't already for logistic regression purposes.
  data[,1]=ifelse(data[,1]>0,1,0)
  #Initialize lasso model
  regression_model=cv.glmnet(data[,-1], data[,1], alpha = 0.95,intercept = FALSE,family='binomial')
  #Extract the coefficients of the variables.
  tmp_coeffs = as.matrix(coef(regression_model,s=regression_model$lambda.1se))
  regression_coeff=as.vector(tmp_coeffs)
  #Don't include the intercept.
  regression_coeff=regression_coeff[2:length(regression_coeff)]
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
#### Bayesian Variable Selection ####

#' Conducts the Bayesian Feature Selection
#'
#' \code{find_bayesian_variables} returns a vector of the selected features/indices of the curve
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
find_bayesian_variables=function(data,param=0.5,radius=0){
  #Convert the response to 0 or 1, for binary classification in the varbvs model
  data[,1]=ifelse(data[,1]>0,1,0)
  #Fit the model
  fit=varbvs(X = data[,-1],Z = NULL, y= data[,1],family='binomial',verbose = FALSE)
  #Find the pips for each feature.
  bayesian_probs=fit$pip
  #Only keep the features that are above the threshold (pip).
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
#' \code{find_elastic_variables} returns a vector of the selected features/indices of the curve
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
find_elastic_variables=function(data,radius=0){
  #Transforming the data to -1, and 1 if it isn't already for logistic regression purposes.
  data[,1]=ifelse(data[,1]>0,1,0)
  #Fit the model
  regression_model=cv.glmnet(data[,-1], data[,1], alpha = 0.5,intercept = FALSE,family='binomial')
  #Extract the coefficients
  tmp_coeffs = as.matrix(coef(regression_model,s=regression_model$lambda.1se))
  regression_coeff=as.vector(tmp_coeffs)
  #Don't include the intercept
  regression_coeff=regression_coeff[2:length(regression_coeff)]
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


#### Classification/Accuracy Functions ####
#' Conducts Gaussian Process Classification
#'
#' \code{gpc} returns a float representing test accuracy.
#'
#' The function uses a Gaussian Process for classification; the response variables are assumed to be -1,1.
#' An rbf kernel is used, with fixed sigma paramter of 0.015. The Gaussian Process regression is from the package 'kernlab'.
#' Data is split into train and test splits; the indices provided indicate the test indices. All the other non-test indices are used to indicate training data.
#' This function is typically used within a cross validated wrapper function.
#' @param gp_data (nxm matrix) A nxm matrix containing the covariates and responses.
#' @param test_indices (vector of integers). A vector of integers specifying which observations are to be used as test data.
#' @return A float (0-1) of accuracies.
#'
# Gaussian Process Classification - data is assumed to have classes 1,-1 on the first column
gpc = function(gp_data,test_indices){
  ind=test_indices
  #Set the response to -1,1 in case it isn't.
  gp_data[,1]=ifelse(gp_data[,1]>0,1,-1)
  #Split the data
  train_data <- gp_data[-ind,-1]
  train_labels <- gp_data[-ind,1]
  test_data <- gp_data[ind,-1]
  test_labels <- gp_data[ind,1]
  # fit the GP
  gp <- gausspr(x = train_data,y = train_labels, kernel = 'rbfdot', kpar = list(sigma = 0.015))
  label_probabilities <- predict(gp, test_data, type = 'probabilities')
  predicted_labels <- ifelse(label_probabilities > 0, 1,-1)
  #Assess Accuracy
  correct=sum(predicted_labels==test_labels)/length(test_labels)
  correct
}

#'
#' \code{kfoldcvgp} Cross Validated Gaussian Process Classification; returns the cross validated accuracy and standard deviation
#'
#' The function performs a k-fold cross validated Gaussian Process classification for the data.
#' @param data (nxm matrix) A nxm matrix containing the covariates and responses.
#' @param k (integer). How many folds to cross validate on.
#' @return A vector of the average accuracy, and the standard deviation of the errors.
#'
# Gaussian Process Classification - data is assumed to have classes 1,-1 on the first column
kfoldcvgp = function(k = 10,data){
  #Make sure k divides the length of the data.
  if ((dim(data)[1]%%k) != 0){
    k=k-1
    return(kfoldcvgp(k=k,data))
  }
  cv.errors = rep(0,k)
  indices = sample(1:dim(data)[1])
  inc = dim(data)[1]/k
  #Initialize k fold cross validated accuracy.
  for (i in 1:k){
    cv.errors[i] = gpc(data,indices[((i-1)*inc+1):(i*inc)])
  }
  c(mean(cv.errors),std(cv.errors))

}
#### For accuracy with multiple methods ####

#' Conducts Classification for multiple methods (Elastic net, Lasso, Bayesian, Gaussian Process).
#'
#' \code{test_accuracy} returns a vector of floats representing accuracy for each different method.
#'
#' The function uses a Lasso, Elastic Net, Gaussian Processes, and varbvs for classification; the response variables are 0,1 for all methods except Gaussian Process.
#' For Gaussian Processes an rbf kernel is used, with fixed sigma paramter of 0.015. The Gaussian Process regression is from the package 'kernlab'.
#' For lasso and elastic net, a cross validated model is fit with the package 'glmnet'.
#' Bayesian method is fit using varbvs.
#' Data is split into train and test splits; the indices provided indicate the test indices. All the other non-test indices are used to indicate training data.
#' This function is typically used within a cross validated wrapper function.
#' @param data (nxm matrix) A nxm matrix containing the covariates and responses.
#' @param test_indices (vector of integers). A vector of integers specifying which observations are to be used as test data.
#' @return A vector of floats of accuracies
#'
test_accuracy = function(data,test_indices){
  ind=test_indices
  #Training data
  train_data <- data[-ind,-1]
  train_labels <- data[-ind,1]
  train_labels_gpc=ifelse(train_labels>0,1,-1)
  #Test data
  test_data <- data[ind,-1]
  test_labels <- data[ind,1]
  test_labels_gpc=ifelse(test_labels>0,1,-1)
  # fit the GP
  gp <- gausspr(x = train_data,y = train_labels_gpc, kernel = 'rbfdot', kpar = list(sigma = 0.015))
  #Fit lasso
  lasso=cv.glmnet(train_data, train_labels, alpha = 0.95,intercept = FALSE,family='binomial')
  #Fit varbvs
  bayesian=varbvs(X = train_data,Z = NULL, y= train_labels,family='binomial',verbose = FALSE)
  #Fit Elastic Net
  elastic=cv.glmnet(train_data, train_labels, alpha = 0.5,intercept = FALSE,family='binomial')

  #Predict and assess Accuracies
  label_probabilities_gp <- predict(gp, test_data, type = 'probabilities')
  label_probabilities_lasso <- predict(lasso, test_data)
  label_probabilities_bayesian <- predict(bayesian, test_data)
  label_probabilities_elastic= predict(elastic, test_data)
  predicted_labels_gp <- ifelse(label_probabilities_gp > 0, 1,-1)
  predicted_labels_lasso <- ifelse(label_probabilities_lasso > 0, 1,0)
  predicted_labels_bayesian <- ifelse(label_probabilities_bayesian > 0, 1,0)
  predicted_labels_elastic <- ifelse(label_probabilities_elastic > 0, 1,0)

  correct_gp=sum(predicted_labels_gp==test_labels_gpc)/length(test_labels)
  correct_lasso=sum(predicted_labels_lasso==test_labels)/length(test_labels)
  correct_bayesian=sum(predicted_labels_bayesian==test_labels)/length(test_labels)
  correct_elastic=sum(predicted_labels_elastic==test_labels)/length(test_labels)
  #Return everything.
  return(c(correct_gp,correct_lasso,correct_bayesian,correct_elastic))
}

#' Performs
#'
#' \code{kfoldcv_iteration} Elastic Net, Lasso, Bayesian, and Gaussian Process Classification over user-specified iterations with specified train, test splits.
#' Returns a vector of accuracies
#'
#' The function performs a random train test split with the  aforementioned methods over n interations.
#' @param iter (integer (1-\infinity)) How many iterations of train test validation to perform.
#' @param data (nxm matrix) A nxm matrix containing the covariates and responses.
#' @param k (integer). How to determine the train test splits.
#' @return A vector of the average accuracies across the methods with labels for the methods.
#'
kfoldcv_iteration = function(iter=100,k = 5,data){
  #Make sure 'k', the number of folds divide the length of data
  if ((dim(data)[1]%%k) != 0){
    k=k-1
    return(kfoldcvgp(k=k,data))
  }
  cv.errors = matrix(NA,nrow=iter,ncol=4)
  colnames(cv.errors)=c('GP','Lasso','PPA','Elastic')
  inc = dim(data)[1]/k
  #Initialize the accuracy evaluation; over each iteration, randomly sample a portion of the data to use as test data.
  for (i in 1:iter){
    indices = sample(1:dim(data)[1],inc)
    cv.errors[i,] = test_accuracy(data,indices)
  }
  mean(cv.errors)
  summary=apply(X = cv.errors,MARGIN = 2,FUN = mean)
  #Return accuracies.
  return(summary)
}
