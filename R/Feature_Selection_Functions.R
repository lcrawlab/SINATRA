#Feature Selection Functions
library(Rcpp)
library(kernlab)
library(doParallel)
sourceCpp("BAKRGibbs.cpp")
source('RATE.R')
library(truncnorm)
library(varbvs)
library(matrixcalc)
library(glmnet)
library(svd)
library(mvtnorm)
### Rate ###

find_rate_variables=function(gp_data,radius=0,bandwidth = 0.01){
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
  #Dimension of Res should be 50 samples
  res = RATEv2(X=X,f.draws=fhat.samples,prop.var = 1,snp.nms = colnames(X),cores = cores)

  #What we want
  rates=res$RATE
  #print(res$RATE)
  fhat.samples=0
  B=0
  W=0
  L=0
  f=0
  Kn=0
  want=rates>1/length(rates)
  numeric_want=as.numeric(want)
  want_indices=which(1==numeric_want)
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
find_lasso_variables=function(data,radius=0){
  data[,1]=ifelse(data[,1]>0,1,0)
  regression_model=cv.glmnet(data[,-1], data[,1], alpha = 0.95,intercept = FALSE,family='binomial')
  #Extract the coefficients
  tmp_coeffs = as.matrix(coef(regression_model,s=regression_model$lambda.1se))
  regression_coeff=as.vector(tmp_coeffs)
  regression_coeff=regression_coeff[2:length(regression_coeff)]
  real_indices=c()
  want_indices_lasso=which(0<abs(regression_coeff))
  for (i in 1:(length(want_indices_lasso))){
    for (j in 0:radius){
      real_indices=c(real_indices,want_indices_lasso[i]+j)
      real_indices=c(real_indices,want_indices_lasso[i]-j)
    }
  }
  real_indices=unique(real_indices)
  return(real_indices)
}
#### PPA ####
find_bayesian_variables=function(gp_data,param=0.5,radius=0){
  gp_data[,1]=ifelse(gp_data[,1]>0,1,0)
  fit=varbvs(X = gp_data[,-1],Z = NULL, y= gp_data[,1],family='binomial',verbose = FALSE)
  bayesian_probs=fit$pip
  want_indices_bayesian=which(param<bayesian_probs)
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
find_elastic_variables=function(data,radius=0){
  data[,1]=ifelse(data[,1]>0,1,0)
  regression_model=cv.glmnet(data[,-1], data[,1], alpha = 0.5,intercept = FALSE,family='binomial')
  #Extract the coefficients
  tmp_coeffs = as.matrix(coef(regression_model,s=regression_model$lambda.1se))
  regression_coeff=as.vector(tmp_coeffs)
  regression_coeff=regression_coeff[2:length(regression_coeff)]
  real_indices=c()
  want_indices_elastic=which(0<abs(regression_coeff))
  for (i in 1:(length(want_indices_elastic))){
    for (j in 0:radius){
      real_indices=c(real_indices,want_indices_elastic[i]+j)
      real_indices=c(real_indices,want_indices_elastic[i]-j)
    }
  }
  real_indices=unique(real_indices)
  return(real_indices)
}
#### Other Functions ####
# Gaussian Process Classification - data is assumed to have classes 1,-1 on the first column
gpc = function(gp_data,test_indices){
    ind=test_indices
    gp_data[,1]=ifelse(gp_data[,1]>0,1,-1)
    train_data <- gp_data[-ind,-1]
    train_labels <- gp_data[-ind,1]
    test_data <- gp_data[ind,-1]
    test_labels <- gp_data[ind,1]
    # fit the GP
    gp <- gausspr(x = train_data,y = train_labels, kernel = 'rbfdot', kpar = list(sigma = 0.015))
    label_probabilities <- predict(gp, test_data, type = 'probabilities')
    predicted_labels <- ifelse(label_probabilities > 0, 1,-1)

    # get the results of classification
    #cm = table(predicted_labels, test_labels)
    #correct = (cm[1,1] + cm[2,2])/sum(cm)
    #correct
    correct=sum(predicted_labels==test_labels)/length(test_labels)
    correct
}

# assumption is that k divides the length of data

### Need to implement stratified
kfoldcvgp = function(k = 10,data){
    if ((dim(data)[1]%%k) != 0){
        k=k-1
        return(kfoldcvgp(k=k,data))
    }
    cv.errors = rep(0,k)
    indices = sample(1:dim(data)[1])
    inc = dim(data)[1]/k
    for (i in 1:k){
        cv.errors[i] = gpc(data,indices[((i-1)*inc+1):(i*inc)])
    }
    c(mean(cv.errors),std(cv.errors))

}
#### For accuracy with multiple methods ####
gpc_multiple = function(gp_data,test_indices){
  ind=test_indices
  train_data <- gp_data[-ind,-1]
  train_labels <- gp_data[-ind,1]
  train_labels_gpc=ifelse(train_labels>0,1,-1)
  test_data <- gp_data[ind,-1]
  test_labels <- gp_data[ind,1]
  test_labels_gpc=ifelse(test_labels>0,1,-1)
  # fit the GP
  gp <- gausspr(x = train_data,y = train_labels_gpc, kernel = 'rbfdot', kpar = list(sigma = 0.015))
  lasso=cv.glmnet(train_data, train_labels, alpha = 0.95,intercept = FALSE,family='binomial')
  bayesian=varbvs(X = train_data,Z = NULL, y= train_labels,family='binomial',verbose = FALSE)
  elastic=cv.glmnet(train_data, train_labels, alpha = 0.5,intercept = FALSE,family='binomial')

  label_probabilities_gp <- predict(gp, test_data, type = 'probabilities')
  label_probabilities_lasso <- predict(lasso, test_data)
  label_probabilities_bayesian <- predict(bayesian, test_data)
  label_probabilities_elastic= predict(elastic, test_data)
  predicted_labels_gp <- ifelse(label_probabilities_gp > 0, 1,-1)
  predicted_labels_lasso <- ifelse(label_probabilities_lasso > 0, 1,0)
  predicted_labels_bayesian <- ifelse(label_probabilities_bayesian > 0, 1,0)
  predicted_labels_elastic <- ifelse(label_probabilities_elastic > 0, 1,0)
  # get the results of classification
  #cm = table(predicted_labels, test_labels)
  #correct = (cm[1,1] + cm[2,2])/sum(cm)
  correct_gp=sum(predicted_labels_gp==test_labels_gpc)/length(test_labels)
  correct_lasso=sum(predicted_labels_lasso==test_labels)/length(test_labels)
  correct_bayesian=sum(predicted_labels_bayesian==test_labels)/length(test_labels)
  correct_elastic=sum(predicted_labels_elastic==test_labels)/length(test_labels)
  return(c(correct_gp,correct_lasso,correct_bayesian,correct_elastic))
}

# assumption is that k divides the length of data
kfoldcvgp_multiple = function(iter=100,k = 5,data){
  if ((dim(data)[1]%%k) != 0){
    k=k-1
    return(kfoldcvgp(k=k,data))
  }
  cv.errors = matrix(NA,nrow=iter,ncol=4)
  colnames(cv.errors)=c('GP','Lasso','PPA','Elastic')
  inc = dim(data)[1]/k
  for (i in 1:iter){
    indices = sample(1:dim(data)[1],inc)
    cv.errors[i,] = gpc_multiple(data,indices)
  }
  mean(cv.errors)
  summary=apply(X = cv.errors,MARGIN = 2,FUN = mean)
  return(summary)
  # output std?
}



