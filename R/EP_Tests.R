# Test EP
library(KRLS) # for gaussian kernel
library(mvtnorm) # for multivariate normal
source('RATE.R')

# Toy example - the class will be determined by the last three predictors, if their sum is postive or negative.
# All other variables are not relevant

#### Generate Data
num_data_points = 200
data <- matrix(rnorm(20*num_data_points,0,1),nrow = num_data_points,ncol = 20)
class_labels <- ifelse(rowSums(data[,18:20]) > 0,-1,1)
data <- cbind(class_labels,data)


#### Set up the GPC model
n <- dim(data)[1]
X <- data[,-1]
h <- 0.1 # choose the bandwidth
K <- gausskernel(X,1/(2*h^2))
diag(K) <- 1

#### Sample from the latent posterior
params <- ExpectationPropagation(K,class_labels = data[,1])
mu <- params[[1]]
sigma <- params[[2]]

posterior_samples <- rmvnorm(1e4,mu, sigma)

# Uncomment the line below if MCMC sampler is desired
#posterior_samples <- Elliptical_Slice_Sampling(K,class_labels,1e4)

#### Run RATE
cores = cores=detectCores()
res = RATEv2(X=X,f.draws=posterior_samples,prop.var = 1,snp.nms = colnames(X),cores = cores)
rates=res$RATE
