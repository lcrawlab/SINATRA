### Predictor Prioritization via RelATive cEntrality (RATE) centrality measures ###

### NOTE: This function assumes that one has already obtained (posterior) draws/estimates of a nonparametric or nonlinear function as suggested in Crawford et al. (2018) ###

### Review of the function parameters ###
#'X' is the nxp design matrix (e.g. genotypes) where n is the number of samples and p is the number of dimensions. This is the original input data
#'f.draws' is the Bxn matrix of the nonparametric model estimates (i.e. f.hat) with B being the number of sampled (posterior) draws;
#'nullify' is an optional vector specifying a given predictor variable effect that user wishes to "remove" from the data set. An example of this corresponds to Figures 2(b)-(d) in Crawford et al. (2018);
#'snp.nms' is an optional vector specifying the names of the predictor variables;
#'cores' is a parameter detailing the number of cores to parallelize over. If too many are assigned, RATE will set this variable to maximum number of cores that are available on the operating machine.

#NOTE: As decribed in the Supplementary Material of Crawford et al. (2018), we implement different matrix factorizations
#for two cases: (i) n > p, and (ii) n < p. Again, n is the number of samples and p is the number of predictors.

######################################################################################
######################################################################################
######################################################################################

RATEv2 = function(X, f.draws = NULL, pre.specify = FALSE, beta.draws = NULL, prop.var =1, rank.r = min(nrow(X),ncol(X)), nullify = NULL,snp.nms = NULL, cores = 1){
  
  ### Install the necessary libraries ###
  usePackage("doParallel")
  usePackage("MASS")
  usePackage("Matrix")
  usePackage("svd")
  
  ### Determine the number of Cores for Parallelization ###
  if(cores > 1){
    if(cores>detectCores()){warning("The number of cores you're setting is larger than detected cores!");cores = detectCores()}
  }
  
  ### Register those Cores ###
  registerDoParallel(cores=cores)
  
  if(pre.specify == FALSE){
    
    if(is.null(f.draws)){stop("The function draws have not been specified!")}
    
    ### Take the SVD of the Design Matrix for Low Rank Approximation ###
    svd_X = propack.svd(X,rank.r); 
    dx = svd_X$d > 1e-10
    px = cumsum(svd_X$d^2/sum(svd_X$d^2)) < prop.var
    r_X = dx&px 
    u = with(svd_X,(1/d[r_X]*t(u[,r_X])))
    v = svd_X$v[,r_X]
    
    # Now, calculate Sigma_star
    SigmaFhat = cov(f.draws)
    Sigma_star = u %*% SigmaFhat %*% t(u)
    
    # Now, calculate U st Lambda = U %*% t(U)
    svd_Sigma_star = propack.svd(Sigma_star,rank.r)
    r = svd_Sigma_star$d > 1e-10
    U = t(ginv(v)) %*% with(svd_Sigma_star, t(1/sqrt(d[r])*t(u[,r])))
    
    ### Create Lambda ###
    Lambda = tcrossprod(U)
    
    ### Compute the Kullback-Leibler divergence (KLD) for Each Predictor ###
    #mu = c(ginv(X)%*%colMeans(f.draws))
    mu = v%*%u%*% colMeans(f.draws)
    int = 1:length(mu); l = nullify;
    
    if(length(l)>0){int = int[-l]}
    
    KLD = foreach(j = int, .combine='c')%dopar%{
      q = unique(c(j,l))
      m = mu[q]
      U_Lambda_sub = qr.solve(U[-q,],Lambda[-q,q,drop=FALSE])
      kld = crossprod(U_Lambda_sub%*%m)/2
      names(kld) = snp.nms[j]
      kld
    }
  }
  else{
    ### Check to Make sure the Effect Size Draws have been Pre-specified ###
    if(is.null(beta.draws)){stop("The effect size draws have not been specified!")}
    
    ### Specify the Effect Size Information ###
    Sigma_star = cov(beta.draws)
    svd_Sigma_star = propack.svd(Sigma_star,rank.r)
    r = svd_Sigma_star$d > 1e-10
    U = with(svd_Sigma_star, t(1/sqrt(d[r])*t(u[,r])))
    
    ### Create Lambda ###
    Lambda = tcrossprod(U)
    
    ### Compute the Kullback-Leibler divergence (KLD) for Each Predictor ###
    mu = colMeans(beta.draws)
    int = 1:length(mu); l = nullify;
    
    if(length(l)>0){int = int[-l]}
    
    KLD = foreach(j = int, .combine='c')%dopar%{
      q = unique(c(j,l))
      m = mu[q]
      U_Lambda_sub = qr.solve(U[-q,],Lambda[-q,q,drop=FALSE])
      kld = crossprod(U_Lambda_sub%*%m)/2
      names(kld) = snp.nms[j]
      kld
    }
  }
  
  ### Compute the corresponding “RelATive cEntrality” (RATE) measure ###
  RATE = KLD/sum(KLD)
  
  ### Find the entropic deviation from a uniform distribution ###
  Delta = sum(RATE*log((length(mu)-length(nullify))*RATE))
  
  ### Calibrate Delta via the effective sample size (ESS) measures from importance sampling ###
  #(Gruber and West, 2016, 2017)
  ESS = 1/(1+Delta)*100
  
  ### Return a list of the values and results ###
  return(list("KLD"=KLD,"RATE"=RATE,"Delta"=Delta,"ESS"=ESS))
}

### Define the Package ###
usePackage <- function(p){
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}