
require(doMC)
registerDoMC(cores=10)

source("kbmtl.R")
source("load_CTRPv2.R")

sens=y
x=scale(x)
x[is.na(x) | is.nan(x)]=0.0

results=foreach(cv.fold=1:10) %dopar% {
  
  set.seed(cv.fold)
  n.cell.lines=nrow(sens)
  perm=sample(n.cell.lines)[1:82]
  test=logical(n.cell.lines)
  test[perm]=T
  sens[perm,]=NA
  
  x.train=x[!test,]
  y.train=scale(sens[!test,])
  #y.train[is.na(y.train)]=0. # kbtml can handle NAs
  
  d=dist(x)
  defaultLengthScale=mean(d)
  
  gram=exp( -.5 * as.matrix(d^2)/defaultLengthScale^2 )
  
  #initalize the parameters of the algorithm
  parameters <- list()
  
  #set the hyperparameters of gamma prior used for projection matrix
  parameters$alpha_lambda <- 1
  parameters$beta_lambda <- 1
  
  #set the hyperparameters of gamma prior used for output noise
  parameters$alpha_epsilon <- 1
  parameters$beta_epsilon <- 1
  
  ### IMPORTANT ###
  #For gamma priors, you can experiment with three different (alpha, beta) values
  #(1, 1) => default priors
  #(1e-10, 1e+10) => good for obtaining sparsity
  #(1e-10, 1e-10) => good for small sample size problems
  
  #set the number of iterations
  parameters$iteration <- 200
  
  #set the subspace dimensionality
  parameters$R <- 20
  
  #set the seed for random number generator used to initalize random variables
  parameters$seed <- 1606
  
  #set the standard deviation of hidden representations
  parameters$sigma_h <- 0.1
  
  #set the standard deviation of weight parameters
  parameters$sigma_w <- 1.0
  
  #initialize the kernel and target outputs for training
  Ktrain <- gram[!test,!test] #should be an Ntra x Ntra matrix containing similarity values between training samples
  Ytrain <- y.train #should be an Ntra x T matrix containing target outputs (contains only real values and NaNs)
  
  #perform training
  state <- kbmtl_regression_train(Ktrain, Ytrain, parameters)
  
  #initialize the kernel for testing
  Ktest <- gram[!test,test] #should be an Ntra x Ntest matrix containing similarity values between training and test samples
  
  #perform prediction
  prediction <- kbmtl_regression_test(Ktest, state)
  
  #display the predicted probabilities
  print(prediction$Y$mu)
  
  #y.test=sweep( sweep( sens[test,], 2, attr(y.train, "scaled:center") ), 2, attr(y.train, "scaled:scale"), FUN ="/" )
  prediction
}

saveRDS(results, "kbtml_CTRPv2.rds")
