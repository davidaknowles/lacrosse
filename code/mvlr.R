# Bayesian Multi-view Multi-task Linear Regression
# Copyright (R), all rights reserved by authors.
#
# Authors: 
# Suleiman Khan (suleiman.khan@helsinki.fi)
# Muhammad Ammad (muhammad.ammad-ud-din@helsinki.fi )

require("rstan")

init_prior_model_nreg <- function(data,beta = NULL, lambda = NULL, W = NULL, tau = NULL,sigma = NULL)
{
  if(length(tau) == 0) { 
    if(length(beta) == 0) { tau <- abs(rcauchy(data$M,location=0,scale=1)); tau = tau/sum(tau); 
    } else { tau <- abs(runif(data$M,0.5,1)); tau = tau/sum(tau); tau = array(tau,dim=length(tau));  }
  }
  if(length(lambda) == 0) { 
    if(length(beta) == 0) { lambda <- abs(rcauchy(data$dX,location=0,scale=1)); inds = lambda>M; if(length(inds)>0) lambda[inds] = M-1e-3; 
    } else { lambda <- 1 + beta; inds = lambda>M; if(length(inds)>0) lambda[inds] = M-1e-3; }
  }
  if(length(sigma) == 0) { sigma <- runif(data$dY,0.5,1); sigma = array(sigma,dim=length(sigma)); } #abs(rcauchy(data$dY,location=0,scale=1));
  if(length(beta) == 0) beta <- abs(rnorm(data$dX,0,(data$lambda_v*lambda)*(tau*data$tau_v) ));
  if(length(W) == 0) { W <- abs(rnorm(data$dY,0,1)); W = as.vector(W/sum(W)); W = array(W,dim=length(W)); }
  list(beta = beta, lambda = lambda, tau = tau,sigma = sigma,W=W)
}

# data$Y : [samples] x [tasks]
# data$X : list of [samples] x [features]
runMVSTAN <- function(data, tol_rel_obj = 1e-3,iter = 50000, eval_elbo=1000)
{

  dX = ncol(data$X)
  bp = rep(1e-3,dX)*sample(c(1,-1),dX,replace=T)
  cYX = cor(data$Y,data$X)
  inds = order(apply(abs(cYX),2,max),decreasing=T)[1:round(dX*0.25,2)]
  for(i in 1:length(inds))
    bp[inds[i]] = cYX[which.max(cYX[,inds[i]]),inds[i]]
  init.dat = init_prior_model_nreg(data,beta=bp)

  print("Compiling Stan...")
  model = rstan::stan_model('nreg.stan');
  
  print("Running Stan...")
  out=1
  counter=1
  class(out)="try-error"
  while(class(out)=="try-error") {
    out <- try(vb(model, 
                  data = data, 
                  algorithm = "meanfield",
                  seed=counter,
                  init = init.dat,
                  tol_rel_obj=tol_rel_obj,
                  iter=iter,
                  eval_elbo=eval_elbo  ))
    counter=counter+1
    if (counter > 5) return(out)
  }
  out
}

predictMVSTAN <- function(out, xTest)
{

  post = list()
  post$W = rstan::extract(out,"W",permuted=TRUE)
  post$beta = rstan::extract(out,"beta",permuted=TRUE)
  ynew = array(NA,dim=c(nrow(xTest),ncol(post$W[[1]]),nrow(post$W[[1]])))
  for(i in 1:nrow(xTest))
  {
    ynew[i,,] = t(post$W[[1]] * matrix(xTest[i,] %*% t(post$beta[[1]]),nrow(post$W[[1]]),ncol(post$W[[1]])))   # analogous to for loop computation - verified empirically  
  }

  return(apply(ynew,c(1,2),mean))
}

getBeta.MV.STAN <- function(out)
{
  beta = rstan::extract(out,"beta",permuted=TRUE)
  beta = t(apply(beta[[1]],c(2),mean))
  return(beta)
}

getPosteriorEV.MV.STAN.CV <- function(run.stan)
{
  for(i in 1:length(run.stan$res))
  {
    tmp = getPosteriorEV.MV.STAN(run.stan$res[[i]]$out)
    if(i == 1) { post = tmp; next; }
    for(j in 1:length(post))
      post[[j]] = post[[j]] + tmp[[j]]
  }
  for(j in 1:length(post))
    post[[j]] = post[[j]]/length(run.stan$res)
  return(post)
}

getPosteriorEV.MV.STAN <- function(out)
{
  post = list()
  post$sigma = rstan::extract(out,"sigma",permuted=TRUE)
  post$sigma = apply(post$sigma[[1]],c(2),mean)
  #estim.noise = post$sigma/apply(Y,2,sd)
  
  post$W = rstan::extract(out,"W",permuted=TRUE)
  post$W = apply(post$W[[1]],c(2),mean)
  post$W = matrix(post$W,1,length(post$W))
  
  post$beta = rstan::extract(out,"beta",permuted=TRUE)
  post$beta = t(apply(post$beta[[1]],c(2),mean))
  
  sc = sum(post$W);
  post$W = post$W * 1./sc
  post$beta = post$beta * sc
  
  post$lambda = rstan::extract(out,"lambda",permuted=TRUE)
  post$lambda = t(apply(post$lambda[[1]],c(2),mean))
  post$tau = rstan::extract(out,"tau",permuted=TRUE)
  post$tau = t(apply(post$tau[[1]],c(2),mean))
  for(iM in 1:M) {
    if(iM>1) {
      sc = sd(post$lambda[,1:dXm[iM]+ sum(dXm[1:(iM-1)])])
      post$lambda[,1:dXm[iM]+ sum(dXm[1:(iM-1)])] = post$lambda[,1:dXm[iM]+ sum(dXm[1:(iM-1)])] / sc
      post$tau[,iM] = post$tau[,iM] * sc
    } else {
      sc = sd(post$lambda[,1:dXm[1]])
      post$lambda[,1:dXm[iM]] = post$lambda[,1:dXm[iM]] / sc
      post$tau[,iM] = post$tau[,iM] * sc
    }
  }
  return(post)
}