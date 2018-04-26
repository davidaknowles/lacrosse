

require(dnsfa)

test.dnsfa <- function(D=30,N=100,K=5,P=300,iterations=5000,Bsparsity=0.1,Gsparsity=0.3,noise.var=0.03,make.plot=F,lasso=F){
  cat("tjhere2")
  obsF=matrix(rnorm(P*N), nrow = P, ncol = N)
  B=matrix(rnorm(K*P)* (runif(K * P) < Bsparsity), nrow = K, ncol = P)
  X=B %*% obsF
  G=matrix(rnorm(D*K)* (runif(D * K) < Gsparsity), nrow = D, ncol = K)
  Y = G %*% X + sqrt(noise.var)*matrix(rnorm(D*N), nrow = D, ncol = N)
  mv = matrix(runif(D * N) < 0.1, nrow = D, ncol = N)
  settings=default.settings()
  settings$verbose=T
  settings$isotropic.noise=F
  settings$hierarchical.noise.prior=T
  settings$sample.noise.variance=T
  settings$sample.mixing.matrix.hypers=T
  settings$sample.alpha=T
  settings$debug=F
  settings$lasso=lasso
  settings$sample.Bprior=T
  settings$predictive.performance=1
  settings$iterations=iterations
  settings$truncate=T
#  K=0
  init.params=initial.param(N,D,P,settings)
  cat("enter")
  res=run.nsfa(Y = t(Y), X=t(obsF), settings=settings, missing.values = t(mv), init.param = init.params)
  cat("exit")
  # plot some diagnostics
  if (make.plot){
      par(mfrow=c(3,1))
      plot(res$loglikelihood,type='b',xlab="iteration",ylab="log likelihood")
      plot(res$test.loglike,type='b',xlab="iteration",ylab="test log likelihood")
      plot(res$num.features,type='b',xlab="iteration",ylab="# features")
      # plot the best factor loading matrix			      
      g=res$best.params$factor.loadings
      col.map=if (require(RColorBrewer)) colorRampPalette(c("red","black","blue"))(256) else topo.colors(256)
      x11()
      heatmap(as.matrix(g),xlab="features",ylab="dimensions",Rowv=NA,Colv=NA,col=col.map,scale="none",main="Best factor loading matrix")
  }
  list(res=res,true.B=B,true.G=G)
}

    args=if (interactive()) c(1,1) else commandArgs(trailingOnly=T)
    cat(args[1],"\n")
    cat(args[2],"\n")
    lasso=as.logical(as.integer(args[1]))
    cat("here")
    seed=as.integer(args[2])
    cat("there")
    set.seed(seed)
    cat("Test DNFSA seed ",seed," lasso? ",lasso,"\n")
    res=test.dnsfa(lasso=lasso)
    save(res,file=paste("~/scratch/dnsfa_tests/test",seed,if (lasso) "lasso" else "",".RData"))
