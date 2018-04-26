

require(lacrosse)

D=30
N=60
K=3
P=5
iterations=1000
Bsparsity=0.5
Gsparsity=0.7
noise.var=0.1
plot.flag=F
debug=T
  obsF=matrix(rnorm(P*N), nrow = P, ncol = N)
  B=matrix(rnorm(K*P)* (runif(K * P) < Bsparsity), nrow = K, ncol = P)
  X=B %*% obsF
  G=matrix(rnorm(D*K)* (runif(D * K) < Gsparsity), nrow = D, ncol = K)
  Y = G %*% X + sqrt(noise.var)*matrix(rnorm(D*N), nrow = D, ncol = N)
  mv = matrix(runif(D * N) < 0.1, nrow = D, ncol = N)
  settings=default.settings()
  settings$verbose=T
  settings$isotropic.noise=T
  settings$hierarchical.noise.prior=T
  settings$sample.noise.variance=T
  settings$sample.mixing.matrix.hypers=T
  settings$sample.alpha=F
  settings$debug=debug
  settings$lasso=F
  settings$sample.Bprior=T
  settings$predictive.performance=1
  settings$iterations=iterations
  settings$truncate=T
  settings$initial_size=5
#  K=0
  init.params=initial.param(N,D,P,settings)

  mrf=matrix(runif(D*D),ncol=D,nrow=D)<0.1
  mrf=mrf | t(mrf)
  diag(mrf)=F
  mrf.connections=apply(mrf,1,which)
  mrf.strengths=lapply(mrf.connections,function(x) as.vector(array(1.0,length(x))))
  settings$mrf.connections=mrf.connections
  settings$mrf.strengths=mrf.strengths
  res=run.nsfa(Y = t(Y), X=t(obsF), settings=settings, missing.values = t(mv), init.param = init.params)
#  res=run.nsfa(Y = t(Y), X=matrix(0,nrow=N,ncol=K), settings=settings, missing.values = t(mv), init.param = init.params)
  
  # plot some diagnostics
  if (plot.flag){
      par(mfrow=c(3,1))
      plot(res$loglikelihood,type='b',xlab="iteration",ylab="log likelihood")
      plot(res$test.loglike,type='b',xlab="iteration",ylab="test log likelihood")
      plot(res$num.features,type='b',xlab="iteration",ylab="# features")
	# plot the best factor loading matrix	
  g=res$best.params$factor.loadings	
  if (require(RColorBrewer)) {	
    col.map=colorRampPalette(c("red","black","blue"))(256)
  } else {
    col.map=topo.colors(256)
  }
  #heatmap(g,xlab="features",ylab="dimensions",Rowv=NA,Colv=NA,col=col.map,scale="none",main="Best factor loading matrix")
}	
