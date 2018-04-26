package.name="lacrosse"

default.settings = function(){
  list(debug=F,
    a=1.0, # noise variance hyperparameters
    b=1.0,
    a0=1.0, # hyperparameters on b
    b0=1.0,
    d=1.0, # factor variance hyperparameters
    c0=1.0,
    d0=1.0,
    e=1.0, # IBP alpha hyperparamters
    f=1.0,
    alpha=1.0, # IBP hyperparameters
    beta=1.0,
    iterations=1000,
    verbose=F,
    lasso=F,
    store.samples=F, 
    sample.noise.variance=T, 
    lambda.e=10.0, # noise precision
    sample.mixing.matrix.hypers=T, 
    lambda.x=1.0, # factor precision
    sample.beta=F,
    sample.Bprior=F,
    # Interesting settings
    sample.alpha=F, # whether to sample the IBP concentration parameter
    isotropic.noise=T, # enforce isotropic noise?
    perfactor.variance=T, # Per factor covariances? With c=1 corresponds to ARD prior
    learn.scale=T, # Hierarchical prior on the factor covariances? Only valid if perfactor is T
    c=1.0, # 1 corresponds to an ARD prior. 2 corresponds to not. 
    hierarchical.noise.prior=T, # share power across noise dimensions? 
    predictive.performance=F,
    truncate=F,
    initial_size=10,
    mrf.connections=list(),
    mrf.strengths=list(),
    feature.mrf.connections=list(),
    feature.mrf.strengths=list(),
    log.file="",
    samples.dir="",
    thin.out=5,
    burn.in=500) 
}

initial.param <- function(N,D,P,settings){
  .Call("initial_param", N, D, P, settings, PACKAGE=package.name)
}

test.bprior <- function(N=1000,its=100){
    .Call("test_Bprior", N, its, PACKAGE=package.name)
}

add.names <- function(params,Y){
  rownames(params$factor.loadings)=rownames(Y)
  colnames(params$factors)=colnames(Y)
  names(params$lambda.e)=rownames(Y)
  params
}

check.consistency <- function(Y, missing.values, init.param, mrf.connections, mrf.strengths)
{
  stopifnot(nrow(Y)==nrow(missing.values) && ncol(Y)==ncol(missing.values))
  stopifnot(length(mrf.connections)==nrow(Y) || length(mrf.connections)==0)
  stopifnot(length(mrf.strengths)==nrow(Y) || length(mrf.strengths)==0)
}

run.lacrosse <- function(Y,X,settings,missing.values,init.param){
  Y=as.matrix(t(Y)) # convert if data.frame, and use R standard: cols as variables
  X=as.matrix(t(X))
  if (missing(settings))
    settings=default.settings()
  if (settings$store.samples)
    dir.create(settings$samples.dir)
  if (missing(init.param))
    init.param=initial.param(ncol(Y),nrow(Y),nrow(X),settings)
  # may need to convert from sparse to dense B/G
  init.param$B=as.matrix(init.param$B)
  init.param$factor.loadings=as.matrix(init.param$factor.loadings)
  if (missing(missing.values)){
    if (any(is.na(Y))|any(is.nan(Y))){
      missing.values=is.na(Y)|is.nan(Y)
      settings$predictive.performance=F
      Y[missing.values]=0.0
      missing.values=matrix(as.double(missing.values),nrow=nrow(missing.values),ncol=ncol(missing.values))
      check.consistency(Y,missing.values,init.param,settings$mrf.connections,settings$mrf.strengths)
      res=.Call( "run_lacrosse_missing", Y, X, missing.values, settings, init.param, PACKAGE = package.name )
      res$imputed.Y=t(res$imputed.Y)
    }
    else {
      settings$predictive.performance=F
      check.consistency(Y,missing.values,init.param,settings$mrf.connections,settomgs$mrf.strengths)
      res=.Call("run_lacrosse", Y, X, settings, init.param, PACKAGE = package.name )
    }
  }
  else {
    missing.values=t(matrix(as.double(missing.values),nrow=nrow(missing.values),ncol=ncol(missing.values)))
    settings$predictive.performance = !(any(is.na(Y))|any(is.nan(Y)));
    check.consistency(Y,missing.values,init.param,settings$mrf.connections,settings$mrf.strengths)
    res=.Call("run_lacrosse_missing", Y, X, missing.values, settings, init.param, PACKAGE = package.name )
    res$imputed.Y=t(res$imputed.Y)
  }
  res$best.params=add.names(res$best.params,Y)
  res$final.params=add.names(res$best.params,Y)
  res
}

test.lacrosse <- function(D=30,N=60,K=3){
  G=matrix(rnorm(D*K)* (runif(D * K) < 0.7), nrow = D, ncol = K)
  X=matrix(rnorm(K*N), nrow = K, ncol = N)
  Y = G %*% X + 0.1*matrix(rnorm(D*N), nrow = D, ncol = N)
  mv = matrix(runif(D * N) < 0.1, nrow = D, ncol = N)
  settings=default.settings()
  settings$verbose=T
  #settings$isotropic.noise=F
  settings$hierarchical.noise.prior=T
  settings$sample.noise.variance=T
  settings$sample.mixing.matrix.hypers=F
  settings$sample.alpha=F
  settings$debug=0
  settings$predictive.performance=1
  settings$iterations=100
  init.params=initial.param(ncol(Y),nrow(Y),K,settings)

  mrf=matrix(runif(D*D),ncol=D,nrow=D)<0.1
  mrf=mrf | t(mrf)
  diag(mrf)=F
  mrf.conn=apply(mrf,1,which)
  mrf.strengths=lapply(mrf.conn,function(x) as.vector(array(1.0,length(x))))
  settings$mrf.connections=mrf.connections
  settings$mrf.strengths=mrf.strengths
  res=run.lacrosse(Y = t(Y), X=t(X), settings=settings, missing.values = t(mv), init.param = init.params)
  
  # plot some diagnostics
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
  res
}
