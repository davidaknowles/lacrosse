
source("load_CTRPv2.R")
sens=y
x=scale(x)
x[is.na(x) | is.nan(x)]=0.0

isotropic=F
n.drugs=ncol(sens)

require(glmnet)

n.reps=10
methods=4
results=data.frame(matrix(NA,ncol=methods,nrow=n.reps))
train.rmse=data.frame(matrix(NA,ncol=methods,nrow=n.reps))
rmse=function(err) sqrt(mean(err*err))

lacrosse_dir=paste0(scratch_dir,"/lacrosse/") 

for (cv.fold in 1:n.reps){
    meth=1

    print(cv.fold)

    set.seed(cv.fold)
    n.cell.lines=nrow(sens)
    perm=sample(n.cell.lines)[1:82]
    test=logical(n.cell.lines)
    test[perm]=T
    x.train=x[!test,]
    y.train=scale(sens[!test,])
    x.test=x[test,]
    y.mean=attr(y.train,"scaled:center")
    y.scale=attr(y.train,"scaled:scale")

    adjust.y=function(to.adjust) sweep(sweep(to.adjust,2,y.mean,"-"),2,y.scale,"/")
    y.adjusted=adjust.y(sens)
    train.calc=function(y.pred) rmse((y.pred-y.adjusted[!test,])[!is.na(y.adjusted[!test,])])
    test.calc=function(y.pred) rmse((y.pred-y.adjusted[test,])[!is.na(y.adjusted[test,])])
    alpha=1
    use.mrf=T
    truncate=T
    suffix=paste("_dnsfa",cv.fold,"alpha",alpha,"mrf",use.mrf,"truncate",truncate,sep="")

    base.dir=paste(lacrosse_dir,"samples",suffix,"/",sep="")
    
    load(paste("results/results",suffix,".RData",sep=""))
    p=r$best.params

    y_preds=foreach (s=dir(base.dir)) %do% {
      load(paste(base.dir,s,sep=""))
      if (ncol(param$B)>0) { as.matrix(t(param$factor.loadings %*% param$B %*% t(x.test))) } else {NULL}
    }
    y.pred=Reduce("+",y_preds) / length(y_preds)
    
    per.sample=foreach(y_pred_here=y_preds, .combine = c) %do% { rmse((y_pred_here-sens[test,])[!is.na(sens[test,])]) }
    
    results[cv.fold,meth]=test.calc(y.pred)
    
    y.pred=p$factor.loadings %*% p$B %*% t(x.train)
    y.pred=as.matrix(t(y.pred))
    train.rmse[cv.fold,meth]=train.calc(y.pred)
    meth.name="dfa"

    colnames(results)[meth]=meth.name
    meth=meth+1
    
    # glmnet per drug
    print("glmnet")
    load(paste0(lacrosse_dir,"glmnet_single_",cv.fold,".RData"))
    
    co_matrix=lapply(res, as.numeric) %>% as.data.frame() %>% as.matrix() %>% set_colnames(NULL)
    
    y.pred=cbind(1, x.test) %*% co_matrix
    y.pred.train=cbind(1, x.train) %*% co_matrix
    
    results[cv.fold,meth]=test.calc(y.pred)
    train.rmse[cv.fold,meth]=train.calc(y.pred.train)
    
    colnames(results)[meth]="L1"
    meth=meth+1
    
    print("glmnet mgaussian")
    load(paste0(lacrosse_dir,"glmnet_mgaussian/",cv.fold,".RData"))
    res=. # well that's unfortunate
    co_matrix=lapply(res, as.numeric) %>% as.data.frame() %>% as.matrix() %>% set_colnames(NULL)
    
    y.pred=cbind(1, x.test) %*% co_matrix
    y.pred.train=cbind(1, x.train) %*% co_matrix
    
    results[cv.fold,meth]=test.calc(y.pred)
    train.rmse[cv.fold,meth]=train.calc(y.pred.train)
    
    colnames(results)[meth]="mgaussian"
    meth=meth+1

    results[cv.fold,meth]=test.calc(y.adjusted[test,]*0)
    train.rmse[cv.fold,meth]=train.calc(y.adjusted[!test,]*0)
    colnames(results)[meth]="mean"
    meth=meth+1
    
}

colnames(train.rmse)=colnames(results)

#save(train.rmse,results,file="results.RData")
