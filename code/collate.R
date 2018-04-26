
load("all_gdsc.RData")
x=scale(t(as.matrix(features)))
x[is.na(x) | is.nan(x)]=0.0

isotropic=F
n.drugs=ncol(sens$IC_50)

require(glmnet)

n.reps=10
methods=4
results=data.frame(matrix(NA,ncol=methods,nrow=n.reps))
train.rmse=data.frame(matrix(NA,ncol=methods,nrow=n.reps))
rmse=function(err) sqrt(mean(err*err))
# per.sample=list()
sens=sens$IC_50
for (cv.fold in 1:n.reps){
    meth=1
      #  per.sample[[cv.fold]]=list() 
    print(cv.fold)

    set.seed(cv.fold)

    set.seed(cv.fold)
    n.cell.lines=nrow(sens)
    perm=sample(n.cell.lines)[1:64]
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
    use.mrf=F
    truncate=T
    suffix=paste("_dnsfa",cv.fold,"alpha",alpha,"mrf",use.mrf,"truncate",truncate,sep="")

    load(paste("results",suffix,".RData",sep=""))
    p=r$best.params
    
    base.dir=paste("samples",suffix,"/",sep="")
    samples=dir(base.dir)
    y.pred=as.matrix(sens[test,]*0)
    per.sample=numeric(length(samples))
    counter=1
    for (s in samples){
        load(paste(base.dir,s,sep=""))
        if (ncol(param$B)>0){
            y.pred.1=as.matrix(t(param$factor.loadings %*% param$B %*% t(x.test)))
            y.pred=y.pred+y.pred.1
            
            per.sample[counter]=rmse((y.pred.1-sens[test,])[!is.na(sens[test,])])
            counter=counter+1
        }
    }
    y.pred=y.pred/length(samples)
    results[cv.fold,meth]=test.calc(y.pred)
    
    y.pred=p$factor.loadings %*% p$B %*% t(x.train)
    y.pred=as.matrix(t(y.pred))
    train.rmse[cv.fold,meth]=train.calc(y.pred)
    meth.name="dfa"
                                        #plot(1:length(samples),per.sample,main=meth.name,xaxt="n",lty=2,ylim=c(0,1))
                                        #axis(1,at=1:length(samples),labels=samples)
    colnames(results)[meth]=meth.name
    meth=meth+1
                                        # glmnet per drug
    print("glmnet")
    load(paste("glmnet_single",cv.fold,"IC_50.RData",sep=""))
    y.pred=sens[test,]*NA
    y.pred.train=sens[!test,]*NA
    for (i in 1:n.drugs){
        y.pred[,i]=predict(res[[i]],newx=x.test)
        y.pred.train[,i]=predict(res[[i]],newx=x.train)
    }
    results[cv.fold,meth]=test.calc(y.pred)
    train.rmse[cv.fold,meth]=train.calc(y.pred.train)
    
    colnames(results)[meth]="L1"
    meth=meth+1

    results[cv.fold,meth]=test.calc(y.adjusted[test,]*0)
    train.rmse[cv.fold,meth]=train.calc(y.adjusted[!test,]*0)
    colnames(results)[meth]="mean"
    meth=meth+1
    
}

colnames(train.rmse)=colnames(results)

save(train.rmse,results,file="results.RData")
