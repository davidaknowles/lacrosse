

set.seed(cv.fold)
n.cell.lines=nrow(sens)
perm=sample(n.cell.lines)[1:64]
test=logical(n.cell.lines)
test[perm]=T
sens[perm,]=NA
x=scale(t(as.matrix(features)))
x[is.na(x) | is.nan(x)]=0.0
x.train=x[!test,]
y.train=scale(sens[!test,])

require(glmnet)
res=list()

for (i in 1:ncol(sens)){
    print(i)
    y=y.train[,i]
    to.keep=! (is.nan(y) | is.na(y))
    y.nonan=y[to.keep]
    x.nonan=x.train[to.keep,]
    alpha=.8
    cv.res=cv.glmnet(x.nonan,y.nonan,alpha=alpha)
    res[[i]]=glmnet(x.nonan,y.nonan,lambda=cv.res$lambda.min,alpha=alpha)
}

save(res,file=paste("glmnet_single",cv.fold,measure,".RData",sep=""))


