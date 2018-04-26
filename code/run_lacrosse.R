if (interactive()){
   cv.fold=2
   use.mrf=F
   alpha=1
   truncate=T
} else {
    args <- commandArgs(trailingOnly = TRUE)
    cv.fold=as.integer(args[1])
    alpha=as.numeric(args[2])
    use.mrf=as.logical(as.integer(args[3]))
    truncate=as.logical(as.integer(args[4]))
}

source("load_CTRPv2.R")

sens=t(y)
set.seed(cv.fold)
n.cell.lines=nrow(sens)
perm=sample(n.cell.lines)[1:55]
test=logical(n.cell.lines)
test[perm]=T
sens[perm,]=NA
x=scale(x)
x[is.na(x) | is.nan(x)]=0.0
x.train=x[!test,]
y.train=scale(sens[!test,])

suffix=paste("_dnsfa",cv.fold,"alpha",alpha,"mrf",use.mrf,"truncate",truncate,sep="")

require(lacrosse)
s=default.settings()
s$iterations=2000
s$burn.in=1000
s$thin.out=100
s$verbose=T
s$debug=T
s$truncate=truncate
s$log.file=paste("log",suffix,".txt",sep="")
s$lambda.e=1.0
s$isotropic.noise=T
s$sample.alpha=(alpha<0)
s$alpha=abs(alpha)
s$store.samples=T
s$truncate=T
s$initial_size=10
s$samples.dir=paste("samples",suffix,sep="")
if (use.mrf) {
    source("calculate.mrfs.R")
    s$mrf.connections=drug.mrf(colnames(sens),0)
    s$mrf.strengths=lapply(s$mrf.connections, function(temp) as.numeric(array(1.0,length(temp))))
    s$feature.mrf.connections=gene.mrf()
    s$feature.mrf.strengths=lapply(s$feature.mrf.connections, function(temp) as.numeric(array(1.0,length(temp))))
}

r=run.lacrosse(y.train,x.train,settings=s)

save(r,file=paste("results",suffix,".RData",sep=""))
