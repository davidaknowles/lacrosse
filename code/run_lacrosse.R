if (interactive()){
   cv.fold=2
   use.mrf=T
   alpha=1
   truncate=T
   initial_size=50
   burn.in=10
} else {
    args <- commandArgs(trailingOnly = TRUE)
    cv.fold=as.integer(args[1])
    alpha=1
    use.mrf=T
    truncate=T
    initial_size=50
    burn.in=100
}

source("load_CTRPv2.R")

sens=y
x=scale(x)
x[is.na(x) | is.nan(x)]=0.0

set.seed(cv.fold)
n.cell.lines=nrow(sens)
perm=sample(n.cell.lines)[1:82]
test=logical(n.cell.lines)
test[perm]=T
sens[perm,]=NA

x.train=x[!test,]
y.train=scale(sens[!test,])

# y_temp=y.train
# y_temp[is.na(y_temp)]=0
# beta_ls=solve( t(x.train) %*% x.train + diag(ncol(x.train)), t(x.train) %*% y_temp )
# s=svd(beta_ls)
# plot(s$d)

suffix=paste("_dnsfa",cv.fold,"alpha",alpha,"mrf",use.mrf,"truncate",truncate,sep="")

require(lacrosse)
s=default.settings()
s$iterations=2000
s$burn.in=burn.in
s$thin.out=100
s$verbose=T
s$debug=F
s$log.file=paste("log",suffix,".txt",sep="")
s$lambda.e=1.0
s$isotropic.noise=T
s$sample.alpha=(alpha<0)
s$alpha=abs(alpha)
s$store.samples=T
s$truncate=truncate
s$initial_size=initial_size
s$samples.dir=paste(scratch_dir,"/lacrosse/samples",suffix,sep="")

require(stringr)
if (use.mrf) {
  targets_list=strsplit(PSets$CTRPv2@drug$gene_symbol_of_protein_target,";")
  s$mrf.connections=foreach(i=seq_along(targets_list)) %dopar% {
    if (length(targets_list[[i]])==0) return(numeric(0))
    setdiff( foreach(j=seq_along(targets_list), .combine=c) %do% { # need to ignore self
      if (length(intersect(targets_list[[i]],targets_list[[j]])) > 0) j else NULL
    }, i )
  }
  s$mrf.strengths=lapply(s$mrf.connections, function(temp) as.numeric(array(1.0,length(temp))))
  
  targets_list=str_split_fixed(colnames(x),":",2)[,2]
  s$feature.mrf.connections=foreach(i=seq_along(targets_list)) %dopar% {
    if (length(targets_list[i])==0 | targets_list[i]=="NA") return(numeric(0))
    setdiff( which(targets_list==targets_list[i]), i )
  }
  s$feature.mrf.strengths=lapply(s$feature.mrf.connections, function(temp) as.numeric(array(1.0,length(temp))))
}

r=run.lacrosse(y.train,x.train,settings=s)

save(r,file=paste("results",suffix,".RData",sep=""))
