source("load_CTRPv2.R")
require(glmnet)

require(doMC)
#if (!interactive())  { registerDoMC(20) }
if (!interactive()) {
  registerDoMC(10)
  args <- commandArgs(trailingOnly = TRUE)
  cv.fold=as.integer(args[1])
} else { cv.fold=1 }

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

resdir=paste0(scratch_dir,"/lacrosse/glmnet_mgaussian/")
dir.create(resdir)

y.train[is.na(y.train)]=0
cv.glmnet( x.train, y.train, family="mgaussian", parallel=!interactive(), nfolds = 10) %>% 
	   coef() %>% 
	   save(file=paste0(resdir,cv.fold,".RData"))


