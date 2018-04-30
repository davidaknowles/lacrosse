source("load_CTRPv2.R")
args <- commandArgs(trailingOnly = TRUE)
cv.fold=as.integer(args[1])

require(doMC)
if (!interactive())  { registerDoMC(20) }

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

resdir=paste0(scratch_dir,"/lacrosse/glmnet_single/")
dir.create(resdir)

require(glmnet)

for(alpha in c(0,1)) {
  foreach (i=1:ncol(sens)) %dopar% {
      print(i)
      y=y.train[,i]
      to.keep=! (is.nan(y) | is.na(y))
      cv.glmnet(x.train[to.keep,],y[to.keep],alpha=alpha) %>% coef()
  } %>% saveRDS( file=paste0(resdir,cv.fold,"_alpha",alpha,".RData" ) )
}

