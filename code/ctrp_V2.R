
require(doMC)
registerDoMC(6)
require(glmnet)

source("load_CTRPv2.R")

x= scale(x)
x[is.na(x)]=0
y=scale(y)
y[is.na(y)]=0
cv = cv.glmnet( x[common_cells,], y, family="mgaussian", parallel=T, nfolds = 6)

require(rrpack)
rrr_results=rrr(y,x[common_cells,])

# crashes
sofar_results=sofar(y,x[common_cells,], nrank = 10)
