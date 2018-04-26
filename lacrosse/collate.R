lasso=T
n.reps=1
df=data.frame(errs=numeric(n.reps)*NA,gb=numeric(n.reps))
for (seed in 1:n.reps){
	load(paste("~/scratch/dnsfa_tests/test",seed,if (lasso) "lasso" else "",".RData"))
	gb=res$true.G %*% res$true.B
	gbi=res$res$best.params$factor.loadings %*% res$res$best.params$B
	its=length(res$res$test.rmse)
	df$errs[seed]=mean(res$res$test.rmse[(its-100):its])
	df$gb[seed]=mean(abs(gb-gbi))
}