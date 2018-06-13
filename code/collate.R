
source("load_CTRPv2.R")
sens=y
x=scale(x)
x[is.na(x) | is.nan(x)]=0.0

isotropic=F
n.drugs=ncol(sens)

require(glmnet)

n.reps=10
methods=10
results=data.frame(matrix(NA,ncol=methods,nrow=n.reps))
train.rmse=data.frame(matrix(NA,ncol=methods,nrow=n.reps))
rmse=function(err) sqrt(mean(err*err))

lacrosse_dir=paste0(scratch_dir,"/lacrosse/") 
# 
# for (cv.fold in 1:n.reps){
#   fn=paste0(lacrosse_dir,"glmnet_mgaussian/",cv.fold,"_alpha1.RData")
#   a=load(fn)
#   saveRDS(., file=fn)
# }

kbtml_results=readRDS("kbtml_CTRPv2.rds")

foreach (cv.fold=1:n.reps) %do% {
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
    # correlation
    #test.calc=function(y.pred) { 
    #  mean( foreach(i=seq_len(ncol(y.pred)), .combine=c) %do% cor(y.pred[,i],y.adjusted[test,i], use="pairwise"), na.rm=T )
    #}
    
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

    colnames(results)[meth]="dfa"
    meth=meth+1
    
    # glmnet per drug
    print("glmnet")
    
    foreach(alpha=c(0,0.5,1)) %do% {
      res=readRDS(paste0(lacrosse_dir,"glmnet_single/",cv.fold,"_alpha",alpha,".RData"))
      co_matrix=lapply(res, as.numeric) %>% as.data.frame() %>% as.matrix() %>% set_colnames(NULL)
    
      y.pred=cbind(1, x.test) %*% co_matrix
      y.pred.train=cbind(1, x.train) %*% co_matrix
    
      results[cv.fold,meth]=test.calc(y.pred)
      train.rmse[cv.fold,meth]=train.calc(y.pred.train)
      
      colnames(results)[meth]=paste0("L1_",alpha)
      meth=meth+1
    } 
    
    print("glmnet mgaussian")
    foreach(alpha=c(0,0.5,1)) %do% {
    
      res=readRDS(paste0(lacrosse_dir,"glmnet_mgaussian/",cv.fold,"_alpha",alpha,".RData"))
      co_matrix=lapply(res, as.numeric) %>% as.data.frame() %>% as.matrix() %>% set_colnames(NULL)
    
      y.pred=cbind(1, x.test) %*% co_matrix
      y.pred.train=cbind(1, x.train) %*% co_matrix
    
      results[cv.fold,meth]=test.calc(y.pred)
      train.rmse[cv.fold,meth]=train.calc(y.pred.train)
    
      colnames(results)[meth]=paste0("mgaussian_",alpha)
      meth=meth+1
    }
    
    print("mvlr")
    y.pred=readRDS(paste0(lacrosse_dir,"mvlr/",cv.fold,".RData"))
    results[cv.fold,meth]=test.calc(y.pred)
    colnames(results)[meth]="MVLR"
    meth=meth+1
  
    print("kbtml")
    results[cv.fold,meth]=test.calc(kbtml_results[[cv.fold]]$Y$mu)
    colnames(results)[meth]="KBMTL"
    meth=meth+1

    results[cv.fold,meth]=test.calc(y.adjusted[test,]*0)
    train.rmse[cv.fold,meth]=train.calc(y.adjusted[!test,]*0)
    colnames(results)[meth]="mean"
    meth=meth+1
    
}

colnames(train.rmse)=colnames(results)

foreach(i=seq_len(ncol(results)), .combine = c) %do% { t.test(results[,i], results$dfa, paired = T)$p.value }

my_gsub=function(g,a) {
  for(i in seq_along(a)) { g=gsub(names(a)[i],a[i],g) }
  g
}

theme_set(theme_bw(base_size = 14))

results %>% 
  gather(method, cor, -mean) %>% 
   ggplot(aes(method, cor)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

to_plot = results %>% set_colnames(colnames(results) %>% 
                           my_gsub(c("dfa"="Lacrosse",
                                     "kbmtl"="KBMTL", 
                                     "mvlr"="MVLR", 
                                     "L1_0.5"="Elastic Net",
                                     "L1_1"="Lasso",
                                     "L1_0"="Ridge", 
                                     "mgaussian_0.5"="Group ENet",
                                     "mgaussian_0"="Ignore", # should be same as Ridge, is worse for some reason? 
                                     "mgaussian_1"="Group Lasso"))) %>% 
  mutate(fold=1:n()) %>%
  gather(method, rmse, -mean, -fold) %>% 
  filter(method!="Ignore") %>% 
  mutate(pve=1-rmse/mean, 
         method=factor(method, c("Ridge", "Elastic Net", "Lasso", "Group ENet", "Group Lasso", "KBMTL", "MVLR", "Lacrosse"))) 
  
to_plot %>%   
  ggplot(aes(method, pve * 100, label=fold)) + geom_boxplot(outlier.shape = NA) + 
  geom_text(position = position_jitter(width=0.2)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylab("Percent variance explained") + xlab(NULL)
ggsave("../figures/CTRPv2_pve.pdf", width=6, height=4)

to_plot %>% mutate(pve=100*pve) %>% group_by(method) %>% summarize(m=mean(pve),s=sd(pve)) %>% 
  ggplot(aes(method, m, ymin=m-s, ymax=m+s)) + geom_bar(stat="identity") + geom_errorbar(width=0.5) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylab("Percent variance explained")  + scale_y_continuous(breaks=c(2.5,5,7.5,10))+ xlab(NULL) + coord_cartesian(ylim=c(2.5,10)) 
ggsave("../figures/CTRPv2_pve_clean.pdf", width=5, height=4)

t.test(results$L1_0, results$mgaussian_0, paired = T)
t.test(results$L1_0.5, results$mgaussian_0, paired = T)

