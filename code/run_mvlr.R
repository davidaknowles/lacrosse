if (interactive()){
  cv.fold=2
} else {
  args <- commandArgs(trailingOnly = TRUE)
  cv.fold=as.integer(args[1])
}

scratch_dir="/scratch/users/dak33/"
#scratch_dir="~"
print(scratch_dir)

resdir=paste0(scratch_dir,"/lacrosse/mvlr/")
dir.create(resdir)
resfile=paste0(resdir,cv.fold,".RData")

if (file.exists(resfile)) {
  cat(resfile," already exists, exiting\n")
  quit()
}

source("load_CTRPv2.R")

source("mvlr.R")

sens=y
set.seed(cv.fold)
n.cell.lines=nrow(sens)
perm=sample(n.cell.lines)[1:82]
test=logical(n.cell.lines)
test[perm]=T
sens[perm,]=NA
x=scale(x)
x[is.na(x) | is.nan(x)]=0.0
x.train=x[!test,]
y.train=scale(sens[!test,])

remove_na=function(g) { g[is.na(g)]=0.0; g }
#foreach(dat=genomic_data_sub) %do% { dat[,!test] %>% t() %>% scale() %>% remove_na()  }

M=3
list(X=x.train,
     dX=ncol(x.train) %>% as.integer(),
     Y=remove_na(y.train), 
     nY = nrow(x.train) %>% as.integer(),
     M=M %>% as.integer(),
     dY=ncol(y.train) %>% as.integer(),
     lambda_v = 1e-1,
     tau_v=1,
     wp=array(1,dim=ncol(y.train)),
     taup=array(1,dim=M),
     dXm=array(1000,dim=M) %>% as.integer()) %>% # dimension in each view
  runMVSTAN() %>% 
  predictMVSTAN(x[test,]) %>% 
  saveRDS(file=resfile)
