library(Biobase)
library(PharmacoGx)

#PSet_names = availablePSets() %>% filter(Dataset.Type=="sensitivity") %>% .$PSet.Name %>% as.character()
PSet_names=c("CCLE","GDSC1000","CTRPv2")

PSets=foreach(PSet_name=PSet_names) %do% {
  downloadPSet(PSet_name, saveDir="~/PSets/")
} %>% set_names(PSet_names)

GDSC=downloadPSet("GDSC")

pset_info=foreach(pset_name=names(PSets), .combine=bind_rows) %do% {
  pset=PSets[[pset_name]]
  data.frame(pset=pset_name, ndrug=nrow(pset@drug), ncell=nrow(pset@cell), nexp=nrow(pset@sensitivity$info), stringsAsFactors = F)
}

# 93% of the CTRPv2 cell lines are in CCLE
length(intersect(PSets$CCLE@cell$cellid, PSets$CTRPv2@cell$cellid))
length(unique(PSets$CTRPv2@cell$cellid))

my_sub=function(x, pat, repl="") { gsub(pat, repl, x) }
clean_drugs=function(n) {
  n %>% my_sub("-") %>% my_sub(" ") %>% toupper()
}
  gsub("-","",gsub(" ","",PSets$CTRPv2@drug$drugid)) %>%  toupper()

ccle_drugs=PSets$CCLE@drug$drugid %>% my_sub("drugid_") %>%
  clean_drugs()
CTRPv2_drugs = clean_drugs(PSets$CTRPv2@drug$drugid)
GDSC_drugs=clean_drugs(PSets$GDSC1000@drug$drugid)

gdsc=data.frame(gdsc=PSets$GDSC1000@drug$drugid, stringsAsFactors = F) %>% 
  mutate(clean=clean_drugs(gdsc))
data.frame(ctrp=PSets$CTRPv2@drug$drugid, stringsAsFactors = F) %>% 
  mutate(clean=clean_drugs(ctrp)) %>% 
  inner_join(gdsc) %>% 
  filter(gdsc != ctrp)

length(intersect(ccle_drugs,CTRPv2_drugs)) # 13
length(intersect(GDSC_drugs,CTRPv2_drugs)) # 74
length(union(GDSC_drugs,CTRPv2_drugs)) # 722
length(intersect(GDSC@cell$cellid,PSets$GDSC1000@cell$cellid)) # >1000
length(intersect(PSets$CTRPv2@cell$cellid,PSets$GDSC1000@cell$cellid)) # >1000

ip=intersectPSet( PSets[c(1,3)],  intersectOn=c("cell.lines", "drugs")) #, strictIntersect=TRUE)

commonGenes <- Reduce(intersect, foreach(pset=PSets[1:2]) %do% { fNames(pset, "rna") })
common <- intersectPSet( PSets[1:2],  intersectOn=c("cell.lines", "drugs"), strictIntersect=TRUE)

aucs <- foreach(pset=common) %do% { 
  summarizeSensitivityProfiles(
    pSet=pset,
    sensitivity.measure='auc_recomputed',
    summary.stat="median",
    verbose=FALSE)
} #%>% set_names(names(PSets))

ic50s <- foreach(pset=common) %do% { 
  summarizeSensitivityProfiles(
    pSet=pset,
    sensitivity.measure='ic50_published',
    summary.stat="median",
    verbose=FALSE)
}

expressions <- foreach(pset=common) %do% { summarizeMolecularProfiles(pset,
                                                 cellNames(pset),
                                                 mDataType="rna",
                                                 features=commonGenes,
                                                 verbose=FALSE) } %>% set_names(names(common))

gg <- fNames(common[[1]], 'rna')

get_cors=function(to_cor, method="pearson") {
  foreach(x=cellNames(common[[1]]), .combine = c) %do% { stats::cor(to_cor[[1]][ , x], to_cor[[2]][ , x], method=method , use="pairwise.complete.obs") } 
}

ge_cor = get_cors(foreach(e=expressions) %do% { exprs(e) })
hist(ge_cor)

ic50_cor = get_cors(ic50s)
auc_cor= get_cors(aucs)

boxplot(list("GE"=ge_cor,
                 "AUC"=auc_cor,
                 "IC50"=ic50_cor),
            main="Concordance between cell lines",
            ylab=expression(R[s]),
            sub=ss,
            ylim=yylim,
            col="lightgrey",
            pch=20,
            border="black")

CTRPv2_auc=summarizeSensitivityProfiles(
  pSet=PSets$CTRPv2,
  sensitivity.measure='auc_recomputed',
  summary.stat="median",
  verbose=FALSE)

names(PSets$CCLE@molecularProfiles)
mDataTypes=c("rna", "mutation", "cnv")

genomic_data=foreach(mDataType=mDataTypes) %dopar% { # "rnaseq"
  summarizeMolecularProfiles(PSets$CCLE, cellNames(PSets$CCLE), mDataType=mDataType, features=fNames(PSets$CCLE, mDataType), summary.stat=if(mDataType %in% c("mutation","fusion")) "or" else "median", verbose=FALSE) %>% exprs()
} %>% set_names(mDataTypes)

common_cells=intersect(colnames(CTRPv2_auc), colnames(genomic_data$rna))

y=CTRPv2_auc[,common_cells] %>% t()
genomic_data_sub = foreach(gd=genomic_data) %do% {
  class(gd)="numeric"
  s=apply(gd,1,function(g) sd(g,na.rm=T))
  gd[ order(-s)[1:1000], ]
}
x=Reduce(rbind,genomic_data_sub) %>% t() %>% scale()
x[is.na(x)]=0
y=scale(y)
y[is.na(y)]=0
cv = cv.glmnet( x[common_cells,], y, family="mgaussian", )

