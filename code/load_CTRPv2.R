library(PharmacoGx)
library(Biobase)
require(tidyverse)
require(magrittr)
require(doMC)

#PSet_names = availablePSets() %>% filter(Dataset.Type=="sensitivity") %>% .$PSet.Name %>% as.character()
PSet_names=c("CCLE","CTRPv2")

#scratch_dir=Sys.getenv("SCRATCH")
scratch_dir="/scratch/users/dak33/"
#scratch_dir="~"
print(scratch_dir)

PSets=foreach(PSet_name=PSet_names) %do% {
  downloadPSet(PSet_name, saveDir=paste0(scratch_dir,"/PSets/"))
} %>% set_names(PSet_names)

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

rownames( genomic_data$rna ) =PSets$CCLE@molecularProfiles$rna@featureData@data$Symbol

y=CTRPv2_auc[,common_cells] %>% t()
genomic_data_sub = foreach(gdn=names(genomic_data)) %do% {
  gd=genomic_data[[gdn]]
  rownames(gd)=paste(gdn,rownames(gd),sep=":")
  class(gd)="numeric"
  s=apply(gd,1,function(g) sd(g,na.rm=T))
  gd[ order(-s)[1:1000], ]
} %>% set_names(mDataTypes)
x=Reduce(rbind,genomic_data_sub) %>% t()
x=x[common_cells,]


