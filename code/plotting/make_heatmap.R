source("load_CTRPv2.R")

sens=y
feature.names=colnames(x)

require(RColorBrewer)
if (1){
    load("results_dnsfa_all1_alpha1mrfTRUEtruncateTRUE.RData")
    f=r$best.params$factor.loadings
    to.keep=colSums(f!=0.0)>1
    f=f[,to.keep]
    b=r$best.params$B[to.keep,]
} else {
    load("samples_dnsfa_long_alpha1mrfTRUEtruncateTRUE/sample3900.RData")
    f=as.matrix(param$factor.loadings)
    to.keep=colSums(f!=0.0)>0
    f=f[,to.keep]
    b=as.matrix(param$B[to.keep,])
}
col.map = colorRampPalette(c("black", "blue"))(256)


#f[,2:ncol(f)]=-f[,2:ncol(f)]
#to.keep=colSums(f!=0)>1
#f=f[,to.keep]
#
#b=b[to.keep,]
b=b[ncol(f):1,]
#swap.f=f[,ncol(f):1]
swap.f=f
colnames(swap.f)=1:ncol(f)
for (i in 1:ncol(f)){
    true.i=ncol(f)-i+1
    cat("--------",true.i,"---\n")
    bf=b[i,b[i,]!=0]
    s=sort.int(abs(bf),index.return = T,decreasing=T)
    bf=bf[s$ix]
    fn=feature.names[b[i,]!=0]
    fn=fn[s$ix]
    mystr=paste(true.i,": ",sep="")
    to.print=min(6,length(fn))
    for (j in 1:to.print)
        {
            cat(fn[j],":",bf[j],"\n")
            if (length(grep("orf",fn[j]))==0){
                mystr=paste(mystr,fn[j],sep="")
                if (j<to.print)
                    mystr=paste(mystr,", ",sep="")
            }
       }
    colnames(swap.f)[[i]]=mystr
}
#colnames(swap.f)=ncol(f):1
den=as.dendrogram(hclust(dist(abs(swap.f)),method="ward"))
                                        #oldm=par("mar")
num_drugs=ncol(sens)
p.values=numeric(num_drugs)
ests=numeric(num_drugs)
for (i in 1:num_drugs) {
    test=cor.test(sens[,i],ge[,"FBXW8"],method="spearman")
    p.values[i]=test$p.value
    ests[i]=test$estimate
}

from_0_pv=foreach(i=1:ncol(f), .combine = c) %do% {
  t.test(f[,i]) %$% p.value
}
most_non_0=which.min(from_0_pv)

clustersDir="clusters/"
dir.create(clustersDir)
maxmax=0
require(lattice)
colnames(x)=feature.names
maxlog10p=10
#niceCols <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
niceCols <- colorRampPalette(c("blue", "white", "red"), space = "rgb")(100)
meta=foreach (i=1:ncol(swap.f), .combine = bind_rows) %do% {
    bf=b[i,b[i,]!=0]
    s=sort.int(abs(bf),index.return = T,decreasing=T)
    bf=bf[s$ix]
    fn=feature.names[b[i,]!=0]
    fn=fn[s$ix]
    print(fn)
    drugsInCluster=rownames(swap.f)[ abs(swap.f[,i]) >= 1]
    pvals=matrix(NA,length(fn),length(drugsInCluster))
    rownames(pvals)=fn
    colnames(pvals)=drugsInCluster
    ests=pvals
    for (feat in fn)
        for (drug in drugsInCluster){
            ct=cor.test(x[,feat],sens[,drug],method="spearman") 
            pvals[feat,drug]=ct$p.value
            ests[feat,drug]=ct$estimate
        }
    if (nrow(pvals)==0) return(NULL)
    pdf(paste0(clustersDir,"plot",i,".pdf"),width=2+.2*length(drugsInCluster),height=10.5+.25*length(fn))
    y=-log10(t(pvals))
    ests=t(ests)
    yscaled=pmin(y,maxlog10p)
    yscaled[ests<0]=-yscaled[ests<0]
    rownames(yscaled)=rownames(y)
    colnames(yscaled)=colnames(y)
    latticePlot=levelplot(yscaled,col.regions = niceCols,scales=list(x=list(rot=90)),xlab="",ylab="",at=seq(-maxlog10p,maxlog10p,length.out = 100))
    print(latticePlot)
    dev.off()
    data.frame(max_nlp=max(-log10(pvals)), num_drugs=length(drugsInCluster), num_features=length(fn))
}
rownames(swap.f)=colnames(sens)

m=read.table("~/Dropbox/ccle/names.tex",sep=" ")
m$V1=as.character(m$V1)
m$V2=as.character(m$V2)
for (i in 1:nrow(m)){
    rownames(swap.f)[rownames(swap.f)==m$V1[i]]=m$V2[i]
    colnames(sens)[colnames(sens)==m$V1[i]]=m$V2[i]
}

                                        #m=oldm
#m[3]=m[3]+2.0
#par(mar=m)
fig.dir="~/Dropbox/ccle/figures/"
#pdf(paste(fig.dir,"temp.pdf",sep=""),width=10,height=7)
heatmap((abs(t(swap.f)))^(1/5),col=col.map,scale="none",Rowv=NA,Colv=den,margins=c(16,12),cexRow=1.3,cexCol=1.1)
#dev.off()
#system("pdfcrop ~/Dropbox/ccle/figures/temp.pdf")
#system("mv ~/Dropbox/ccle/figures/temp-crop.pdf factor_loadings.pdf")

a=read.table("target.tab",sep="\t",header=T)
rownames(a)=a$drugnm
a$drugnm=NULL
a$Gene=NULL
a["Chemo","Topotecan"]=1
a=as.matrix(a)
a=a[rowSums(a)>1,]
colnames(a)[1]="17AAG"
a=a[,rownames(f)]

b=a[,p.values<0.002]
rowSums(b)/rowSums(a)


sensible.order=c(  "MEK",   "BRAF", "PDGFR",   "EGFR",   "MET",   "ALK",   "VEGFR",  "CKIT" , "TOP1", "ABL") 
pdf("~/Dropbox/ccle/figures/targets.pdf",width=3.5,height = 7) 
require(plotrix)
i=order.dendrogram(den)
b.ord=b[sensible.order,i]
color2D.matplot(t(1-b.ord),axes=F,ylab="",xlab="")
axis(1,(1:nrow(b))-.5,rownames(b.ord),las=2,cex.axis=1.3)
dev.off()

col.map = colorRampPalette(c("white", "black"))(2)
#pdf("~/Dropbox/ccle/figures/targets.pdf")
heatmap(a,Colv=den,col=col.map,Rowv=NA,margins=c(18,5))
#dev.off()

                               
