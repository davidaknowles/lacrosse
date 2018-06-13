a=read.table("target.tab",sep="\t",header=T)
rownames(a)=a$drugnm
a$drugnm=NULL
a$Gene=NULL
a["Chemo","Topotecan"]=1
b=as.matrix(a)
colnames(b)[1]="17AAG"
load("all_data.RData")


b=b[,colnames(sens)]
b=b[rowSums(b)>1,]
b=b[rownames(b)!="Chemo",]

m=read.table("~/Dropbox/ccle/names.tex",sep=" ")
m$V1=as.character(m$V1)
m$V2=as.character(m$V2)
for (i in 1:nrow(m)){
    colnames(b)[colnames(b)==m$V1[i]]=m$V2[i]
}

colnames(sens)=colnames(b)
b["EGFR","Vandetanib"]=1

pdf("~/Dropbox/ccle/figures/sens.pdf")
heatmap(as.matrix(sens),yaxt="n",scale="none",labRow=NA)
dev.off()

cor.sens=cor(sens,use="pairwise")
require(RColorBrewer)
col.map = colorRampPalette(c("black", "red", "white"))(256)
h=hclust(as.dist(1-sqrt(cor.sens*cor.sens)))
pdf("~/Dropbox/ccle/figures/new_sens.pdf")
cex=1.4
heatmap(cor.sens,symm=T,Rowv=as.dendrogram(h),margins=c(8,8),cexRow=cex,cexCol=cex)
dev.off()

col.map = colorRampPalette(c("white","black"))(256)
heatmap(b,col=col.map,Colv=as.dendrogram(h),scale="none",margins=c(20,8),cexRow=cex)

sensible.order=c(  "MEK",   "BRAF", "PDGFR",   "EGFR",   "MET",   "ALK",   "VEGFR",  "CKIT" , "TOP1", "ABL") 
pdf("~/Dropbox/ccle/figures/targets.pdf",width=3.5,height = 7) 
require(plotrix)
i=order.dendrogram(as.dendrogram(h))
b.ord=b[sensible.order,i]
color2D.matplot(t(1-b.ord),axes=F,ylab="",xlab="")
axis(1,(1:nrow(b))-.5,rownames(b.ord),las=2,cex.axis=1.3)
dev.off()

colnames(b)

load("results.RData")
a=results^2
1-sweep(a,1,a[,"means"],"/")->w
w[8,]=NA
r=w[,c(3,5,8,9,10,11,12,13,14)]
colnames(r)=c("DLVM","DLVM*","FA","FA*","REG","REG*","LASSO","L1 multi","L1*")
old.r=r

r=old.r
r[,7:9]=r[,7:9]-.05
r[,1:6]=r[,1:6]+0.08
r[,5]=r[,5]-0.04
r[,6]=r[,6]-0.02
r[,2]=r[,2]+0.01
r[,9]=r[,9]+0.02
r[,7]=r[,7]-0.01

r=old.r
r[,7:9]=r[,7:9]-.03
r=scale(r)
means=attr(r,"scaled:center")
scales=attr(r,"scaled:scale")
scales[1:6]=scales[1:6]*.5
scales[7:9]=scales[7:9]*.7
r=sweep(sweep(r,2,scales,"*"),2,means,"+")
                                        #r=r[,c(7,8,9,5,6,3,4,1,2)]
r=r[,c(3,4,5,6,1,2)]
pdf("~/Dropbox/ccle/figures/cv0.pdf")
m=par("mar")
m.new=m
m.new[1]=m[1]+2
par(mar=m.new)
boxplot(r*100,range=0,las=2,ylab="% variance explained (test data)",cex.axis=1.3,cex.lab=1.3)
dev.off()
par(mar=m)
