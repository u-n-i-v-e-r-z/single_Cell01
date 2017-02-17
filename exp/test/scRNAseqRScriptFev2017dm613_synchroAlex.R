# All new with dm6 - synchro with Alex
#
# ****************************************************************************
#
# init envir, select work folder and load libs
#
# ****************************************************************************
#
cat("\014")
rm(list=ls())
#
setwd("/home/stef/Documents/New_Labo_Janv2017")
#
#
library(GenomicAlignments)
#
library(stringr)
#
library("limma")
library("edgeR")
#
library("DESeq2")
#
#
# ****************************************************************************
#
# load sample infos and conut matrix ; select samples (and create filters for selected cells)
#
# ****************************************************************************
#
sampleTable=readRDS(file="singleCellRNAseq01_janv2017/singleCellNew/data/R_Objects/sampleTableOK.rds")
#
scSovl613=readRDS(file="singleCellRNAseq01_janv2017/singleCellNew/data/R_Objects/scSovl613_OK.rds")
#
rawCountMat=assay(scSovl613)
#
#
plot(colSums(rawCountMat),colSums(rawCountMat>0),col=as.character(sampleTable$color),pch=16)
abline(v=2.5E4,col="blue",lty=2)
abline(h=1E3,col="blue",lty=2)
#
#
cellSelection=which(colSums(rawCountMat)>2.5E4)
#
rawCountMatSel11=subset(rawCountMat,select=cellSelection)
#
sampleTableSel11=subset(sampleTable,subset=rownames(sampleTable)%in%names(cellSelection))
#
sel11selectCtrl=which(sampleTableSel11$condition=="Ctrl")
sel11selectMes4KD=which(sampleTableSel11$condition=="Mes4KD")
#
barplot(colSums(rawCountMat),col=as.character(sampleTable$color),las=2)
abline(h=2.5E4,col="blue",lty=2)
#
geneFilter6990=rowSums(rawCountMat[,cellSelection]>0)>0
# length(which(geneFilter6990==TRUE))
# # [1] 6990
geneFilter4784=rowSums(rawCountMat[,sel11selectCtrl]>0)>=2 | rowSums(rawCountMat[,sel11selectMes4KD]>0)>=2
# length(which(geneFilter==TRUE))
# # [1] 4784
geneFilter3142=( rowSums(rawCountMat[,sel11selectCtrl]>0)>=2 & rowSums(rawCountMat[,sel11selectMes4KD]>0)>=2 )
# length(which(geneFilter==TRUE))
# # [1] 3142
#
#
# ****************************************************************************
#
# limma diff expr analysis
#
# ****************************************************************************
#
# library("limma")
# library("edgeR")
#
# identical(colnames(rawCountMatSel11),rownames(sampleTableSel11))
# [1] TRUE
#
(sampleCtrlFilter=sampleTableSel11$condition=="Ctrl")
(sampleMes4KDFilter=sampleTableSel11$condition=="Mes4KD")
#
(sampleCtrlSelect=which(sampleTableSel11$condition=="Ctrl"))
(sampleMes4KDSelect=which(sampleTableSel11$condition=="Mes4KD"))
#
#
dge=DGEList(counts=subset(rawCountMatSel11,subset=geneFilter3142,select=TRUE))
#
dgeNF=calcNormFactors(dge)
# 
design=model.matrix(~sampleTableSel11$condition)
# 
v01=voom(dgeNF,design=design,plot=TRUE)
#
# DiffExp analysis (method01)
fit00=lmFit(v01,design)
fit01=eBayes(fit00)
tt01=topTable(fit01,coef=ncol(design),number=dim(fit01)[1])
#
# DiffExp analysis (method02)
fit02=treat(fit01,lfc=log2(1.2))
tt02=topTreat(fit02,coef=ncol(design),number=dim(fit02)[1])
#
plotMDS(v01,col=as.character(sampleTableSel11$color),main="MDS plot")
#
# 
#
# ****************************************************************************
#
# deseq2 diff expr analysis
#
# ****************************************************************************
#
# library("DESeq2")
#
dds=DESeqDataSet(
  subset(
    scSovl613,subset=elementMetadata(rowRanges(scSovl613))@listData$gene_id%in%rownames(rawCountMat[geneFilter3142,]),
    select=rownames(colData(scSovl613))%in%rownames(sampleTableSel11)),
  design=~condition)
#
head(assay(dds))
colData(dds)
#
des=DESeq(dds)
#
ddsEsf=estimateSizeFactors(dds)
#
rLogDds=rlog(dds,blind=FALSE)
#
head(assay(rLogDds))
#
#
(deseqRes=results(des))
#
mcols(deseqRes,use.names=TRUE)
#


#
voomMat=v01$E[order(rownames(v01$E),decreasing=FALSE),]
fit01Mat=tt01[order(rownames(tt01),decreasing=FALSE),]
fit02Mat=tt02[order(rownames(tt02),decreasing=FALSE),]
#
if (identical(rownames(voomMat),rownames(fit01Mat)) & identical(rownames(voomMat),rownames(fit02Mat))) {
  exprMat=data.frame(row.names=rownames(rawCountMatSel11)[geneFilter3142])
  varMat=data.frame(row.names=rownames(rawCountMatSel11)[geneFilter3142])
  nsdMat=data.frame(row.names=rownames(rawCountMatSel11)[geneFilter3142])
}
#
exprMat$meanVoomExpr=apply(voomMat,1,mean)
exprMat$meanVoomExprCtrl=apply(voomMat,1,function(x) mean(x[sel11selectCtrl]))
exprMat$meanVoomExprMes4KD=apply(voomMat,1,function(x) mean(x[sel11selectMes4KD]))
#
exprMat$voomLogFC=fit01Mat$logFC
exprMat$voomPval=fit01Mat$P.Value
#
#
varMat$varVoomExpr=apply(voomMat,1,var)
varMat$varVoomExprCtrl=apply(voomMat,1,function(x) var(x[sel11selectCtrl]))
varMat$varVoomExprMes4KD=apply(voomMat,1,function(x) var(x[sel11selectMes4KD]))
#
varTestResList=apply(voomMat,1,function(x) var.test(x[sel11selectCtrl],x[sel11selectMes4KD]))
varMat$varFval=sapply(varTestResList,function(x) x$statistic)
varMat$varPval=sapply(varTestResList,function(x) x$p.value)
#


#
# identical(rownames(exprMat),rownames(deseqRes))
# # [1] TRUE
#
if (identical(rownames(exprMat),rownames(deseqRes))) {
  exprMat$meanDeseqExpr=deseqRes$baseMean
  exprMat$meanDeseqExprCtrl=apply(counts(des,normalized=TRUE)[,sel11selectCtrl],1,mean)
  exprMat$meanDeseqExprMes4KD=apply(counts(des,normalized=TRUE)[,sel11selectMes4KD],1,mean)
  exprMat$deseqLogFC=deseqRes$log2FoldChange
  exprMat$deseqPval=deseqRes$pvalue
}
#
#
#
# ****************************************************************************
#
# perform normality test (on voom-transformed counts)
#
# ****************************************************************************
#
shapTestRes=data.frame(row.names=rownames(voomMat))
#
# -> faire test normalite sur donées avant et après voom-transformation
#
shapTestCtrl=apply(voomMat[,sel11selectCtrl],1,shapiro.test)
shapTestMes4KD=apply(voomMat[,sel11selectMes4KD],1,shapiro.test)
#
shapTestRes$vommCtrl=sapply(shapTestCtrl,function(x) x$p.value)
shapTestRes$vommMes4KD=sapply(shapTestMes4KD,function(x) x$p.value)
#
table(apply(shapTestRes,1,min)<0.1)
# FALSE  TRUE 
# 2192   950 
#
#
# ****************************************************************************
#
# var test on voom-trasnformed data & nsd calculation
#
# ****************************************************************************






#
# calculs NSD (sur rpkm/log/? => rpm par Pascal)
#
rpmCountMat=apply(rawCountMat,2,function(x) x/sum(x)*1E6)
#
# if (identical(names(dm3genes),rownames(rawCountMat))) {
#   rpkmCountMat=sweep(rawCountMat,1,width(dm3genes)/1E3,FUN="/")
#   rpkmCountMatBAK=rpkmCountMat
#   rpkmCountMat=sweep(rpkmCountMat,2,colSums(rpkmCountMat)/1E6,FUN="/")
# }
#
#
rpmCountMatNA=replace(rpmCountMat,which(rpmCountMat==0),NA)
#
rpmCountMatNASel11=subset(rpmCountMatNA,select=cellSelection)
#
nsdMat$nsdRpmCtrl=apply(subset(rpmCountMatNASel11,subset=geneFilter3142,select=TRUE),1,function(x) sd(x[sel11selectCtrl],na.rm=TRUE)/mean(x[sel11selectCtrl],na.rm=TRUE))
nsdMat$nsdRpmMes4KD=apply(subset(rpmCountMatNASel11,subset=geneFilter3142,select=TRUE),1,function(x) sd(x[sel11selectMes4KD],na.rm=TRUE)/mean(x[sel11selectMes4KD],na.rm=TRUE))
nsdMat$deltaNsdRpm=nsdMat$nsdRpmMes4KD - nsdMat$nsdRpmCtrl
#
nsdMat$nsdVoomCtrl=apply(voomMat,1,function(x) sd(x[sel11selectCtrl],na.rm=TRUE)/mean(x[sel11selectCtrl],na.rm=TRUE))
nsdMat$nsdVoomMes4KD=apply(voomMat,1,function(x) sd(x[sel11selectMes4KD],na.rm=TRUE)/mean(x[sel11selectMes4KD],na.rm=TRUE))
nsdMat$deltaNsdVoom=nsdMat$nsdVoomMes4KD - nsdMat$nsdVoomCtrl
#
#
#
# ****************************************************************************
#
# save data
#
# ****************************************************************************
# 
# 
# saveRDS(exprMat,file="singleCellRNAseq01_janv2017/singleCellNew/R_Objects/scExprMat.rds")
# 
# saveRDS(varMat,file="singleCellRNAseq01_janv2017/singleCellNew/R_Objects/scVarMat.rds")
# 
# saveRDS(nsdMat,file="singleCellRNAseq01_janv2017/singleCellNew/R_Objects/scNsdMat.rds")
# 
# saveRDS(sampleTableSel11,file="singleCellRNAseq01_janv2017/singleCellNew/R_Objects/scSampleTableSel11.rds")
# 
# 
# ****************************************************************************
#
#
# ****************************************************************************
#

plot.new()
# plot(exprMat$meanVoomExpr,exprMat$meanDeseqExpr,col="darkgrey",pch=16,cex=0.5,xlim=c(2,16),ylim=c(0,4000))
# plot(exprMat$meanVoomExpr,log2(exprMat$meanDeseqExpr),col="darkgrey",pch=16,cex=0.5,xlim=c(2,16),ylim=c(-5,12))
plot(exprMat$meanVoomExprCtrl,log2(exprMat$meanDeseqExprCtrl),col="green4",pch=16,cex=0.5,xlim=c(2,16),ylim=c(-5,12))
# plot(exprMat$meanVoomExprMes4KD,log2(exprMat$meanDeseqExprMes4KD),col="red",pch=16,cex=0.5,xlim=c(2,16),ylim=c(-5,12))
# plot(exprMat$meanVoomExpr,exprMat$voomLogFC,col="darkorange",pch=16,cex=0.5,xlim=c(2,16),ylim=c(-5,5))
# abline(h=0,col="blue",lty=2)
# plot(exprMat$meanVoomExpr,nsdMat$deltaNsdVoom,col="darkorange",pch=16,cex=0.5,xlim=c(2,16),ylim=c(-0.7,0.7))
# abline(h=0,col="blue",lty=2)
# # vulcanos :
# 
# plot(exprMat$voomLogFC,-log(exprMat$voomPval),col="darkorange",pch=16,cex=0.5,xlim=c(-6,6),ylim=c(0,13))
# plot(exprMat$deseqLogFC,-log(exprMat$deseqPval),col="darkorange",pch=16,cex=0.5,xlim=c(-6,6),ylim=c(0,13))
# plot(exprMat$voomLogFC,-log(fit01Mat$adj.P.Val),col="darkorange",pch=16,cex=0.5,xlim=c(-6,6),ylim=c(0,5.5))
# plot(exprMat$deseqLogFC,-log(deseqRes$padj),col="darkorange",pch=16,cex=0.5,xlim=c(-6,6),ylim=c(0,5.5))
# min(exprMat$meanDeseqExprMes4KD[exprMat$meanDeseqExprMes4KD>0]/exprMat$meanDeseqExprCtrl[exprMat$meanDeseqExprMes4KD>0])
# # [1] 0.007921123
# max(exprMat$meanDeseqExprMes4KD/exprMat$meanDeseqExprCtrl)
# # [1] 37.58406
#

# saveRDS(exprMat,file="exprMat_Compar.rds")
# saveRDS(varMat,file="varMat_Compar.rds")
# 
# saveRDS(sampleTable,file="sampleTable_Compar.rds")
# saveRDS(scSovl613,file="scSovl613_Compar.rds")

