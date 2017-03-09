'''
David DEPIERRE
02/03/2017



data :  DHMG sur k27 (deseq2 without replciates)
			-2000/0 proximal promoter
			-1000/0
			-500/0 core pormoter
			-2000/1500
			0/1000 genebody
		DE analysis on RNA seq

'''

library(rtracklayer)
library(GenomicRanges)
library(pROC)
library(ggplot2)
require(grid)

k9_Mes4_m500 = readRDS("/home/ddepierre/work/David_M2_K9K27DE/DHMG/DHMG_analysis_k9me3/DHMG_limma_FROMm500b_TOp0b/output_matrix_Limma_chipseq_H3K9me3_FROMm500b_TOp0b_control_Mes4_KD.rds")
de_Mes4 = readRDS("/home/ddepierre/work/David_M2_K9K27DE/DEAnalysis/DE_analysis/output_matrix_Limma_rnaseq_Diff_Exp_RNAseq_control_Mes4_KD.rds")

k9_Mes4_m500 = readRDS("/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/DHMG/DHMG_analysis_k9me3/DHMG_limma_FROMm500b_TOp0b/output_matrix_Limma_chipseq_H3K9me3_FROMm500b_TOp0b_control_Mes4_KD.rds")
de_Mes4 = readRDS("/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/DEAnalysis/DE_analysis/output_matrix_Limma_rnaseq_Diff_Exp_RNAseq_control_Mes4_KD.rds")

de_mat = de_Mes4
dhmg_mat = k9_Mes4_m500

de_mat = de_mat[de_mat[,"logFC"]<log2(1),]

de_matsub = de_mat[rownames(de_mat)%in%rownames(dhmg_mat),]

matTruePos = de_matsub[,"logFC"]
names(matTruePos) = rownames(de_matsub)
matTruePos = matTruePos[order(matTruePos)]

fc_dhmg = dhmg_mat[names(matTruePos),"logFC"]

matTruePos = cbind(matTruePos, fc_dhmg)
colnames(matTruePos) = c("de_fc","dhmg_fc")
#~ vCursor = dhmg_mat[dhmg_mat[,"logFC"]> log2(1), "logFC"] 
#~ names(vCursor) = rownames(dhmg_mat[dhmg_mat[,"logFC"]> log2(1), ])
#~ vCursor = vCursor[order(vCursor, decreasing=T)]

#~ AUCmat = matrix(NA,ncol = 100, nrow= 100, 
#~ 				dimnames = list(seq(max(matTruePos[,"dhmg_fc"]), 0, length.out = 11)[-1],
#~ 								seq(min(matTruePos[,"de_fc"]), 0, length.out = 11)[-1]))

AUCmat = matrix(NA,ncol = dim(matTruePos)[1], nrow= dim(matTruePos)[1], 
				dimnames = list(matTruePos[order(matTruePos[,"dhmg_fc"], decreasing=T),"dhmg_fc"],
								matTruePos[,"de_fc"]))



sapply(1:dim(AUCmat)[2], function(idx_fcDE){
#~ 	cat(paste0(idx_fcDE," = ", as.double(colnames(AUCmat)[idx_fcDE]), " and "))
	sapply(1:dim(AUCmat)[1], function(idx_fcDHMG) {
#~ 		cat(paste0(idx_fcDHMG, " = ",as.double(rownames(AUCmat)[idx_fcDHMG]),"\n"))
		TFvector = sapply(matTruePos[matTruePos[,"de_fc"] < as.double(colnames(AUCmat)[idx_fcDE]),"dhmg_fc"], function(x) {
			as.numeric(x  >= as.double(rownames(AUCmat)[idx_fcDHMG]))
		})
		if(length(table(TFvector))-1){
			rocC = roc(TFvector, as.numeric(matTruePos[c(1:length(TFvector)),"de_fc"]), direction=">")
			AUCmat[idx_fcDHMG, idx_fcDE] = rocC$auc
		}
	})
})



TFvector = sapply(matTruePos[matTruePos[,"de_fc"] < as.double(colnames(AUCmat)[1]),"dhmg_fc"], function(x) {
			as.numeric(x  >= as.double(rownames(AUCmat)[1]))})



TFvector = sapply(matTruePos[matTruePos[,"de_fc"] < as.double(colnames(AUCmat)[2]),"dhmg_fc"], function(x) {
			as.numeric(x  >= as.double(rownames(AUCmat)[5]))})




test = sapply(vCursor, function(cursor){
	sapply(matTruePos[,"dhmg_fc"], function(x) {
	as.numeric(x  >= cursor)
	})
})

saveRDS(test, file="/home/ddepierre/work/David_M2_K9K27DE/AUC_foldchange/matTP.rds")

test = readRDS("/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/AUC_foldchange/matTP.rds")

pdf("/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/AUC_foldchange/numberOfTrue.pdf")
plot(seq(1:length(colSums(test))),colSums(test))
plot(seq(1:length(vCursor)), vCursor)
plot(roc(as.numeric(test[,1]), as.numeric(matTruePos[,"de_fc"]), direction=">"),col="yellow", lwd=3, main="Validation ROC")
plot(roc(as.numeric(test[,10]), as.numeric(matTruePos[,"de_fc"]), direction=">"),col="yellow", lwd=3, main="Validation ROC")
plot(roc(as.numeric(test[,50]), as.numeric(matTruePos[,"de_fc"]), direction=">"),col="yellow", lwd=3, main="Validation ROC")
plot(roc(as.numeric(test[,200]), as.numeric(matTruePos[,"de_fc"]), direction=">"),col="yellow", lwd=3, main="Validation ROC")
plot(roc(as.numeric(test[,500]), as.numeric(matTruePos[,"de_fc"]), direction=">"),col="yellow", lwd=3, main="Validation ROC")
plot(roc(as.numeric(test[,1000]), as.numeric(matTruePos[,"de_fc"]), direction=">"),col="yellow", lwd=3, main="Validation ROC")

dev.off()


