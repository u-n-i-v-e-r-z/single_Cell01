# PROJET SINGLE CELL 1
# Alexandre Heurteau
# Cuvier's team
# 02/02/17


#Alignement realized with STAR with default opts 
# RNAseq analysis realized with rnaseqGene workflow :
# http://www.bioconductor.org/help/workflows/rnaseqGene/


#####################################################################################
#
#          LOAD LIBRARIES ------------
#
#####################################################################################
source("http://bioconductor.org/workflows.R")
workflowInstall("rnaseqGene")
biocLite("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(GenomicAlignments)
library("BiocParallel")
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library("Rsamtools")
#####################################################################################
#
#          SET PATH & LOAD DATA ------------
#
#####################################################################################
TxDb.Dmelanogaster.UCSC.dm6.ensGene
dir <- "/media/alexandre/Data/Recherche/LBME/Projects/sCell_RNAseq_clone/exp/2_Alignment/results/bam_files/star_alignment_sorted/"
xp_list <- unlist(lapply(read.csv("/media/alexandre/Data/Recherche/LBME/Projects/sCell_RNAseq_clone/exp/2_Alignment/results/bam_files/star_alignment_sorted/exp_list.csv",h=F,sep="\n"),as.character))
filenames <- file.path(dir,paste0(xp_list))
file.exists(filenames)

bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
#####################################################################################
#
#          TREATMENT ------------
#
#####################################################################################
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
ebg <- exonsBy(txdb, by="gene")

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE)


count_matrix <- assay(se)

# DO SOME PLOTS TO VISUALIZE DATA
nb_reads_xp <- colSums(count_matrix)
plot(nb_reads_xp,apply(count_matrix>0,2,sum))
# KEEP ONLY GOOD BAM
count_matrix_sel=count_matrix[,apply(count_matrix>0,2,sum)>1000]





which(rowSums(assay(se)) == min(rowSums(assay(se))))

apply(assay(se),2,function(x)x/mean(x))

