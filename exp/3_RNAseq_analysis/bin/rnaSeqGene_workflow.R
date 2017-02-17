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
#DM6
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
dir <- "/media/alexandre/Data/Recherche/LBME/Projects/sCell_RNAseq_clone/exp/2_Alignment/results/str_dmel613_sted_bam/"
xp_list <- unlist(lapply(read.csv("/media/alexandre/Data/Recherche/LBME/Projects/sCell_RNAseq_clone/exp/2_Alignment/results/str_dmel613_sted_bam/exp_list.csv",h=F,sep="\n"),as.character))
filenames <- file.path(dir,paste0(xp_list))
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
#seqinfo(bamfiles[1])
xp_type <- readRDS("/media/alexandre/Data/Recherche/LBME/Projects/sCell_RNAseq_clone/exp/3_RNAseq_analysis/data/sampleTable.rds")
ebg <- genes(txdb)
seqlevelsStyle(ebg) <- "Ensembl"
#Remove the ".1" suffixe on dm6 gene IDs
names(ebg) <- gsub("[.]1$","",names(ebg))
elementMetadata(ebg)$gene_id <- names(ebg)
#Create assay of #reads/gene
se_dm6 <- summarizeOverlaps(features=ebg, reads=bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
sum(unlist(lapply(bamfiles,function(x)length(readGAlignments(x)))))
#[1] 5514072

#DM3
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
dir <- "/media/alexandre/Data/Recherche/LBME/Projects/sCell_RNAseq_clone/exp/2_Alignment/results/bam_files/star_alignment_sorted/"
xp_list <- unlist(lapply(read.csv("/media/alexandre/Data/Recherche/LBME/Projects/sCell_RNAseq_clone/exp/2_Alignment/results/bam_files/star_alignment_sorted/exp_list.csv",h=F,sep="\n"),as.character))
filenames <- file.path(dir,paste0(xp_list))
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
#seqinfo(bamfiles[1])
xp_type <- readRDS("/media/alexandre/Data/Recherche/LBME/Projects/sCell_RNAseq_clone/exp/3_RNAseq_analysis/data/sampleTable.rds")
ebg <- genes(txdb)
seqlevelsStyle(ebg) <- "Ensembl"
#Create assay of #reads/gene
#se_dm3 <- summarizeOverlaps(features=ebg, reads=bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
sum(unlist(lapply(bamfiles,function(x)length(readGAlignments(x)))))
#4283967

# ACTIVE GENES LISTS
ac_genes.vec <- readRDS("/media/alexandre/Data/Recherche/LBME/Projects/sCell_RNAseq_clone/raw/Lists/Active_genes/xpaus.rds")

#####################################################################################
#
#          TREATMENT ------------
#
#####################################################################################
count_matrix_dm3 <- assay(se_dm3)
count_matrix_dm6 <- assay(se_dm6)

cnt_dm6_srted.mtx <- count_matrix_dm6[-(which())]






saveRDS(object = se_dm6,file = "/home/alexandre/Bureau/se_dm6.rds")
#log cpm (voom)
l_cpm.mtx <- log((count_matrix_dm6+0.5/sum(colSums(count_matrix_dm6))+1)*10e6)
hc <- hclust(dist(l_cpm.mtx,method="euclide"))




sum(colSums(count_matrix_dm6))
sum(colSums(count_matrix_dm3))


plot(colSums(count_matrix_dm6),colSums(count_matrix_dm3))
abline(a=1,b=1)

apply(count_matrix>0,2,sum)


# DO SOME PLOTS TO VISUALIZE DATA
nb_reads_xp <- colSums(count_matrix)
plot(nb_reads_xp,apply(count_matrix>0,2,sum),xlim=c(0,300000))

dm6_nb_gene <- nb_reads_xp
# KEEP ONLY GOOD BAM
count_matrix_sel=count_matrix[,apply(count_matrix>0,2,sum)>1000]

which(rowSums(assay(se)) == min(rowSums(assay(se))))

apply(assay(se),2,function(x)x/mean(x))

