#
# Copyright (C) 2017 LBME
# Schaak - Heurteau
# Single Cell project RNAseq


#First we need to index the genome


g_dir="/home/aheurteau/work/LBME/single_Cell01/exp/2_Alignment/raw/dmel_541/dmel-all-chromosome-r5.41.fasta"
o_dir="/home/aheurteau/work/LBME/single_Cell01/exp/2_Alignment/raw/dmel_541/"

STAR  --runMode genomeGenerate --runThreadN 4 --genomeDir ${o_dir}  --genomeFastaFiles ${g_dir}


# then create GTF file from transcript file .GFF

gffFile="/home/aheurteau/work/LBME/single_Cell01/exp/2_Alignment/raw/dmel_541/GFF/dmel-all-r5.41.gff"
gtfFile="/home/aheurteau/work/LBME/single_Cell01/exp/2_Alignment/raw/dmel_541/GTF/dmel-all-r5.41.gtf"

gffread ${gffFile} -T -o ${gtfFile}
