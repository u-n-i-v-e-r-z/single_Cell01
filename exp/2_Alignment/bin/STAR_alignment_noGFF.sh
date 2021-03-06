#
# Copyright (C) 2017 LBME
# Schaak - Heurteau
# Single Cell project RNAseq

#!/bin/bash

# Tools path
star="/usr/local/bioinfo/src/STAR/STAR-2.5.2b/bin/Linux_x86_64/STAR"



# GET OPTION
while getopts "t:i:o:l:s:h" option
do
case $option in
    t)
        trimFQ_path=$OPTARG
        ;;
    i)
        idx_path=$OPTARG
        ;;
    o)
        out_path=$OPTARG
        ;;
    l)
        qarray_out_path=$OPTARG
        ;;
    s)
        sub_path=${OPTARG%%.*}
        ;;
    \?)
        echo "-t emplacement of trimmed reads input path eg:/my/trim/directory"
        echo "-i emplacement of the indexed reference eg:/my/idx_ref/directory"
        echo "-o where should be written samfiles eg: /my/sam/directory (created by the script if doesn't exists)"
        echo "-l log directory where error and output of qarray will be written (created by the script if doesn't exists) eg: /my/log/directory"
        echo "-s where the script that will be launched by qarray should be written eg: /my/out/directory"
        echo "-h Display this help"
        exit 1
        ;;
    h)
        echo "-t emplacement of trimmed reads input path eg:/my/trim/directory"
        echo "-i where should be written the indexed reference (created by the script if doesn't exists) eg:/my/idx_ref/directory"
        echo "-o where should be written samfiles eg: /my/sam/directory (created by the script if doesn't exists)"
        echo "-l log directory where error and output of qarray will be written (created by the script if doesn't exists) eg: /my/log/directory"
        echo "-s where the script that will be launched by qarray should be written eg: /my/out/directory"
        echo "-h Display this help"
        ;;
esac
done



# Errors
mkdir -p ${qarray_out_path}/{out_`date +%F_%H-%M`,err_`date +%F_%H-%M`}

# Alignment
mkdir -p ${out_path}
cd ${trimFQ_path}
if [ -s ${out_path}/trace.cmd ];then rm ${out_path}/trace.cmd;fi
script_name=$0
sub_dir=${sub_path}/${script_name%.*}_sub.sh
if [ -s ${sub_dir} ];then rm ${sub_dir};fi
for r1 in `ls *R1.fastq.gz`
do
  header=${r1%%_L001*} #COLXXXX
  replaceBy="R2"
  r2="${r1/R1/$replaceBy}"
  myCmd="${star} --genomeDir ${idx_path} --readFilesCommand zcat --outFileNamePrefix ${out_path}/${header} --runThreadN 24 --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimFQ_path}/${r1} ${trimFQ_path}/${r2}"
  echo -e ${myCmd} >> ${out_path}/trace.cmd
  echo -e ${myCmd} >> ${sub_dir}
done
