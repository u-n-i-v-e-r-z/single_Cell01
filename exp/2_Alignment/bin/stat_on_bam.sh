#
# Copyright (C) 2017 LBME
# Schaak - Heurteau
# Single Cell project RNAseq


#!/bin/bash
#Toolspath
samTools="/usr/local/bioinfo/src/samtools/samtools-1.3/samtools"

# GET OPTION
while getopts "b:o:l:s:h" option
do
case $option in
    b)
        bam_path=$OPTARG
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
        echo "-b emplacement of .bam input files eg:/my/bamfiles/directory"
        echo "-o where should be written sorted bamfiles eg: /my/sorted_bamfiles/directory (created by the script if doesn't exists)"
        echo "-l log directory where error and output of qarray will be written (created by the script if doesn't exists) eg: /my/log/directory"
        echo "-s where the script that will be launched by qarray should be written eg: /my/out/directory"
        echo "-h Display this help"
        exit 1
        ;;
    h)
        echo "-b emplacement of .bam input files eg:/my/bamfiles/directory"
        echo "-o where should be written sorted bamfiles eg: /my/sorted_bamfiles/directory (created by the script if doesn't exists)"
        echo "-l log directory where error and output of qarray will be written (created by the script if doesn't exists) eg: /my/log/directory"
        echo "-s where the script that will be launched by qarray should be written eg: /my/out/directory"
        echo "-h Display this help"
        ;;
esac
done

#Log path
mkdir -p ${qarray_out_path}/{out_`date +%F_%H-%M`,err_`date +%F_%H-%M`}
mkdir -p ${out_path}
script_name=$0
sub_dir=${sub_path}/${script_name%.*}_sub.sh
if [ -s ${sub_dir} ];then rm ${sub_dir};fi
cd ${bam_path}
for bamFile in $(ls *.bam)
do
  runName=${bamFile%%.*}
  echo ${samTools} flagstat ${bam_path}/${bamFile} \> ${out_path}/${runName}.stat >> ${out_path}/trace
  echo ${samTools} flagstat ${bam_path}/${bamFile} \> ${out_path}/${runName}.stat >> ${sub_dir}
done
