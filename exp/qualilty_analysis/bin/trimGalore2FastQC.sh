# Single Cell RNA seq project
# Cuvier's Team
# Schaak - Heurteau
# 2017


# Apply trimGalore and Fastqc on new trimmed sequences recursively on all directories under the current one

#!/bin/bash

# Tools path
trimG="/usr/local/bioinfo/src/Trim_Galore/trim_galore_v0.4.0/trim_galore"

# GET OPTION
while getopts "f:g:o:l:s:h" option
do
case $option in
    f)
        fastq_path=$OPTARG
        ;;
    g)
        out_path=$OPTARG
        ;;
    o)
        out_trim_path=$OPTARG
        ;;
    l)
        qarray_out_path=$OPTARG
        ;;
    s)
        sub_path=${OPTARG%%.*}
        ;;
    \?)
        echo "-f emplacement of fastq files eg:/myfastq/directory"
        echo "-o path to output trimmed reads (created by the script if doesn't exists) eg:/my/path_to_trimmed"
        echo "-g path to graphics fastqc output (created by the script if doesn't exists) eg:/my/path_to_fastqc_graph"
        echo "-l log directory where error and output of qarray will be written (created by the script if doesn't exists) eg: /my/log/directory"
        echo "-s where the script that will be launched by qarray should be written eg: /my/out/directory"
        echo "-h Display this help"
        exit 1
        ;;
    h)
        echo "-f emplacement of fastq files eg:/myfastq/directory"
        echo "-o path to output trimmed reads (created by the script if doesn't exists) eg:/my/path_to_trimmed"
        echo "-g path to graphics fastqc output (created by the script if doesn't exists) eg:/my/path_to_fastqc_graph"
        echo "-l log directory where error and output of qarray will be written (created by the script if doesn't exists) eg: /my/log/directory"
        echo "-s where the script that will be launched by qarray should be written eg: /my/out/directory"
        echo "-h Display this help"
        ;;
esac
done


# Qarray Errors
mkdir -p ${qarray_out_path}/{out_`date +%F_%H-%M`,err_`date +%F_%H-%M`}

for pool in $(ls ${fastq_path})
do
if [ -d ${fastq_path}/${pool} ]
then
	mkdir -p ${out_path}/${pool}
	mkdir -p ${out_trim_path}/${pool}
	cd ${fastq_path}/${pool}
	if [ -s ${out_path}/${pool}/trace.cmd ];then rm ${out_path}/${pool}/trace.cmd;fi
	for r1 in `ls *R1*`
	do
	replaceBy="R2"
	r2="${r1/R1/$replaceBy}"
	script_name=$0
	sub_dir=${sub_path}/${script_name%.*}_sub.sh
	echo ${trimG} -q 20 -stringency 5 --fastqc_args \"--nogroup --outdir ${out_path}/${pool}\" -o ${out_trim_path}/${pool} --paired ${fastq_path}/${pool}/${r1} ${fastq_path}/${pool}/${r2} >> ${out_path}/${pool}/trace.cmd
	echo ${trimG} -q 20 -stringency 5 --fastqc_args \"--nogroup --outdir ${out_path}/${pool}\" -o ${out_trim_path}/${pool} --paired ${fastq_path}/${pool}/${r1} ${fastq_path}/${pool}/${r2} >> ${sub_dir}
	done
	cd ../
fi
done
