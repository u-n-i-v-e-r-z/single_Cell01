#  Copyright (C) 2015 INRA
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

#  Date: 08.03.2016 15:05:18 CET

#Realised with samtools version 1.3

#!/bin/bash

#Toolspath
samTools="/usr/local/bioinfo/src/samtools/samtools-1.3/samtools"

# GET OPTION
while getopts "f:o:l:s:h" option
do
case $option in
    f)
        sam_path=$OPTARG
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
        echo "-f emplacement of .sam input files eg:/my/samfiles/directory"
        echo "-o where should be written bamfiles eg: /my/bamfiles/directory (created by the script if doesn't exists)"
        echo "-l log directory where error and output of qarray will be written (created by the script if doesn't exists) eg: /my/log/directory"
        echo "-s where the script that will be launched by qarray should be written eg: /my/out/directory"
        echo "-h Display this help"
        exit 1
        ;;
    h)
        echo "-s emplacement of .sam input files eg:/my/samfiles/directory"
        echo "-o where should be written bamfiles eg: /my/bamfiles/directory (created by the script if doesn't exists)"
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
cd ${sam_path}
for samFile in $(ls *.sam)
do
	runName=${samFile%%.*}
	echo ${samTools} view -Sb ${sam_path}/${pool}/${samFile} \> ${out_path}/${runName}.bam >> ${out_path}/trace
	echo ${samTools} view -Sb ${sam_path}/${pool}/${samFile} \> ${out_path}/${runName}.bam >> ${sub_dir}
done
