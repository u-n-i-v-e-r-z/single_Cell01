#
# Copyright (C) 2017 LBME
# Schaak - Heurteau - Depierre
# Single Cell project RNAseq

#!/bin/bash

# Tools path
star="/usr/local/bioinfo/src/STAR/STAR-2.5.2b/bin/Linux_x86_64/STAR"

while getopts "i:o:" option
do
case $option in
    i)
        g_dir=$OPTARG
        ;;
    o)
        o_dir=$OPTARG
        ;;
esac
done


${star} --runMode genomeGenerate --runThreadN 4 --genomeDir ${o_dir}  --genomeFastaFiles ${g_dir}
