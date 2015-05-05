#!/bin/bash
#$ -cwd
#$ -notify
#$ -P gentech-rqc.p
#$ -N contam
#$ -pe pe_slots 2
#$ -l high.c
#$ -l h_rt=12:00:00
#$ -l ram.c=5G

wdir=`pwd`
pdir=`dirname $wdir`


all_libs=''
while read lib; do
    all_libs=${all_libs},$lib
done < ${pdir}/libs_all
all_libs=${all_libs#,}

/global/projectb/sandbox/rqc/sliang/sag/python_code/contamRatio.py $all_libs > plate1_contam.txt
