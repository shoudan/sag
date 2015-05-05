#!/bin/bash
#$ -cwd
#$ -notify
#$ -P gentech-rqc.p
#$ -N ks.oc1
#$ -pe pe_slots 2
#$ -l high.c
#$ -l h_rt=12:00:00
#$ -l ram.c=5G
#$ -o logs
#$ -e logs

wdir=`pwd`
pdir=`dirname $wdir`
mapping_dir=$wdir/mapping
if [ ! -e $mapping_dir ]
then
exit
fi

ks_dir=$wdir/ks
if [ ! -e $ks_dir ]
then
mkdir $ks_dir
fi



all_libs=''
while read lib; do
    all_libs=${all_libs},$lib
done < ${wdir}/libs_all
all_libs=${all_libs#,}

cd $mapping_dir
while read lib; do
/global/projectb/sandbox/rqc/sliang/sag/python_code/ks-test.py $lib $all_libs > $ks_dir/${lib}-ks.txt
done < ${wdir}/libs_all

