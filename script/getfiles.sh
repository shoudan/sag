#!/bin/sh

module load jamo

#prefetch them
while read lib
do
#filtered fastq file ; copied from RQC page
file=`jamo fetch raw_normal library $lib | cut -f1 -d' '`
#make sure file is fetched
done < libs_all


while read lib
do
refdir=`/global/dna/projectdirs/PI/rqc/prod/jgi-rqc-pipeline/tools/run_folder.py -q -lib $lib |grep fastq.gz|grep "JIGSAW Single Cell"|cut -d, -f5` 
ln -s $refdir/spades/contigs.fasta ${lib}.fasta
done < libs_all

# check if file exist and link them to current dir
while read lib
do
#filtered fastq file ; copied from RQC page
file=`jamo fetch raw_normal library $lib | cut -f1 -d' '`
#make sure file is fetched
while ! test -e $file
do
sleep 10
done

ln -s $file ${lib}.fastq.gz

done < libs_all
