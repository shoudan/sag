#!/bin/sh
# needed  libs


wdir=`pwd`

module load bbtools; 

cd $TMPDIR
lib=$1
ref=$2
input=$wdir/${lib}.fastq.gz
output=${lib}.${ref}.strand#.txt
bbmap.sh sam=1.4 ref=$wdir/${ref}.fasta in=$input basecov=$output strandedcov startcov concisecov minid=0.98 trimq=10 saa=f fast int=f cigar=f usejni nodisk ow 

if [ ! -e $wdir/mapping ]
then
mkdir $wdir/mapping
fi

cp ${lib}.${ref}.strand{1,2}.txt $wdir/mapping
if [ ! $? -eq 0 ]; then 
sleep 10
cp ${lib}.${ref}.strand{1,2}.txt $wdir/mapping
fi
