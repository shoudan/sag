#!/bin/sh
# needed  libs


wdir=`pwd`

module load bbtools; 

cd $TMPDIR
lib=$1
input=$wdir/${lib}.fq.gz
#bbmap.sh sam=1.4 ref=$wdir/${ref}.fasta in=$input basecov=$output strandedcov startcov concisecov minid=0.98 trimq=10 saa=f fast int=f cigar=f usejni nodisk ow 
bbmapskimmer.sh path=$wdir/combined_fasta in=$input basecov=${lib}.strand#.txt ambiguous=all strandedcov startcov concisecov sssr=0.98 minid=0.98 trimq=10 mintl=140 maxindel=0 strictmaxindel=t fast int=f usejni ow

if [ ! -e $wdir/mapping ]
then
mkdir $wdir/mapping
fi

mv ${lib}.strand{1,2}.txt $wdir/mapping
if [ ! $? -eq 0 ]; then 
sleep 10
mv ${lib}.strand{1,2}.txt $wdir/mapping
fi
