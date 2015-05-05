allCtg=all_contigs.fasta
if test -e $allCtg
then
rm $allCtg
fi

while read lib
do
cat ${lib}.fasta >> $allCtg
done < libs_all

module load bbtools

bbmapskimmer.sh ref=$allCtg path=combined_fasta

if ! test -e mapping
then
mkdir mapping
fi

while read lib
do
bbmapskimmer.sh path=combined_fasta in=${lib}.fastq.gz basecov=mapping/${lib}.strand#.txt ambiguous=all strandedcov startcov concisecov sssr=0.98 minid=0.98 trimq=10 mintl=140 maxindel=0 strictmaxindel=t fast int=f usejni ow
done < libs_all

dedupe.sh in=$allCtg out=null maxedits=1 minoverlap=900 removecycles=f fixmultijoins=t processclusters=t findoverlap=t absorbcontainment=f numbergraphnodes=f absorbmatch=f outdot=mapping/graph.dot ow

all_libs=''
while read lib
do
    all_libs=${all_libs},$lib
done < libs_all
all_libs=${all_libs#,}

#qsub -hold_jid $jobid -cwd -notify -P gentech-rqc.p -N ${jobname}_snk_src -pe pe_slots 2 -l h_rt=12:00:00 -l ram.c=50G -o logs -e logs "/projectb/sandbox/rqc/sliang/sag/bbsplit/findContamContigs.py $all_libs > sink_source_contam_libs.txt"

cd mapping
/projectb/sandbox/rqc/sliang/sag/bbsplit/findContamContigs.py $all_libs > sink_source_contam_libs.txt
