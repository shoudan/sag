module load bbtools

# add a prefix to contig names ; the prefix is the library name
# libs_all lists all library names
while read lib
do
bbrename.sh in=${lib}.fasta out=${lib}t.fasta prefix=$lib addprefix
mv ${lib}t.fasta ${lib}.fasta
done < libs_all

wdir=`pwd`

# prepare lib separated by comma
/global/projectb/sandbox/rqc/sliang/sag/prod_code/make_task.sh > $wdir/run.config
# number of lines in run.config
njobs=`wc ${wdir}/run.config|awk '{print $1}'`
jobname=`basename $wdir`
jobid=`qsub -cwd -P gentech-rqc.p -N $jobname -pe pe_slots 2 -l h_rt=12:00:00 -l ram.c=5G -t 1-$njobs -o logs -e logs /global/projectb/sandbox/rqc/sliang/sag/prod_code/run.sh | cut -d' ' -f3 | cut -d'.' -f1`

if test -e all_contigs.fa
then
rm all_contigs.fa
fi

while read lib
do
cat ${lib}.fasta >>all_contigs.fa
done < libs_all

if ! test -e mapping
then
mkdir mapping
fi

dedupe.sh in=all_contigs.fa out=null maxedits=1 minoverlap=900 removecycles=f fixmultijoins=t processclusters=t findoverlap=t absorbcontainment=f numbergraphnodes=f absorbmatch=f outdot=mapping/graph.dot ow

all_libs=''
while read lib
do
    all_libs=${all_libs},$lib
done < libs_all
all_libs=${all_libs#,}

cd mapping

qsub -hold_jid $jobid -cwd -notify -P gentech-rqc.p -N ${jobname}_snk_src -pe pe_slots 2 -l h_rt=12:00:00 -l ram.c=50G -o logs -e logs "/projectb/sandbox/rqc/sliang/sag/prod_code/findContamContigs.py $all_libs > sink_source_contam_libs.txt"

cd ..

echo "cd mapping" > contam.sh
echo "/projectb/sandbox/rqc/sliang/sag/prod_code/main.py $all_libs > ../contamCtg.txt" >> contam.sh
echo "cd .." >> contam.sh
echo "while read lib; do /global/homes/c/copeland/local/x86_64/bin/faUnGetList \${lib}.fasta contamCtg.txt \${lib}-clean.fasta; /global/homes/c/copeland/local/x86_64/bin/faGetList \${lib}.fasta contamCtg.txt \${lib}-dirty.fasta; done < libs_all" >> contam.sh

chmod +x contam.sh

#qsub -hold_jid $jobid -cwd -notify -P gentech-rqc.p -N ${jobname}_contamCtg -pe pe_slots 2 -l h_rt=12:00:00 -l ram.c=50G logs -e logs "/projectb/sandbox/rqc/sliang/sag/python_code/main.py $all_libs > contamCtg.txt"
qsub -hold_jid $jobid -cwd -P gentech-rqc.p -N ${jobname}_contamCtg -pe pe_slots 2 -l h_rt=12:00:00 -l ram.c=50G -e logs ./contam.sh
#/projectb/sandbox/rqc/sliang/sag/python_code/main.py $all_libs > contamCtg.txt
