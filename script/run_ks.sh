
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
done < ${pdir}/libs_all
all_libs=${all_libs#,}

cd $mapping_dir
while read lib; do
qsub -cwd -b yes -j yes -now no -w e -N $lib -l high.c,h_rt=12:00:00,ram.c=5G -pe pe_slots 1 -o ks.log -js 500 "/global/projectb/sandbox/rqc/sliang/sag/python_code/ks-test.py $lib $all_libs > $TMPDIR/${lib}-ks.txt; cp $TMPDIR/${lib}-ks.txt $ks_dir; if [ ! $? -eq 0 ]; then sleep 10; cp $TMPDIR/${lib}-ks.txt $ks_dir; fi"
done < ${wdir}/libs_all

