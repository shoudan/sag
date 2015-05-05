wdir=`pwd`
pdir=`dirname $wdir`

if [ ! -e logs ]; then
mkdir logs
fi

while read lib; do
qsub -cwd -b yes -j yes -now no -w e -N $lib -l high.c,h_rt=12:00:00,ram.c=5G -pe pe_slots 1 -o logs -js 500 "/global/projectb/sandbox/rqc/sliang/sag/python_code/contamRatioOne.py ${lib} "
done < ${pdir}/libs_all

