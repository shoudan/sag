wdir=`pwd`
pdir=`dirname $wdir`

if [ ! -e logs ]; then
mkdir logs
fi

while read lib; do
/global/projectb/sandbox/rqc/sliang/sag/python_code/contamContigLax.py ${lib}-ks.txt 
done < ${pdir}/libs_all
