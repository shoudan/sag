wdir=`pwd`
pdir=`dirname $wdir`

while read lib; do
/global/homes/c/copeland/local/x86_64/bin/faUnGetList $pdir/${lib}.fasta $pdir/ks_lax/${lib}-contig.txt ${lib}-clean.fasta
done < ${pdir}/libs_all
