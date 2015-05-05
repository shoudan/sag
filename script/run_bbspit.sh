
all_libs=''
while read lib
do
    all_libs=${all_libs},${lib}.fasta
done < libs_all
all_libs=${all_libs#,}

module load bbtools

bbsplit.sh ref=$all_libs path=combined_fasta

if ! test -e mapping
then
mkdir mapping
fi

while read lib
do
bbsplit.sh path=combined_fasta in=${lib}.fq.gz out=mapping/${lib}.strand#.txt ambiguous=all ambiguous2=best
done < libs_all
