module load bbtools

dedupe.sh in=all_contam.fa out=null maxedits=1 minoverlap=900 removecycles=f fixmultijoins=t processclusters=t findoverlap=t absorbmatch=f absorbcontainment=f numbergraphnodes=f absorbmatch=f outdot=graph.dot ow
#dedupe.sh in=/projectb/sandbox/rqc/sliang/sag/plate1/dedupe/all_lib.fasta out=null maxedits=1 fo minoverlap=500 fmj pc ac=f am=f ngn=f ac=f am=f outdot=graph.dot ow
#dedupe.sh in=all_lib.fasta out=null maxedits=1 fo minoverlap=500 fmj pc ngn=f ac=f am=f outdot=graph.dot ow
#dedupe.sh in=all_lib.fasta out=null maxedits=1 fo minoverlap=500 fmj pc ngn=f outdot=graph.dot ow
#maxsubs=10
#maxedits=10
# use original name
#ngn=f
# retain overlap
#ac=f am=f
#
# fmj pc
# output traph
# outdot=graph.dot
# require min overlap to 500
#fo
#mo=500
