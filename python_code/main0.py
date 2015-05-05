#!/usr/bin/env python

#!/usr/bin/env python

import sys
from read_dedupe import read_dot
import coverageProfile 
import kstest 
import string
import collections

from overlap import *
from contig import *
from source_sink import compute_ks_by_contained
from contamCtgHighRatio import contamCtgByHighRatio


lst_contig_overlap=read_dot('graph.dot')
#print len(lst_contig_overlap)
lst_of_libs=string.split(sys.argv[1], ',')
set_of_contigs=set()
profile_fw=dict()
profile_rv=dict()

coverageProfile.get_profile(lst_of_libs, set_of_contigs, profile_fw, profile_rv)
contigs_by_lib_name={}
fill_contigs(lst_contig_overlap, set_of_contigs, profile_fw, profile_rv, contigs_by_lib_name)

#for lib in contigs_by_lib_name:
#    print lib, len(contigs_by_lib_name[lib])

contamCtg={}
contamCtgSrc={}
contamCtgByHighRatio(contigs_by_lib_name, profile_fw, profile_rv, contamCtg, contamCtgSrc)

#pairLibSameOrg=set()
#overlap (contigs_by_lib_name, lst_contig_overlap, pairLibSameOrg)

#print pairLibSameOrg
for lib in contamCtg:
    for ctg in contamCtg[lib]:
        print ctg



#sinks=collections.defaultdict(set)
#sources=collections.defaultdict(set)
#compute_ks_by_contained(contigs_by_lib_name, sinks, sources)

'''
for c in set_of_contigs:
    lib1=libname(c)
    for lib in lst_of_libs:
        lib2=lib.replace('_','')
        if lib1 == lib2 or\
           coverageProfile.num_reads(profile_fw[(lib1, c)]) < 20 or \
           coverageProfile.num_reads(profile_rv[(lib1, c)]) < 20 or \
           coverageProfile.num_reads(profile_fw[(lib2, c)]) < 20 or \
           coverageProfile.num_reads(profile_rv[(lib2, c)]) < 20:
            continue
        #print lib1, lib2, c, coverageProfile.num_overlapping_reads_1strand(profile_fw[(lib1, c)], profile_fw[(lib2, c)]), coverageProfile.num_reads(profile_fw[(lib1, c)]), coverageProfile.num_reads(profile_fw[(lib2, c)]), \
        #                     coverageProfile.num_overlapping_reads_1strand(profile_rv[(lib1, c)], profile_rv[(lib2, c)]), coverageProfile.num_reads(profile_rv[(lib1, c)]), coverageProfile.num_reads(profile_rv[(lib2, c)])
        print lib1, lib2, c, coverageProfile.num_overlapping_reads(profile_fw[(lib1, c)], profile_rv[(lib1, c)], profile_fw[(lib2, c)], profile_rv[(lib2, c)]),\
        kstest.ks_overlap_only(profile_fw[(lib1, c)], profile_rv[(lib1, c)], profile_fw[(lib2, c)], profile_rv[(lib2, c)])
                #kstest.ks_overlap_only_1strand(profile_fw[(lib1, c)], profile_fw[(lib2, c)]), kstest.ks_overlap_only_1strand(profile_rv[(lib1, c)], profile_rv[(lib2, c)]),\
        sys.stdout.flush()
                '''

