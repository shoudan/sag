#!/usr/bin/env python

class contig:
    def __init__ (self,cn, mm, nr):
        self.name=cn    # full name of the contig including library prefix
        fields=cn.split('_')
        self.lib=fields[0]     # library from which the contig was assembled from; Obtained from the name
        for i in range(2,len(fields)-1):
            if fields[i] == 'length':
                self.ctg_len = int(fields[i+1])
                
        self.maxmin=mm     # max min of read coverage over contig length
        self.nreads=nr     # number of reads mapped to this contig from self.lib
        self.sinks = [] # a list of sinks libraries; computed from comparing KS test on maxmin
        self.sources = [] # a list of sources libraries; computed from comparing KS profile test
        self.contained_in = [] # a list of sources computed from dedupe
        self.overlapping = [] # a list of contigs from dedupe
        self.containing = [] # a list of sinks computed from dedupe
        self.contam=False


# look over contigs
from maxmin import maxmin
from overlap import *
from coverageProfile import num_reads


def fill_contigs(lst_contig_overlap, set_all_contigs, profile_fw, profile_rv, contigs_by_lib_name):
    '''
    fill in every items in class contig except sink[], source[]
    '''
#    contigs_by_lib_name={}
    for contigname in set_all_contigs:
        lib=libname(contigname)
        #print lib, contigname
        if lib not in contigs_by_lib_name:
            contigs_by_lib_name[lib]={}
        if contigname not in contigs_by_lib_name[lib]:
            contigs_by_lib_name[lib][contigname]=contig(contigname,
                                                        max(maxmin(profile_fw[(lib, contigname)]), maxmin(profile_rv[(lib, contigname)]) ),
                                                        num_reads(profile_fw[(lib, contigname)]) + num_reads(profile_rv[(lib, contigname)]) )

    for co in lst_contig_overlap:
        lib1,contig1=libname(co.contig1name),co.contig1name
        lib2,contig2=libname(co.contig2name),co.contig2name
        # no need to check if lib1,2 and contig1,2 are in contigs_by_lib_name
        
        if lib1 not in contigs_by_lib_name or contig1 not in contigs_by_lib_name[lib1]:
            continue
        if lib2 not in contigs_by_lib_name or contig2 not in contigs_by_lib_name[lib2]:
            continue
            
        if co.contained:    # insert contained
            contigs_by_lib_name[lib2][contig2].contained_in.append(contigs_by_lib_name[lib1][contig1])
            contigs_by_lib_name[lib1][contig1].containing.append(contigs_by_lib_name[lib2][contig2])
        else:               # overlap
            #insert to 1, 2 symmetrically
            contigs_by_lib_name[lib2][contig2].overlapping.append(contigs_by_lib_name[lib1][contig1])
            contigs_by_lib_name[lib1][contig1].overlapping.append(contigs_by_lib_name[lib2][contig2])
            
            
def bymaxmin(contigs_4_lib):
    # conver dict to list; then sort in decreasing order of ctg.maxmin
    ctg_by_maxmin = sorted(list(contigs_4_lib), key=lambda ctg: ctg.maxmin, reverse=True)
