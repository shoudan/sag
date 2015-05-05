#!/usr/bin/env python

# quantities that can help identify source-sink pairs, between a source S and a sink s
# number of contained contigs, and the base pairs of contained contigs
# contamination ratio from source contigs
# contamination ratio from the sink contigs
# variance from the source contigs
# variance from the sink contigs
# GSEA style KS test p-value comparing the maxmin of contained and not contained contigs in the source
# comparing the coverage profiles of source and sink, how many 
import sys

def find_source_and_sink(contigs_by_lib_name, lst_contig_overlap, profile_fw, profile_rv):
    # all following dict are indexed by tuple of (source_library, sink_library)
#    n_contained_contigs={}  # how many contigs are contained in this library
#    nbp_contained_contigs={}  # how many contigs basepairs are contained in this library
#    snkRead_to_snkCtg={}    # sum of sink reads mapped to sink contigs
#    srcRead_to_snkCtg={}    # sum of source reads mapped to sink contigs
#    snkRead_to_srcCtg={}    # sum of sink reads mapped to source contigs containing the sink contig
#    srcRead_to_srcCtg={}    # sum of source reads mapped to source contigs containing the sink contig
#    aveR_to_snkCtg={}   # average of contam ratio; for variance
#    aveRN_to_snkCtg={}   # average of contam ratio weighted by the  of reanumberd; for variance
#    aveR2N_to_snkCtg={}   # average of contam ratio square weighted; for variance
#    aveN_to_snkCtg={}   # average of the  of reanumberd; for variance
#    ncnt_to_snkCtg={}   # count to compute average; for variance
#    aveR_to_srcCtg={}   # average of contam ratio; for variance
#    aveRN_to_srcCtg={}   # average of contam ratio weighted by the  of reanumberd; for variance
#    aveR2N_to_srcCtg={}   # average of contam ratio square weighted; for variance
#    aveN_to_srcCtg={}   # average of the  of reanumberd; for variance
#    ncnt_to_srcCtg={}   # count to compute average; for variance

#    for snkCtg in lst_all_contigs.itervalues():
    
    srcCtg_of_snkCtg=find_sourceContig(contigs_by_lib_name, lst_contig_overlap, profile_fw, profile_rv)
    
    for lib in contigs_by_lib_name:
        for snkCtg in contigs_by_lib_name[lib].itervalues():    # loop through all contigs treated as snk
            
            '''
            containing_sources = [srcCtg.lib for srcCtg in snkCtg.contained_in]
            # source_libs[snkCtg.lib] : list of all source libraries, for snkCtg.lib; initialize to include all libraries
            # skip any contig that is contained in more than one source
            if 1 != len(set(containing_sources) & set(source_libs[snkCtg.lib])):
                continue
            
            srcCtg = srcCtg_by_snkCtg[snkCtg]
            if srcCtg.lib not in source_libs[snkCtg.lib]:
                continue
            '''
            # only register these source contigs that has a sink match as recorded in srcCtg_of_snkCtg
            if snkCtg not in srcCtg_of_snkCtg:
                continue
            
            srcCtg=srcCtg_of_snkCtg[snkCtg]

            n_contained_contigs[(srcCtg.lib, snkCtg.lib)]+=1
            nbp_contained_contigs[(srcCtg.lib, snkCtg.lib)]+=snkCtg.cnt_len
            
            snk, src = num_overlapping_reads (profile_fw[(snkCtg.lib, snkCtg)], profile_rv[(snkCtg.lib, snkCtg)],
                                              profile_fw[(srcCtg.lib, snkCtg)], profile_rv[(srcCtg.lib, snkCtg)])
            snkRead_to_snkCtg[(srcCtg.lib, snkCtg.lib)]+=snk
            srcRead_to_snkCtg[(srcCtg.lib, snkCtg.lib)]+=scr
            if 0 != src:
                x=float(snk)/src
                aveR_to_snkCtg[(srcCtg.lib, snkCtg.lib)]+=x
                aveRN_to_snkCtg[(srcCtg.lib, snkCtg.lib)] += x*snk
                aveR2N_to_snkCtg[(srcCtg.lib, snkCtg.lib)] += x*x*snk
                aveN_to_snkCtg[(srcCtg.lib, snkCtg.lib)] += snk
                ncnt_to_snkCtg[(srcCtg.lib, snkCtg.lib)] += 1
            
            snk, src = num_overlapping_reads(profile_fw[(snkCtg.lib, srcCtg)], profile_rv[(snkCtg.lib, srcCtg)],
                                             profile_fw[(srcCtg.lib, srcCtg)], profile_rv[(srcCtg.lib, srcCtg)])
            snkRead_to_srcCtg[(srcCtg.lib, snkCtg.lib)]+=snk
            srcRead_to_srcCtg[(srcCtg.lib, snkCtg.lib)]+=src
            if 0 != src:
                x=float(snk)/src
                aveR_to_srcCtg[(srcCtg.lib, snkCtg.lib)]+=x
                aveRN_to_srcCtg[(srcCtg.lib, snkCtg.lib)] += x*snk
                aveR2N_to_srcCtg[(srcCtg.lib, snkCtg.lib)] += x*x*snk
                aveN_to_srcCtg[(srcCtg.lib, snkCtg.lib)] += snk
                ncnt_to_srcCtg[(srcCtg.lib, snkCtg.lib)] += 1

    print '\t'.join(['source_lib', 'sink_lib', 'ratio_sink', 'variance_sink', 'ratio_source', 'variance_source'])
    for ctgp in n_contained_contigs:
        if 0 != ncnt_to_snkCtg[ctgp]:
            aveR_to_snk[ctgp] /= ncnt_to_snkCtg[ctgp]
            variance_from_snk[ctgp] = (aveR2N_to_snk[ctgp] - 2.*aveR_to_snk[ctgp]*aveRN_to_snk[ctgp] + aveR_to_snk[ctgp]*aveR_to_snk[ctgp]*aveN_to_snk[ctgp])\
                                    /ncnt_to_snkCtg[ctgp]/(aveR_to_snk[ctgp]*aveR_to_snk[ctgp])
        else:
            variance_from_snk[ctgp] = 0.
        if 0 != ncnt_to_srcCtg[ctgp]:
            aveR_to_scr[ctgp] /= ncnt_to_srcCtg[ctgp]        
            variance_from_scr[ctgp] = (aveR2N_to_scr[ctgp] - 2.*aveR_to_scr[ctgp]*aveRN_to_scr[ctgp] + aveR_to_scr[ctgp]*aveR_to_scr[ctgp]*aveN_to_scr[ctgp])\
                                    /ncnt_to_srcCtg[ctgp]/(aveR_to_scr[ctgp]*aveR_to_scr[ctgp])
        else:
            variance_from_src[ctgp] = 0.

        contamRatio_from_snkCtg[ctgp]=snkRead_to_snkCtg[ctgp]/srcRead_to_snkCtg[ctgp]
        contamRatio_from_srcCtg[ctgp]=snkRead_to_srcCtg[ctgp]/srcRead_to_srcCtg[ctgp]
        
        output=[ctgp[0], ctgp[1], contamRatio_from_snkCtg[ctgp], variance_from_snk[ctgp], contamRatio_from_srcCtg[ctgp], variance_from_src[ctgp]]
        print '\t'.join(map(str,output))
        sys.stdout.flush()

        

'''
    for ctg_in_a_lib in contigs_by_lib_name.itervalues():
        nmap_target=sorted(nmap_self,key=lambda t:t[1],reverse=True)

'''

from scipy import stats
import collections
import numpy as np
from contig import *

def find_sourceContig(contigs_by_lib_name, lst_contig_overlap, profile_fw, profile_rv):
    '''
    for a contig, find its source. Among libraries with matching profile, the the source has the largest mapped reads.
    Profile match when the KS-test p-value computed between sink reads to contig and source read to contig is less than 0.05
    and both sink and source profiles are not flat, defined as KS test to flat distribution
    input: contigs_by_lib_name list all contigs by library, contigs_by_lib_name[lib] is a list of all contigs in lib
    output: (as a return value) srcCtg_of_snkCtg[sink_contig] is the source_contig that is linked to the sink_contig
    '''

    linked_ctg={}
    for co in lst_contig_overlap:
        # the reason for using both sink contig sn_ctg and source library src_lib as key is to make src_ctg unique
        snk_ctg = co.contig1name
        src_ctg = co.contig2name
        src_lib = libname(src_ctg)
        linked_ctg[(snk_ctg, src_lib)]=src_ctg
        snk_ctg = co.contig2name
        src_ctg = co.contig1name
        src_lib = libname(src_ctg)
        linked_ctg[(snk_ctg, src_lib)]=src_ctg

    # loop through all sink contigs, and all source lib. find the best lib (max_lib having the largest max_src_nread)
    srcCtg_of_snkCtg={}
    for lib_snk in contigs_by_lib_name:
        # for a fixed lib_snk; do all source libs together
        # contained_ctg: contig names of all source libraries stored by source library names

        for snkCtg in contigs_by_lib_name[lib_snk].itervalues():
            max_src_nread=contigs_by_lib_name[lib_snk][snkCtg].nreads
            max_lib=''
            for lib_src in contigs_by_lib_name:
                if contigs_by_lib_name[lib_src][snkCtg].nreads > max_src_nread:
#                    snk, src = num_overlapping_reads (profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)],
#                                                      profile_fw[(lib_src, snkCtg)], profile_rv[(lib_src, snkCtg)])
#                    if snk == 0 or src/float(snk) < max_ratio:
#                        continue
                    
                    p_ks_profile, p_ks_flat=ks_overlap_only(profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)], profile_fw[(lib_src, snkCtg)], profile_rv[(lib_src, snkCtg)])
                    if p_ks_profile > 0.05 and p_ks_flat < 0.05:
                        max_lib=lib_src
                        max_src_nread=contigs_by_lib_name[lib_src][snkCtg].nreads
                        
            if max_lib:
                #  the contig from the right is the max_lib that hitted snkCtg
                srcCtg_of_snkCtg[snkCtg]=linked_ctg[(snkCtg, max_lib)]
    return srcCtg_of_snkCtg



def p_val_partition(threshold=1):
    '''
    determine a threshold for maxmin
    '''
    for ctg_name in contigs_by_lib_name[lib_src]:
        mm=contigs_by_lib_name[lib_src][ctg_name].maxmin
        n_above=0
        if mm > threshold:
            n_above+=1
            if ctg_name in source_ctg[(lib_snk, lib_src)]:
                c_above+=1
        n_below=len(contigs_by_lib_name[lib_src][ctg_name])
    

'''
        for srcCtg in snkCtg.contained_in:
            contained_ctg[srcCtg.lib].add(srcCtg.name)
        for lib_src in contigs_by_lib_name:
            if lib_src in contained_ctg:
                contained=[]
                not_contained=[]
                for ctg in contigs_by_lib_name[lib_src]:
                    if ctg in contained_ctg[lib_src]:
                        contained.append(contigs_by_lib_name[lib_src][ctg].maxmin)
                    else:
                        not_contained.append(contigs_by_lib_name[lib_src][ctg].maxmin)
 #               contained=[contigs_by_lib_name[lib_src][ctg].maxmin for ctg in contigs_by_lib_name[lib_src] if ctg in contained_ctg[lib_src]]
 #               not_contained=[contigs_by_lib_name[lib_src][ctg].maxmin for ctg in contigs_by_lib_name[lib_src] if ctg not in contained_ctg[lib_src]]
                ks_pvalue = stats.ks_2samp(contained, not_contained)[1]
                print lib_src, lib_snk, ks_pvalue, sum(contained)/len(contained), sum(not_contained)/len(not_contained)
                if ks_pvalue < 0.05 and np.median(contained) > np.median(not_contained):
                    sources[lib_snk] |= {lib_src}
                    sinks[lib_src] |= {lib_snk}
'''
                    
                    
def compute_ks_by_contained(contigs_by_lib_name, sinks, sources):
    # compute median of maxmin as well as ks p-value of contained maxmin
    for lib_snk in contigs_by_lib_name:
        # for a fixed lib_snk; do all source libs together
        # contained_ctg: contig names of all source libraries stored by source library names
        contained_ctg=collections.defaultdict(set)
        for snkCtg in contigs_by_lib_name[lib_snk].itervalues():
            for srcCtg in snkCtg.contained_in:
                contained_ctg[srcCtg.lib].add(srcCtg.name)
        for lib_src in contigs_by_lib_name:
            if lib_src in contained_ctg:
                contained=[]
                not_contained=[]
                for ctg in contigs_by_lib_name[lib_src]:
                    if ctg in contained_ctg[lib_src]:
                        contained.append(contigs_by_lib_name[lib_src][ctg].maxmin)
                    else:
                        not_contained.append(contigs_by_lib_name[lib_src][ctg].maxmin)
 #               contained=[contigs_by_lib_name[lib_src][ctg].maxmin for ctg in contigs_by_lib_name[lib_src] if ctg in contained_ctg[lib_src]]
 #               not_contained=[contigs_by_lib_name[lib_src][ctg].maxmin for ctg in contigs_by_lib_name[lib_src] if ctg not in contained_ctg[lib_src]]
                ks_pvalue = stats.ks_2samp(contained, not_contained)[1]
                print lib_src, lib_snk, ks_pvalue, sum(contained)/len(contained), sum(not_contained)/len(not_contained)
                if ks_pvalue < 0.05 and np.median(contained) > np.median(not_contained):
                    sources[lib_snk] |= {lib_src}
                    sinks[lib_src] |= {lib_snk}