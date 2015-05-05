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
import math
from collections import defaultdict

def find_source_and_sink(contigs_by_lib_name, lst_contig_overlap, profile_fw, profile_rv):
    # all following dict are indexed by tuple of (source_library, sink_library)
    n_contained_contigs=defaultdict(int)  # how many contigs are contained in this library
    nbp_contained_contigs=defaultdict(int)  # how many contigs basepairs are contained in this library
    snkRead_to_snkCtg=defaultdict(int)      # sum of sink reads mapped to sink contigs
    srcRead_to_snkCtg=defaultdict(int)      # sum of source reads mapped to sink contigs
    snkRead_to_srcCtg=defaultdict(int)      # sum of sink reads mapped to source contigs containing the sink contig
    srcRead_to_srcCtg=defaultdict(int)      # sum of source reads mapped to source contigs containing the sink contig
    aveR_to_snkCtg=defaultdict(float)       # average of contam ratio; for variance
    aveRN_to_snkCtg=defaultdict(float)      # average of contam ratio weighted by the  of reanumberd; for variance
    aveR2N_to_snkCtg=defaultdict(float)     # average of contam ratio square weighted; for variance
    aveN_to_snkCtg=defaultdict(float)       # average of the  of reanumberd; for variance
    ncnt_to_snkCtg=defaultdict(int)       # count to compute average; for variance
    aveR_to_srcCtg=defaultdict(float)       # average of contam ratio; for variance
    aveRN_to_srcCtg=defaultdict(float)      # average of contam ratio weighted by the  of reanumberd; for variance
    aveR2N_to_srcCtg=defaultdict(float)     # average of contam ratio square weighted; for variance
    aveN_to_srcCtg=defaultdict(float)       # average of the  of reanumberd; for variance
    ncnt_to_srcCtg=defaultdict(int)       # count to compute average; for variance

    contamRatio_from_snkCtg=defaultdict(float)
    variance_from_snk=defaultdict(float)
    contamRatio_from_srcCtg=defaultdict(float)
    variance_from_src=defaultdict(float)
#    for snkCtg in lst_all_contigs.itervalues():
    
    srcCtg_of_snkCtg, srcCtg_of_snkCtg_region=find_sourceContig(contigs_by_lib_name, lst_contig_overlap, profile_fw, profile_rv)
    data_maxmin, random_maxmin = monte_carlo(contigs_by_lib_name, srcCtg_of_snkCtg, srcCtg_of_snkCtg_region, profile_fw, profile_rv)
    # debugging
    for k in data_maxmin:
        print k
        print sorted(data_maxmin[k])
        print sorted(random_maxmin[k])

    print 'srcCtg_of_snkCtg'
    for snk in srcCtg_of_snkCtg:
        print snk, srcCtg_of_snkCtg[snk]
    
    for lib_snk in contigs_by_lib_name:
        for snkCtg in contigs_by_lib_name[lib_snk]:    # loop through all contigs treated as snk
            
            '''
            containing_sources = [lib_src for srcCtg in snkCtg.contained_in]
            # source_libs[lib_snk] : list of all source libraries, for lib_snk; initialize to include all libraries
            # skip any contig that is contained in more than one source
            if 1 != len(set(containing_sources) & set(source_libs[lib_snk])):
                continue
            
            srcCtg = srcCtg_by_snkCtg[snkCtg]
            if lib_src not in source_libs[lib_snk]:
                continue
            '''
            # only register these source contigs that has a sink match as recorded in srcCtg_of_snkCtg
            if snkCtg not in srcCtg_of_snkCtg:
                continue
            
            srcCtg=srcCtg_of_snkCtg[snkCtg]
            lib_src=libname(srcCtg)

            n_contained_contigs[(lib_src, lib_snk)]+=1
            nbp_contained_contigs[(lib_src, lib_snk)]+=contigs_by_lib_name[lib_snk][snkCtg].ctg_len
            
            
            snk, src = num_overlapping_reads (profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)],
                                              profile_fw[(lib_src, snkCtg)], profile_rv[(lib_src, snkCtg)])
            snkRead_to_snkCtg[(lib_src, lib_snk)]+=snk
            srcRead_to_snkCtg[(lib_src, lib_snk)]+=src
            if 0 != src:
                x=float(snk)/src
                aveR_to_snkCtg[(lib_src, lib_snk)]+=x
                aveRN_to_snkCtg[(lib_src, lib_snk)] += x*snk
                aveR2N_to_snkCtg[(lib_src, lib_snk)] += x*x*snk
                aveN_to_snkCtg[(lib_src, lib_snk)] += snk
                ncnt_to_snkCtg[(lib_src, lib_snk)] += 1
            
            snk, src = num_overlapping_reads(profile_fw[(lib_snk, srcCtg)], profile_rv[(lib_snk, srcCtg)],
                                             profile_fw[(lib_src, srcCtg)], profile_rv[(lib_src, srcCtg)])
            snkRead_to_srcCtg[(lib_src, lib_snk)]+=snk
            srcRead_to_srcCtg[(lib_src, lib_snk)]+=src
            if 0 != src:
                x=float(snk)/src
                aveR_to_srcCtg[(lib_src, lib_snk)]+=x
                aveRN_to_srcCtg[(lib_src, lib_snk)] += x*snk
                aveR2N_to_srcCtg[(lib_src, lib_snk)] += x*x*snk
                aveN_to_srcCtg[(lib_src, lib_snk)] += snk
                ncnt_to_srcCtg[(lib_src, lib_snk)] += 1

    print '\t'.join(['source_lib', 'sink_lib', 'ratio_sink', 'variance_sink', 'ratio_source', 'variance_source'])
    for ctgp in n_contained_contigs:
        if 0 != ncnt_to_snkCtg[ctgp]:
            aveR_to_snkCtg[ctgp] /= ncnt_to_snkCtg[ctgp]
            variance_from_snk[ctgp] = (aveR2N_to_snkCtg[ctgp] - 2.*aveR_to_snkCtg[ctgp]*aveRN_to_snkCtg[ctgp] + aveR_to_snkCtg[ctgp]*aveR_to_snkCtg[ctgp]*aveN_to_snkCtg[ctgp])\
                                    /ncnt_to_snkCtg[ctgp]/(aveR_to_snkCtg[ctgp]*aveR_to_snkCtg[ctgp])
        else:
            variance_from_snk[ctgp] = 0.
        if 0 != ncnt_to_srcCtg[ctgp]:
            aveR_to_srcCtg[ctgp] /= ncnt_to_srcCtg[ctgp]        
            variance_from_src[ctgp] = (aveR2N_to_srcCtg[ctgp] - 2.*aveR_to_srcCtg[ctgp]*aveRN_to_srcCtg[ctgp] + aveR_to_srcCtg[ctgp]*aveR_to_srcCtg[ctgp]*aveN_to_srcCtg[ctgp])\
                                    /ncnt_to_srcCtg[ctgp]/(aveR_to_srcCtg[ctgp]*aveR_to_srcCtg[ctgp])
        else:
            variance_from_src[ctgp] = 0.

        contamRatio_from_snkCtg[ctgp]=float(snkRead_to_snkCtg[ctgp])/float(srcRead_to_snkCtg[ctgp])
        contamRatio_from_srcCtg[ctgp]=float(snkRead_to_srcCtg[ctgp])/float(srcRead_to_srcCtg[ctgp])
        err=-10.
        if snkRead_to_snkCtg[ctgp] != 0 and snkRead_to_srcCtg[ctgp]!=0:
            err=float(snkRead_to_snkCtg[ctgp]-snkRead_to_srcCtg[ctgp])/math.sqrt(min(snkRead_to_snkCtg[ctgp], snkRead_to_srcCtg[ctgp]))
        
        output=[ctgp[0], ctgp[1], contamRatio_from_snkCtg[ctgp], variance_from_snk[ctgp], contamRatio_from_srcCtg[ctgp], variance_from_src[ctgp], err]
        print '\t'.join(map(str,output))
        sys.stdout.flush()

    pvalue=prob_sink_source(data_maxmin, random_maxmin, contamRatio_from_srcCtg)
    for ctgp in pvalue:
        print pvalue

    
'''
    for ctg_in_a_lib in contigs_by_lib_name.itervalues():
        nmap_target=sorted(nmap_self,key=lambda t:t[1],reverse=True)
'''

from scipy import stats
import collections
import numpy as np
from contig import *
from kstest import *
from maxmin import *

def find_sourceContig(contigs_by_lib_name, lst_contig_overlap, profile_fw, profile_rv):
    '''
    for a contig, find its source. Among libraries with matching profile, the the source has the largest mapped reads.
    Profile match when the KS-test p-value computed between sink reads to contig and source read to contig is less than 0.005
    input: contigs_by_lib_name list all contigs by library, contigs_by_lib_name[lib] is a list of all contigs in lib
    output: (as a return value) srcCtg_of_snkCtg[sink_contig] is the name of the source_contig that is linked to the sink_contig
    '''

    linked_ctg={}
    overlap_len={}
    overlap_reg={}
    for co in lst_contig_overlap:
        # the reason for using both sink contig sn_ctg and source library src_lib as key is to make src_ctg unique
        snk_ctg = co.contig2name
        src_ctg = co.contig1name
        src_lib = libname(src_ctg)
        if (snk_ctg, src_lib) not in overlap_len or overlap_len[(snk_ctg, src_lib)] < co.overlap_len:
            linked_ctg[(snk_ctg, src_lib)]=src_ctg
            overlap_len[(snk_ctg, src_lib)]=co.overlap_len
            overlap_reg[(snk_ctg, src_lib)]=(co.start1,co.end1)
        snk_ctg = co.contig1name
        src_ctg = co.contig2name
        src_lib = libname(src_ctg)
        if (snk_ctg, src_lib) not in overlap_len or overlap_len[(snk_ctg, src_lib)] < co.overlap_len:
            linked_ctg[(snk_ctg, src_lib)]=src_ctg
            overlap_len[(snk_ctg, src_lib)]=co.overlap_len
            overlap_reg[(snk_ctg, src_lib)]=(co.start2,co.end2)

    # loop through all sink contigs, and all source lib. find the best lib (max_lib having the largest max_src_nread)
    srcCtg_of_snkCtg={}
    srcCtg_of_snkCtg_region={}
    for lib_snk in contigs_by_lib_name:
        # for a fixed lib_snk; do all source libs together
        # contained_ctg: contig names of all source libraries stored by source library names

        for snkCtg in contigs_by_lib_name[lib_snk]:
            max_src_nread=contigs_by_lib_name[lib_snk][snkCtg].nreads
            max_lib=''
            for lib_src in contigs_by_lib_name:
                if lib_snk == lib_src or (snkCtg, lib_src) not in linked_ctg:
                    continue
                nreads=num_reads(profile_fw[(lib_src, snkCtg)]) + num_reads(profile_rv[(lib_src, snkCtg)])
#                print snkCtg, lib_src, nreads, 'max', max_src_nread
                if nreads > max_src_nread:
#                    snk, src = num_overlapping_reads (profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)],
#                                                      profile_fw[(lib_src, snkCtg)], profile_rv[(lib_src, snkCtg)])
#                    if snk == 0 or src/float(snk) < max_ratio:
#                        continue
                    
                    p_ks_profile, p_ks_flat=ks_overlap_only(profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)], profile_fw[(lib_src, snkCtg)], profile_rv[(lib_src, snkCtg)])
                    print p_ks_profile, p_ks_flat
#                    if p_ks_profile > 0.005 and p_ks_flat < 0.005:
                    if p_ks_profile > 0.005:
                        max_lib=lib_src
                        max_src_nread=nreads
                        
            if max_lib: # max_lib has been reset
                if (snkCtg, max_lib) in linked_ctg:
                    srcCtg_of_snkCtg[snkCtg]=linked_ctg[(snkCtg, max_lib)]
                    srcCtg_of_snkCtg_region[snkCtg]=overlap_reg[(snkCtg, max_lib)]
                else:
                    print "Error, max_lib not insrcCtg_of_snkCtg : ", {'max_lib':max_lib}, {'snkCtg':snkCtg}, {'max_src_nread':max_src_nread}, contigs_by_lib_name[lib_snk][snkCtg].nreads
                    sys.stdout.flush()
        sys.stdout.flush()

                #  the contig from the right is the max_lib that hitted snkCtg
    return srcCtg_of_snkCtg, srcCtg_of_snkCtg_region

# given a source library lib_src and sink library lib_snk as well as a list of matching source-sink contig pairs, 
# find the minimum coverage (mc), which is defined as contam ratio times the minumum of source read to source contigs
# in region matched by sink contig.
# Given a cutoff for the mc, find how many mc are above and below the cutoff.
# determine the expected ratio by a monte carlo calculation using the same length distribution as for the matching contigs.
# E_> = r*N and E_<=(1-r)*N, N=n_>+n_<; chi square test X^2=(n_> - E_>)^2/E_> + (n_< - E_<)^2/E_<

def above_cutoff(maxmin, cutoff):
    '''
    return the number of items in the list smaller than cutoff
    '''
    return sum(1 for x in maxmin if x > cutoff)

from scipy import stats

def prob_sink_source(data_maxmin, random_maxmin, contamRatio, threshold=5.):
    '''
    assign a probability to the sink-source libraries
    data_maxmin, random_maxmin, contamRatio  are defined on the same set of (src_lib, snk_lib)
    return pvalues for each (src_lib, snk_lib) based on fisher exact test
    '''
    pval_maxmin={}
    for pair in data_maxmin:
        cutoff=threshold/contamRatio[pair]
        d1=above_cutoff(data_maxmin[pair], cutoff)
        r1=above_cutoff(random_maxmin[pair], cutoff)
        contingency_table=[[d1, len(data_maxmin[pair])-d1], [r1, len(random_maxmin[pair])-r1]]
        R, pval= stats.fisher_exact(contingency_table, 'greater')
        pval_maxmin[pair]=pval
    return pval_maxmin


# Given a list of sink contigs and its source library, ...

import numpy

def monte_carlo(contigs_by_lib_name, srcCtg_of_snkCtg, srcCtg_of_snkCtg_region, profile_fw, profile_rv, n_monte_carlo=100):
# computed from srcCtg_of_snkCtg,
# snkCtg_by_srcLib[(lib_src, lib_snk)] is a list of snkCtg that has a good hit to a source
    snkCtg_by_srcLib=collections.defaultdict(list)
    set_src_lib=set() # set of all source libraries
    for lib_snk in contigs_by_lib_name:
        for snkCtg in contigs_by_lib_name[lib_snk]:    # loop through all contigs treated as snk
            if snkCtg not in srcCtg_of_snkCtg:
                continue
            lib_src=libname(srcCtg_of_snkCtg[snkCtg])
            set_src_lib.add(lib_src)
            snkCtg_by_srcLib[(lib_src, lib_snk)].append(snkCtg)
    
# combining forward and reverse mapped profiles by source library

    profile_by_ctg=collections.defaultdict(list)
    for lib in set_src_lib:
        for ctg in contigs_by_lib_name[lib]:
            profile_by_ctg[ctg]=comb_2_profiles(profile_fw[(lib, ctg)], profile_rv[(lib, ctg)])

    random_maxmin=collections.defaultdict(list)
    data_maxmin=collections.defaultdict(list)
    start=0
    ctg=''
    for k in snkCtg_by_srcLib:
        lib_src=k[0]
        # accumulate to n_monte_carlo
        nr=int(1+n_monte_carlo/len(snkCtg_by_srcLib[k]))
        for snkCtg in snkCtg_by_srcLib[k]:
            ctg_len=contig_len(snkCtg)
            data_maxmin[k].append( maxmin_region(profile_by_ctg[srcCtg_of_snkCtg[snkCtg]], srcCtg_of_snkCtg_region[snkCtg]) )
            for _ in range(nr):
                ctg,start=pick_random_region(lib_src, ctg_len, contigs_by_lib_name)
                #print ctg, start
                random_maxmin[k].append( maxmin_region(profile_by_ctg[ctg], (start, start+ctg_len)) )
            
            # now improve statitistics by lenghening random_maxmin until it has same maximun as data_maxmin
        cutoff=numpy.median(data_maxmin[k])-0.001
        for _ in range(n_monte_carlo):
            if max(random_maxmin[k]) >= cutoff:
                break
            for snkCtg in snkCtg_by_srcLib[k]:
                ctg_len=contig_len(snkCtg)
                ctg,start=pick_random_region(lib_src, ctg_len, contigs_by_lib_name)
                random_maxmin[k].append( maxmin_region(profile_by_ctg[ctg], (start, start+ctg_len)) )

        d1=above_cutoff(data_maxmin[k], cutoff)
        r1=above_cutoff(random_maxmin[k], cutoff)
        contingency_table=[[d1, len(data_maxmin[k])-d1], [r1, len(random_maxmin[k])-r1]]
        R, pval= stats.fisher_exact(contingency_table, 'greater')
        print k
        print 'contingency_table', contingency_table
        print 'fisher_exact pvalue=', pval


    return (data_maxmin, random_maxmin)


    
# TODO: functions needed
# (1) compute ratio of above-below cutoff by monte carlo
# (2) randomly select a region of a given length from all contigs
# (3) given a region, find the minumum coverage


import random

def pick_random_region_test(contigs_by_lib_name):
    for lib in contigs_by_lib_name:
        l=random.randint(0, 10000)
        print 'all contig for ', lib, l
        print '\t'.join([ ctg for ctg in contigs_by_lib_name[lib]])
        for i in range(5):
            print 'randomly picked : ', pick_random_region(lib, l, contigs_by_lib_name)
    
    
def pick_random_region(lib_src, matched_len, contigs_by_lib_name):
    '''
    pick a region at random from all contigs of library lib_src; matched_len is the length of the picked region
    '''
    t_contig_len=0
    for srcCtg in contigs_by_lib_name[lib_src]:    # loop through all contigs in lib_src
        l=contig_len(srcCtg)
        if l>matched_len:
            t_contig_len+=l-matched_len
    rp=random.randint(0, t_contig_len) # a random position
    for srcCtg in contigs_by_lib_name[lib_src]:    # loop through all contigs in lib_src
        l=contig_len(srcCtg)
        if l<matched_len:
            continue
        if rp < l-matched_len:
            return (srcCtg, rp)
        rp-=l-matched_len
    return (srcCtg, rp)



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
            contained_ctg[lib_src].add(srcCtg.name)
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