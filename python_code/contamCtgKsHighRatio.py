#!/usr/bin/env python

# combine ks and high ratio Jan, 2015

from coverageProfile import *
from overlap import *
from kstest import ks_overlap_only
from kstest import ks_overlap_only_1strand
from coverageProfile import comb_2_profiles


def find_ks_match_recur(source_group, start, lib_snk, snkCtg, profile_fw, profile_rv, profile_fw_src, profile_rv_src):
# loop through combination of sources in source_group, combine profile and test for matches to sink profile
# return comma separated list of source libs if match is found; otherwise return empty string
    for i in range(start, len(source_group)):
        # test profile match by KS
        lib_src=source_group[i]
        profile_fw_new=comb_2_profiles(profile_fw_src, profile_fw[lib_src, snkCtg])
        profile_rv_new=comb_2_profiles(profile_rv_src, profile_rv[lib_src, snkCtg])
        p_ks_profile, p_ks_flat=ks_overlap_only(profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)], profile_fw_new, profile_rv_new)
        if p_ks_profile > 0.01:
            #match found
            return lib_src, p_ks_profile, p_ks_flat
   
    for i in range(start+1, len(source_group)):
        # combine profile and then call
        lib_src=source_group[i]
        profile_fw_new=comb_2_profiles(profile_fw_src, profile_fw[lib_src, snkCtg])
        profile_rv_new=comb_2_profiles(profile_rv_src, profile_rv[lib_src, snkCtg])
        scrLibs, pv, pf=find_ks_match_recur(source_group, i, lib_snk, snkCtg, profile_fw, profile_rv, profile_fw_new, profile_rv_new)
        if scrLibs:
            return scrLibs+','+lib_src, pv, pf

    return '', 1e-80, 1.   # return false




def contamCtgByKS(contigs_by_lib_name, profile_fw, profile_rv, lst_sourceLib, threshold=10.):
    for lib_snk in lst_sourceLib:
        for snkCtg in contigs_by_lib_name[lib_snk]:    # loop through all contigs treated as snk
            # go through sources
            for source_group in lst_sourceLib[lib_snk]:
                srcLibs=find_ks_match_recur(source_group, 0, lib_snk, snkCtg, profile_fw, profile_rv, [], [])
            #print "from contamCtgByHighRatio", lib_snk, snkCtg, len(profile_fw)
#            if (lib_snk, snkCtg) not in profile_fw:
#                continue

#            max_num_reads_src=0
# criteria Brian used. Should move out to when contig_by_lib_name is created
            if contig_len(snkCtg) < 500:
                store_contam(lib_snk, snkCtg, contamCtg, contigs_by_lib_name)
                continue
            max_ratio=0.
            lib_max=lib_snk
            tot_num_reads_snk = num_reads(profile_fw[(lib_snk, snkCtg)]) +  num_reads(profile_rv[(lib_snk, snkCtg)])
# criteria Brian used. Should move out to when contig_by_lib_name is created
            if tot_num_reads_snk < 20:
                store_contam(lib_snk, snkCtg, contamCtg, contigs_by_lib_name)
                continue

            tot_num_reads_src={}
            for lib_src in contigs_by_lib_name: # just get all lib names
                if lib_src == lib_snk:
                    continue
                
#                if (lib_src, snkCtg) not in profile_fw:
#                    continue
                
#                num = num_reads(profile_fw[(lib_src, snkCtg)]) +  num_reads(profile_rv[(lib_src, snkCtg)])
#                if num > max_num_reads_src:
#                    max_num_reads_src=num
                snk, src = num_overlapping_reads (profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)],
                                                  profile_fw[(lib_src, snkCtg)], profile_rv[(lib_src, snkCtg)])
#                if 2*snk > tot_num_reads_snk and src/float(snk) > max_ratio:
                if 3*snk > tot_num_reads_snk and src/float(snk) > max_ratio:
                    max_ratio=src/float(snk)
                    max_lib=lib_src
#                    lib_max=lib_src
            #print snkCtg, tot_num_reads_snk, max_num_reads_src
                tot_num_reads_src[lib_src] = num_reads(profile_fw[(lib_src, snkCtg)]) +  num_reads(profile_rv[(lib_src, snkCtg)])


#contam1
#if threshold*tot_num_reads_snk < max_num_reads_src:
#contam2
#if contig_len(snkCtg) < 500 or contigs_by_lib_name[lib_snk][snkCtg].nreads < 20 or threshold*tot_num_reads_snk < max_num_reads_src:
#contam3
#            if contig_len(snkCtg) < 500 or contigs_by_lib_name[lib_snk][snkCtg].nreads < 20 or (max_ratio > threshold):
#contam4
#            if contig_len(snkCtg) < 500 or contigs_by_lib_name[lib_snk][snkCtg].nreads < 20 or (max_ratio > threshold and ks_overlap_only(profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)], profile_fw[(lib_max, snkCtg)], profile_rv[(lib_max, snkCtg)])[0] > 0.001):
#contam5
#            if contig_len(snkCtg) < 500 or contigs_by_lib_name[lib_snk][snkCtg].nreads < 20 or (max_ratio > threshold and ks_overlap_only(profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)], profile_fw[(lib_max, snkCtg)], profile_rv[(lib_max, snkCtg)])[0] > 0.0001):
            if max_ratio > threshold:
                store_contam(lib_snk, snkCtg, contamCtg, contigs_by_lib_name)
                contamCtgSrc[snkCtg]=key(lib_snk,max_lib)
            else:
                max_src_read=0
                for lib in tot_num_reads_src:
                    if tot_num_reads_src[lib] > max_src_read:
                        max_src_read = tot_num_reads_src[lib]
                        max_lib=lib
                if max_src_read / float(tot_num_reads_snk) > threshold:
                    store_contam(lib_snk, snkCtg, contamCtg, contigs_by_lib_name)
                    contamCtgSrc[snkCtg]=key(lib_snk,max_lib)
'''
                    else:
                        print snkCtg, ks_overlap_only_1strand(profile_fw[(lib_snk, snkCtg)], prof_fw)[0], tot_num_reads_src, max_src_read
                        for lib in tot_num_reads_src:
                            if tot_num_reads_src[lib]*100 > max_src_read:
                                print lib, tot_num_reads_src[lib]
                        print 'prof_fw_snk=',profile_fw[(lib_snk, snkCtg)]
                        print 'prof_fw_src=',prof_fw
'''

'''
            elif contig_len(snkCtg) < 500 or contigs_by_lib_name[lib_snk][snkCtg].nreads < 20 or (max_ratio > threshold and ks_overlap_only(profile_fw[(lib_snk, snkCtg)], profile_rv[(lib_snk, snkCtg)], profile_fw[(lib_max, snkCtg)], profile_rv[(lib_max, snkCtg)])[0] < 0.0001):
                print lib_snk, lib_max, snkCtg
                print profile_fw[(lib_snk, snkCtg)]
                print profile_fw[(lib_max, snkCtg)]
                print profile_rv[(lib_snk, snkCtg)]
                print profile_rv[(lib_max, snkCtg)]
'''




