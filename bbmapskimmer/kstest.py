#!/usr/bin/env python

from scipy import stats
from coverageProfile import *


def ks_overlap_only_1strand(profile1, profile2):
    '''
    compare similarity of profile1 and profile2 by KS test
    return (ks_p-value, flat_ks-pvalue) where flat_ks-pvalue indicate if profile1 or profile 2 are flat
    '''
    if 0 == len(profile1) or 0 == len(profile2):
        return (1.0,1.0)

    # compute the overlapping region of profile1 and profile2; then count the number of reads in the overlapping region only
    # not count the first and last 5%
    (s1,e1)=range_ptile(profile1)
    (s2,e2)=range_ptile(profile2)
    bgn=max(s1,s2)
    end=min(e1,e2)
    ksprofil1e1=to_profile(profile1, bgn, end)
    ksprofil1e2=to_profile(profile2, bgn, end)
    if 0==len(ksprofil1e1) or 0==len(ksprofil1e1):
        fpv=1.
    else:
        flat1=range(ksprofil1e1[0],ksprofil1e1[-1])
        flat2=range(ksprofil1e2[0],ksprofil1e2[-1])
        fpv=max(stats.ks_2samp(ksprofil1e1, flat1)[1], stats.ks_2samp(ksprofil1e2, flat2)[1])

    return (stats.ks_2samp(ksprofil1e1, ksprofil1e2)[1], fpv)


def to_profile(cprofile, bgn, end):
    '''
    convert profile from bbmap to profile that stats.ks_2samp takes
    '''
    profile=[]
    for elem in cprofile:
        if elem[0] < bgn or elem[0] > end:
            pass
        #print  'cprofile[i] ', cprofile[i]
        for r in range(elem[1]):
            profile.append(elem[0])
    return profile


def ks_overlap_only(profile1_fw, profile1_rv, profile2_fw, profile2_rv):
    '''
    determined the overlap in forward and reverse region, compute the sum of reads, forward and reverse, in overlapping region 
    '''
    pvfw, flat_pvfw = ks_overlap_only_1strand(profile1_fw, profile2_fw)
    pvrv, flat_pvrv = ks_overlap_only_1strand(profile1_rv, profile2_rv)
    return (max(pvfw, pvrv), max(flat_pvfw, flat_pvrv))
