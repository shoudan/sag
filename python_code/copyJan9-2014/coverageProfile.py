#!/usr/bin/env python

# coverageProfile.py

def read_profile(fn, lib4read, set_of_contigs, profile):
    cname=''
    pos=[]
    with open(fn) as f:
        for line in f:
            line=line.strip()
            if line[0] == '#':
                if cname:  # skip the first encounter
                    set_of_contigs |= {cname}
                    profile[(lib4read, cname)]=pos
                cname=line[1:]
                pos=[]
            else:
                #t2=map(int, line.split())
                line=line.split()
                p=int(line[0])
                c=int(line[1])
                #t2=tuple(map(int, line.split()))
                pos.append((p,c))
        set_of_contigs |= {cname}
        profile[(lib4read, cname)]=pos

    #for i in range(len(profile)):
        #print  'cprofile[i] ', profile[i]
#    print fn, len(profile)
#    print ' '.join([str(len(p)) for p in profile.itervalues()])
    

def get_profile(lst_of_libs, set_of_contigs, profile_fw, profile_rv):
    '''
    return the profile of reads mapped to contig in the form of list of (pos, n_read) in profile_fw (and profile_rv)
    profile_fw is a dict() indexed by [(name_of_lib_read_is_from, name_of_contig_read_mapped_to)] 
    '''
    for lib4read in lst_of_libs:
        for lib4ctg in lst_of_libs:
            fn=lib4read+'.'+lib4ctg+'.strand1.txt'
            lib=lib4read.replace('_','')
            read_profile(fn, lib, set_of_contigs, profile_fw)

            fn=lib4read+'.'+lib4ctg+'.strand2.txt'
            read_profile(fn, lib, set_of_contigs, profile_rv)

def num_reads(profile):
    # the 0th element of profile is genome location; the 1st element is # of reads mapped to the location
    return sum([pos_val[1] for pos_val in profile])

'''
alternative: similar
def comb_2_profiles(prof_fw, prof_rv):
    prof=prof_fw+prof_rv
    prof=sorted(prof,key=lambda t:t[0])
    return prof
'''

def comb_2_profiles(prof_fw, prof_rv):
    '''
    combines profiles from forward and reverse mapped reads into one profile
    '''
    ip=im=0
    prof=[]
    while True:
        if ip==len(prof_fw):
            while True:
                if im == len(prof_rv):
                    return prof
                prof.append(prof_rv[im])
                im+=1
        if im==len(prof_rv):
            while True:
                if ip == len(prof_fw):
                    return prof
                prof.append(prof_fw[ip])
                ip+=1
        if prof_fw[ip][0] < prof_rv[im][0]:
            prof.append(prof_fw[ip])
            ip+=1
        elif im<len(prof_rv) and prof_fw[ip][0] > prof_rv[im][0]:
            prof.append(prof_rv[im])
            im+=1
        else:
            prof.append((prof_fw[ip][0], prof_fw[ip][1]+prof_rv[im][1]))
            ip+=1
            im+=1
    return prof



def range_ptile(profile,pct=0.05):
    '''
    return the begining and the ending positions of 1-2*pct percentile
    '''
    num=num_reads(profile)
    sum=int(pct*num)
    s=0
    for pos_val in profile:
        sum-=pos_val[1]
        if sum <=0:
            s=pos_val[0]
            break
    sum=int(pct*num)
    e=profile[-1][0]
    for pos_val in reversed(profile):
        sum-=pos_val[1]
        if sum <=0:
            e=pos_val[0]
            break
    return (s,e)


def num_overlapping_reads_1strand(profile1, profile2):
    # compute the overlapping region of profile1 and profile2; then count the number of reads in the overlapping region only
    # not count the first and last 5%
    if 0 == len(profile1) or 0 == len(profile2):
        return (0,0)
    (s1,e1)=range_ptile(profile1)
    (s2,e2)=range_ptile(profile2)
    bgn=max(s1,s2)
    end=min(e1,e2)
    return (sum([pos_val[1] for pos_val in profile1 if pos_val[0] >= bgn and pos_val[0] <= end]),
            sum([pos_val[1] for pos_val in profile2 if pos_val[0] >= bgn and pos_val[0] <= end]))


def num_overlapping_reads(profile1_fw, profile1_rv, profile2_fw, profile2_rv):
    '''
    determined the overlap in forward and reverse region, compute the sum of reads, forward and reverse, in overlapping region 
    '''
    nfw1, nfw2 = num_overlapping_reads_1strand(profile1_fw, profile2_fw)
    nrv1, nrv2 = num_overlapping_reads_1strand(profile1_rv, profile2_rv)
    return (nfw1+nrv1, nfw2+nrv2)
