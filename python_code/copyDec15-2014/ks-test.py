#!/usr/bin/env python


import sys
import string
from scipy import stats

def read_profile(fn):
    cname=''
    contig_name=[]
    profile=[]
    pos=[]
    with open(fn) as f:
        for line in f:
            line=line.strip()
            if line[0] == '#':
                if cname:
                    contig_name.append(cname)
                    profile.append(pos)
                cname=line[1:]
                pos=[]
            else:
                t2=map(int, line.split())
                pos.append(t2)
        contig_name.append(cname)
        profile.append(pos)

    #for i in range(len(profile)):
        #print  'cprofile[i] ', profile[i]

    return (profile, contig_name)
    

def get_profile(root):
    fn=root+'.strand1.txt'
    profile_fw, contig_name1=read_profile(fn)
    fn=root+'.strand2.txt'
    profile_rv, contig_name2=read_profile(fn)
# check consistency of contig_name
    assert len(contig_name1) == len(contig_name2), 'contig_name1,2 does not match'
    for n1, n2 in zip(contig_name1, contig_name2):
        assert n1 == n2, 'contig_name1,2 does not match {0} : {1}'.format(n1, n2)
    return(profile_fw, profile_rv, contig_name1)


def to_profile(cprofile):
    profile=[]
    for i in range(len(cprofile)):
        pos=[]
        #print  'cprofile ', cprofile
        for j in range(len(cprofile[i])):
            #print  'cprofile[i] ', cprofile[i]
            for r in range(cprofile[i][j][1]):
                pos.append(cprofile[i][j][0])
        profile.append(pos)
    return profile

#compute to flat distribution
def ks_flat(fw, rv):
# make a flat distribution
    flat=range(fw[0],fw[-1])
    #print fw[0], fw[-1], len(fw)
    try:
        f=stats.ks_2samp(fw, flat)[1]
    except RuntimeWarning:
        f=1

    flat=range(rv[0],rv[-1])
    #print rv[0], rv[-1], len(rv)
    try:
        r=stats.ks_2samp(rv, flat)[1]
    except RuntimeWarning:
        r=1
    return (f,r)

# first arg is the library being decontaminated
# the second arg is a list of libs in the plate
lib_decontam=sys.argv[1]
libs_in_plate=string.split(sys.argv[2], ',')
plate_fw=[]
plate_rv=[]
plate_cn=[]
for lib in libs_in_plate:
    #print "reading %s" % lib
    fn=lib+'.'+lib_decontam
    profile_fw, profile_rv, contigname=get_profile(fn)
    #print len(contigname), contigname
    plate_fw.append(profile_fw)
    plate_rv.append(profile_rv)
    plate_cn.append(contigname)
    #print lib, len(profile_fw), len(profile_rv), len(contigname)
    if lib == lib_decontam:
        decontam_fw=to_profile(profile_fw)
        decontam_rv=to_profile(profile_rv)
        decontam_cn=contigname


# check data, len is the number of contigs in lib_decontam
dlen=len(decontam_fw)
assert len(decontam_rv) == dlen, 'rv bad'
for lib, pltcn in zip(libs_in_plate, plate_cn):
    for dcn, pcn in zip (decontam_cn, pltcn):
        assert dcn == pcn, 'lib: {0} ; dcontig: {1} ; contig: {2}'.format(lib, dcn, pcn)
        #print 'lib: {0} ; dcontig: {1} ; contig: {2}'.format(lib, dcn, pcn)

flat_pvdfw={}
flat_pvdrv={}
for dfw, drv, cn in zip(decontam_fw, decontam_rv, decontam_cn):
    if len(dfw) <= 1 or len(drv) <= 1:
        flat_pvdfw[cn]=flat_pvdrv[cn]=1.
    else:
        flat_pvdfw[cn], flat_pvdrv[cn]=ks_flat(dfw,drv)
for lib, clib_fw, clib_rv in zip(libs_in_plate, plate_fw, plate_rv):
    if lib == lib_decontam:
        continue 
    assert len(decontam_fw) == len(clib_fw),  'lib_fw wrong length, decontam %d, clib %d' % (len(decontam_fw), len(clib_fw))
    assert len(decontam_rv) == len(clib_rv),  'lib_rv wrong length, decontam %d, clib %d' % (len(decontam_rv), len(clib_rv))
    lib_fw=to_profile(clib_fw)
    lib_rv=to_profile(clib_rv)
#    assert len(decontam_fw) = len(clib_fw),  'lib_fw wrong length, decontam %d, clib %d' % (len(decontam_fw), len(clib_fw))
#    assert len(decontam_rv) == len(clib_rv),  'lib_rv wrong length, decontam %d, clib %d' % (len(decontam_rv), len(clib_rv))
#    print len(decontam_fw), len(clib_fw)
#    lib_fw=to_profile(clib_fw[:len(decontam_fw)])
#    lib_rv=to_profile(clib_rv[:len(decontam_rv)])
    #print decontam_cn
    for dfw, drv, fw, rv, cn in zip(decontam_fw, decontam_rv, lib_fw, lib_rv, decontam_cn):
        if (len(dfw) == 0 and len(drv) == 0) or (len(fw) == 0 and len(rv) == 0):
            continue
        elif len(dfw) <= 1 or len(drv) <= 1 or len(fw) <= 1 or len(rv) <= 1:
            pvfw=pvrv=1.
            flat_pvfw=flat_pvrv=1.
        else:
            try:
                pvfw=stats.ks_2samp(fw, dfw)[1]
            except RuntimeWarning:
                pvfw=1.
            try:
                pvrv=stats.ks_2samp(rv, drv)[1]
            except RuntimeWarning:
                pvrv=1.
            #flat_pvdfw, flat_pvdrv=ks_flat(dfw,drv)
            flat_pvfw, flat_pvrv=ks_flat(fw,rv)
        print '\t'.join([lib, str(len(fw)),  str(len(rv)),  "{:.2g}".format(flat_pvfw),  "{:.2g}".format(flat_pvrv),\
                lib_decontam, str(len(dfw)), str(len(drv)), "{:.2g}".format(flat_pvdfw[cn]), "{:.2g}".format(flat_pvdrv[cn]), \
                cn, "{:.2g}".format(pvfw), "{:.2g}".format(pvrv)])
        sys.stdout.flush()

