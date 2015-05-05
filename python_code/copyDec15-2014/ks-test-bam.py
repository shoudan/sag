#!/usr/bin/env python


import pysam
import sys
import string
from scipy import stats

def get_profile(fn_bam):
    profile_fw=[]
    profile_rv=[]
    contig_name=[]
    samfile=pysam.Samfile(fn_bam, "rb")
    for i in range(samfile.nreferences):
        fw_pos = [il_read.pos for il_read in samfile.fetch(samfile.getrname(i)) if not il_read.is_reverse]
        rv_pos = [il_read.pos for il_read in samfile.fetch(samfile.getrname(i)) if il_read.is_reverse]
        profile_fw.append(fw_pos)
        profile_rv.append(rv_pos)
        contig_name.append(samfile.getrname(i))
    
    return(profile_fw, profile_rv, contig_name)

#compute to flat distribution
def ks_flat(fw, rv):
# make a flat distribution
    flat=range(fw[0],fw[-1])
    #print fw[0], fw[-1], len(fw)
    f=stats.ks_2samp(fw, flat)[1]

    flat=range(rv[0],rv[-1])
    #print rv[0], rv[-1], len(rv)
    r=stats.ks_2samp(rv, flat)[1]
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
    fn=lib+'.'+lib_decontam+'_sorted.bam'
    profile_fw, profile_rv, contigname=get_profile(fn)
    plate_fw.append(profile_fw)
    plate_rv.append(profile_rv)
    plate_cn.append(contigname)
    #print lib, len(profile_fw), len(profile_rv), len(contigname)
    if lib == lib_decontam:
        decontam_fw=profile_fw
        decontam_rv=profile_rv
        decontam_cn=contigname


# check data, len is the number of contigs in lib_decontam
dlen=len(decontam_fw)
assert len(decontam_rv) == dlen, 'rv bad'
for lib, lib_fw, lib_rv, pltcn in zip(libs_in_plate, plate_fw, plate_rv, plate_cn):
    assert len(lib_fw) == dlen, '%s bad'.format(lib_fw)
    assert len(lib_rv) == dlen, '%s bad'.format(lib_rv)
    for dcn, pcn in zip (decontam_cn, pltcn):
        assert dcn == pcn, 'lib: {0} ; dcontig: {1} ; contig: {2}'.format(lib, dcn, pcn)
        #print 'lib: {0} ; dcontig: {1} ; contig: {2}'.format(lib, dcn, pcn)

for lib, lib_fw, lib_rv in zip(libs_in_plate, plate_fw, plate_rv):
    if lib == lib_decontam:
        continue 
    for dfw, drv, fw, rv, cn in zip(decontam_fw, decontam_rv, lib_fw, lib_rv, decontam_cn):
        if len(dfw) <= 3 or len(drv) <= 3 or len(fw) <= 3 or len(rv) <= 3:
            pvfw=1.
            pvrv=1.
            flat_pvdfw=flat_pvdrv=1.
            flat_pvfw=flat_pvrv=1.
        else:
            pvfw=stats.ks_2samp(fw, dfw)[1]
            pvrv=stats.ks_2samp(rv, drv)[1]
            flat_pvdfw, flat_pvdrv=ks_flat(dfw,drv)
            flat_pvfw, flat_pvrv=ks_flat(fw,rv)
        print '\t'.join([lib, str(len(fw)),  str(len(rv)),  "{:.2g}".format(flat_pvfw),  "{:.2g}".format(flat_pvrv),\
                lib_decontam, str(len(dfw)), str(len(drv)), "{:.2g}".format(flat_pvdfw), "{:.2g}".format(flat_pvdrv), \
                cn, "{:.2g}".format(pvfw), "{:.2g}".format(pvrv)])

