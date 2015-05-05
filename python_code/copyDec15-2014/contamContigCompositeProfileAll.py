#!/usr/bin/env python

#test sorting

# take list of mapped reads

from scipy import stats    
import sys
import logging
#import collections
import string

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
if len(sys.argv)!=3:
    sys.stderr.write('usage: %s lib_contam comma-sep-libs-in-plate\n' % sys.argv[0])
    sys.exit(-1)
# setup log file

lib_decontam=sys.argv[1]
libs_in_plate=string.split(sys.argv[2], ',')
flog=lib_decontam+'.log'
contigfile=lib_decontam+'-contig.txt'
logging.basicConfig(filename=flog,format='%(asctime)s %(message)s',level=logging.INFO)

plate_fw=[]
plate_rv=[]
plate_cn=[]
plate_lib=[]
for lib in libs_in_plate:
    #print "reading %s" % lib
    fn=lib+'.'+lib_decontam
    profile_fw, profile_rv, contigname=get_profile(fn)
    #print len(contigname), contigname
    plate_fw.append(to_profile(profile_fw))
    plate_rv.append(to_profile(profile_rv))
    plate_cn.append(contigname)
    plate_lib.append(lib)
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

#flat_pvdfw={}
#flat_pvdrv={}
#for dfw, drv, cn in zip(decontam_fw, decontam_rv, decontam_cn):
#    if len(dfw) <= 1 or len(drv) <= 1:
#        flat_pvdfw[cn]=flat_pvdrv[cn]=1.
#    else:
#        flat_pvdfw[cn], flat_pvdrv[cn]=ks_flat(dfw,drv)
#for lib, clib_fw, clib_rv in zip(libs_in_plate, plate_fw, plate_rv):
#    if lib == lib_decontam:
#        continue 
#    assert len(decontam_fw) == len(clib_fw),  'lib_fw wrong length, decontam %d, clib %d' % (len(decontam_fw), len(clib_fw))
#    assert len(decontam_rv) == len(clib_rv),  'lib_rv wrong length, decontam %d, clib %d' % (len(decontam_rv), len(clib_rv))
#    assert len(decontam_fw) <= len(clib_fw),  'lib_fw wrong length, decontam %d, clib %d' % (len(decontam_fw), len(clib_fw))
#    assert len(decontam_rv) <= len(clib_rv),  'lib_rv wrong length, decontam %d, clib %d' % (len(decontam_rv), len(clib_rv))
#    print len(decontam_fw), len(clib_fw)
#    lib_fw=to_profile(clib_fw[:len(decontam_fw)])
#    lib_rv=to_profile(clib_rv[:len(decontam_rv)])

ilen=len(plate_lib)
jlen=len(decontam_cn)
Max_R_Contam=0.05
for j in range(jlen):
    sink_fw=decontam_fw[j]
    sink_rv=decontam_rv[j]
    source_fw=[]
    source_rv=[]
    for i in range(ilen):
        if lib_decontam == plate_lib[i]:
            continue
        source_fw.extend(plate_fw[i][j])
        source_rv.extend(plate_rv[i][j])
            
    #pvfw=-0.01        
    #if len(sink_fw) < Max_R_Contam*len(source_fw):
    pvfw=stats.ks_2samp(source_fw, sink_fw)[1]
    #pvrv=-0.01       
    #if len(sink_rv) < Max_R_Contam*len(source_rv):
    pvrv=stats.ks_2samp(source_rv, sink_rv)[1]


    #flat_pvdfw, flat_pvdrv=ks_flat(dfw,drv)
    #flat_pvfw, flat_pvrv=ks_flat(source_fw,source_rv)
    mssg='\t'.join([decontam_cn[j], lib_decontam, str(len(sink_fw)), str(len(sink_rv)), str(len(source_fw)), str(len(source_rv)),\
#                    "{:.2g}".format(flat_pvfw),  "{:.2g}".format(flat_pvrv),\
#                    "{:.2g}".format(flat_pvdfw[cn]), "{:.2g}".format(flat_pvdrv[cn]), \
                     "{:.2g}".format(pvfw), "{:.2g}".format(pvrv)])
    logging.info(mssg)

#with open(ksfile) as f:
#    for line in f:
#        fields=line.strip().split('\t')
#        #print fields
#        lib_fr=fields[0]
#        nmap_fw_fr=int(fields[1])
#        nmap_rv_fr=int(fields[2])
#        # min is more strict than max 
#        ks_flatness_fr=min(float(fields[3]),float(fields[4]))
#        lib_2=fields[5]
#        assert lib_2 == lib_this, 'column 5 must {0}'.format(lib_this)
#        nmap_fw_to=int(fields[6])
#        nmap_rv_to=int(fields[7])
#        # min is more strict than max 
#        ks_flatness_to=min(float(fields[8]),float(fields[9]))
#        contig_name=fields[10]
#        # max is more strict than min 
#        ks_profile=max(float(fields[11]),float(fields[12]))
#        #print 'flatness_fr %g _to %g' % (ks_flatness_fr, ks_flatness_to)
#        #print 'nmap _from _fw %d _rv %d' % (nmap_fw_fr, nmap_rv_fr)
#        # check fw rv is balanced
#        #stats.binom_test (nmap_fw_fr, nmap_fw_fr+nmap_rv_fr, 0.5) > Pvalue_strand_balance:
#            # insert into dict of dict
#        nmap_cross[contig_name][lib_fr]=sink_read(nmap_fw_fr+nmap_rv_fr,ks_profile,ks_flatness_fr)
#        # check if the contig is already stored
#        if  contig_name not in contig_set:
#            contig_set.add(contig_name)
#            nmap_self.append((contig_name, nmap_fw_to+nmap_rv_to, ks_flatness_to))
#    #output
#    
#
#with open(contigfile,'w') as outf:
#    for i in range(len(nmap_self)):
#        contig=nmap_self[i][0]
#        logging.info('processing contig {0}'.format(contig))
#        if contam_contig(nmap_self[i], nmap_cross[contig], contam_libs):
#            outf.write(contig)
#            outf.write('\n')
#            outf.flush()
#
