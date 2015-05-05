#!/usr/bin/env python

#test sorting

# take list of mapped reads

class sink_read:
    def __init__ (self,nmr=0,sim=0.,flat=1.):
        self.n_mapped_reads=nmr
        self.similar=sim
        self.flatness=flat

from scipy import stats

def contam_contig(nmapped_sink, nmapped_source, contam_libs, KS_THRESHOLD=0.001, P_RATIO=0.0001) :
    """
    determine if a contig is contaminated
    """
    # count # of sources
    nsource = 0
    prop=0.
    N=0
    mssg=''
    for key in contam_libs:
        mssg+=key+','
    logging.info(mssg)
    max=0 # record the library from that contributed most weight to prop
    for lib_scr in nmapped_source:
        #print lib_scr
        if lib_scr in contam_libs:
            # found a source lib
            nsource += 1
            f=float(contam_libs[lib_scr][0])/float(contam_libs[lib_scr][1])*nmapped_source[lib_scr].n_mapped_reads
            prop+=f      #proportion
            N+=nmapped_source[lib_scr].n_mapped_reads
            mssg = '%s %d' % (lib_scr, N)
            logging.info(mssg)
            if f > max:
                max=f
                lib=lib_scr
    if N!=0:
        prop /= N
    else:
        return False
    mssg='max lib %s, # of sources %d, P is %f, nmap_sink is %d, nmap_src is %d' % (lib, nsource, prop, nmapped_sink[1], N)
    logging.info(mssg)
    mssg='bin.test %g' % stats.binom_test(nmapped_sink[1], nmapped_sink[1]+N, prop/(1.+prop))
    logging.info(mssg)
    mssg='%g ~ %g, p-value %g' % (prop, nmapped_sink[1]/float(N), stats.binom_test(contam_libs[lib][0], contam_libs[lib][0]+contam_libs[lib][1], prop/(1.+prop)))
    logging.info(mssg)
    if nsource == 1:
            mssg='fisher exact %g' % stats.fisher_exact([[contam_libs[lib][0], contam_libs[lib][1]], [nmapped_sink[1], nmapped_source[lib].n_mapped_reads]])[1]
    else:
            mssg='bin.test %g' % stats.binom_test(nmapped_sink[1], nmapped_sink[1]+N, prop/(1.+prop))
    logging.info(mssg)

    if nsource > 1 and stats.binom_test(contam_libs[lib][0], contam_libs[lib][0]+contam_libs[lib][1], prop/(1.+prop)) > 0.05:
        nsource=1
    if  nsource==0 or (nsource==1 and nmapped_source[lib].similar < KS_THRESHOLD) \
    or (nsource==1 and stats.fisher_exact([[contam_libs[lib][0], contam_libs[lib][1]], [nmapped_sink[1], nmapped_source[lib].n_mapped_reads]])[1] < P_RATIO)\
    or (stats.binom_test(nmapped_sink[1], nmapped_sink[1]+N, prop/(1.+prop)) < P_RATIO):
        # not a contam
        return False

    slist=[]
    for lib_scr in nmapped_source:
        if lib_scr in contam_libs:
            slist.append(lib_scr)

    logging.info(','.join(slist))
    return True

    
import sys
import logging
import collections

if len(sys.argv)==3:
    ksfile=sys.argv[2]
    contam_lib_file=sys.argv[1]
else:
    sys.stderr.write('usage: %s contam_lib_file[all or source] ks-file\n' % sys.argv[0])
    sys.exit(-1)
# setup log file
lib_this=ksfile[:ksfile.find('-')]
flog=lib_this+'.log'
contigfile=lib_this+'-contig.txt'
logging.basicConfig(filename=flog,format='%(asctime)s %(message)s',level=logging.INFO)
logging.debug('%s' % ksfile)

Pvalue_strand_balance=0.0001
Variance_threshold=100

# determine the sources of contamination
variance=collections.defaultdict(dict) #double dict
contam_libs={}
with open(contam_lib_file) as f:
    for line in f:
        fields=line.strip().split('\t')
        assert len(fields) == 10, 'wrong file %s # of fields is %d' % (contam_lib_file, len(fields))
        lib_src=fields[0]
        lib_snk=fields[1]
        nsnk = int(fields[2])
        nsrc = int(fields[3])
        variance[lib_src][lib_snk]=float(fields[7])
        if lib_snk == lib_this and nsnk != 0 and nsrc != 0:
            contam_libs[lib_src]=(nsnk, nsrc)
            
    lib_rm=set()
    for lib in contam_libs:
        if lib in variance and lib_this in variance and lib in variance[lib_this] and lib_this in variance[lib] and\
            variance[lib_this][lib] > Variance_threshold and variance[lib][lib_this] > Variance_threshold:
            lib_rm.add(lib)
    for lib in lib_rm:
        del contam_libs[lib]

# read profile for contam libs
lib_decontam=lib_this
plate_fw=[]
plate_rv=[]
plate_cn=[]
plate_lib=[]
decontam_fw=dict()
decontam_rv=dict()
for lib in contam_libs:
    #print "reading %s" % lib
    fn=lib+'.'+lib_decontam
    profile_fw, profile_rv, contigname=get_profile(fn)
    #print len(contigname), contigname
    plate_fw.append(profile_fw)
    plate_rv.append(profile_rv)
    plate_cn.append(contigname)
    plate_lib.append(lib)
    #print lib, len(profile_fw), len(profile_rv), len(contigname)
    
assert lib_decontam not in contam_libs, '%s is its own contam lib' % lib_decontam
if lib == lib_decontam:
fn=lib_decontam+'.'+lib_decontam
profile_fw, profile_rv, contigname=get_profile(fn)
decontam_fw=to_profile(profile_fw)
decontam_rv=to_profile(profile_rv)
decontam_cn=contigname

for lib, pltcn in zip(plate_lib, plate_cn):
    assert len(decontam_cn) == len(pltcn), 'incorrect contig name, len(decontam_cn): %d; len(plate_cn): %d' % (len(decontam_cn), len(pltcn))
    for dcn, pcn in zip (decontam_cn, pltcn):
        assert dcn == pcn, 'incorrect contig name, lib: {0} ; dcontig: {1} ; contig: {2}'.format(lib, dcn, pcn)

ilen=len(plate_lib)
jlen=len(decontam_cn)
for i in range(ilen):
    r_contam[i]=contam_libs[plate_lib[i]][0]/contam_libs[plate_lib[i]][1]

for j in range(jlen):
    sink_fw=decontam_fw[j]
    sink_rv=decontam_rv[j]
    source_fw=[]
    source_rv=[]
    if (len(sink_fw) == 0 and len(sink_rv) == 0) or (len(source_fw) == 0 and len(source_rv) == 0):
        continue
    elif len(sink_fw) <= 1 or len(sink_rv) <= 1 or len(source_fw) <= 1 or len(source_rv) <= 1:
        pvfw=pvrv=1.
        flat_pvfw=flat_pvrv=1.
    else:
        for i in range(ilen):
            fw=to_profile(plate_fw[i][j])
            # check if profile from this library should be included
            nread=len(fw)
            if len(sink_fw) < Too_Small or nread*r_contam[i] > 1:
                source_fw.extend(fw)
                
            rv=to_profile(plate_rv[i][j])
            # check if profile from this library should be included
            nread=len(rv)
            if len(sink_rv) < Too_Small or nread*r_contam[i] > 1:
                source_rv.extend(rv)

        pvfw=stats.ks_2samp(source_fw, sink_fw)[1]
        pvrv=stats.ks_2samp(source_rv, sink_rv)[1]
        #flat_pvdfw, flat_pvdrv=ks_flat(dfw,drv)
        flat_pvfw, flat_pvrv=ks_flat(source_fw,source_rv)
    mssg='%s ~ %g, p-value %g' % (decontam_cn[j], nmapped_sink[1]/float(N), stats.binom_test(contam_libs[lib][0], contam_libs[lib][0]+contam_libs[lib][1], prop/(1.+prop)))
    logging.info(mssg)




#collections.defaultdict(dict) #double dict
#myhash = collections.defaultdict(lambda : collections.defaultdict(dict))  #triple dict
nmap_cross=collections.defaultdict(dict)
contig_set=set()
nmap_self=list()

with open(ksfile) as f:
    for line in f:
        fields=line.strip().split('\t')
        #print fields
        lib_fr=fields[0]
        nmap_fw_fr=int(fields[1])
        nmap_rv_fr=int(fields[2])
        # min is more strict than max 
        ks_flatness_fr=min(float(fields[3]),float(fields[4]))
        lib_2=fields[5]
        assert lib_2 == lib_this, 'column 5 must {0}'.format(lib_this)
        nmap_fw_to=int(fields[6])
        nmap_rv_to=int(fields[7])
        # min is more strict than max 
        ks_flatness_to=min(float(fields[8]),float(fields[9]))
        contig_name=fields[10]
        # max is more strict than min 
        ks_profile=max(float(fields[11]),float(fields[12]))
        #print 'flatness_fr %g _to %g' % (ks_flatness_fr, ks_flatness_to)
        #print 'nmap _from _fw %d _rv %d' % (nmap_fw_fr, nmap_rv_fr)
        # check fw rv is balanced
        #stats.binom_test (nmap_fw_fr, nmap_fw_fr+nmap_rv_fr, 0.5) > Pvalue_strand_balance:
            # insert into dict of dict
        nmap_cross[contig_name][lib_fr]=sink_read(nmap_fw_fr+nmap_rv_fr,ks_profile,ks_flatness_fr)
        # check if the contig is already stored
        if  contig_name not in contig_set:
            contig_set.add(contig_name)
            nmap_self.append((contig_name, nmap_fw_to+nmap_rv_to, ks_flatness_to))
    #output
with open(contigfile,'w') as outf:
    for i in range(len(nmap_self)):
        contig=nmap_self[i][0]
        logging.info('processing contig {0}'.format(contig))
        if contam_contig(nmap_self[i], nmap_cross[contig], contam_libs):
            outf.write(contig)
            outf.write('\n')
            outf.flush()

