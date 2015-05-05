#!/usr/bin/env python

#test sorting

# take list of mapped reads

class sink_read:
    def __init__ (self,nmr=0,sim=0.,flat=1.):
        self.n_mapped_reads=nmr
        self.similar=sim
        self.flatness=flat

from scipy import stats

def contam_contig(nmapped_sink, nmapped_source, KS_THRESHOLD=0.0001, Contam_max=0.01) :
    """
    determine if a contig is contaminated
    """
    # count # of sources
    nsource = 0
    prop=0.
    max=0
    max=0 # record the library from that contributed most weight to prop
    for lib_scr in nmapped_source:
        #print lib_scr
        if nmapped_source[lib_scr].n_mapped_reads > max:
            max=nmapped_source[lib_scr].n_mapped_reads
            lib=lib_scr
    if max==0:
        return False
    mssg='max lib %s, nmap %d, nmap_sink is %d' % (lib, nmapped_source[lib].n_mapped_reads, nmapped_sink[1])
    logging.info(mssg)

    #print nmapped_source[lib].similar
    #print nmapped_sink[1]
    #print nmapped_source[lib].n_mapped_reads

    if  nmapped_source[lib].similar < KS_THRESHOLD or float(nmapped_sink[1])/nmapped_source[lib].n_mapped_reads > Contam_max:
        # not a contam
        return False
    return True

    
import sys
import logging
import collections

if len(sys.argv)==2:
    ksfile=sys.argv[1]
else:
    sys.stderr.write('usage: %s ks-file\n' % sys.argv[0])
    sys.exit(-1)
# setup log file
lib_this=ksfile[:ksfile.find('-')]
flog=lib_this+'.log'
contigfile=lib_this+'-contig.txt'
logging.basicConfig(filename=flog,format='%(asctime)s %(message)s',level=logging.INFO)
logging.debug('%s' % ksfile)

Variance_threshold=100


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
        if contam_contig(nmap_self[i], nmap_cross[contig]):
            outf.write(contig)
            outf.write('\n')
            outf.flush()
