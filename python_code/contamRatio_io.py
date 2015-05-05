#!/usr/bin/env python

#test sorting

# take list of mapped reads

class sink_read:
    def __init__ (self,nmr=0,sim=0.,flat=1.):
        self.n_mapped_reads=nmr
        self.similar=sim
        self.flatness=flat

from scipy import stats
# looping over the list so it only works well when the number of contigs is less than a few thousand
def contam_ratio(n50, nmapped_source, nmapped_sink, KS_THRESHOLD=0.001, KS_FLAT_THRESHOLD=0.01, P_RARIO=0.05, niter=3) :
    # n_contig_weighted_by_read: 90% of the reads are mapped to the first n50 contigs of the source
    ns=n50
    for i in range(n50):
        if nmapped_source[i][0] not in nmapped_sink:    # this means n_mapped_reads == 0
            ns=i    
            break
    
    # compute intial ratio
    nsource=0
    nsink=0        
    for i in range (ns):
        contig=nmapped_source[i][0]
        if nmapped_sink[contig].similar > KS_THRESHOLD and nmapped_sink[contig].flatness < KS_FLAT_THRESHOLD:
            nsource+=nmapped_source[i][1]
            nsink+=nmapped_sink[contig].n_mapped_reads
    print n50, ns, nsink, nsource
    if nsource + nsink == 0:  # now relax flatness 
        for i in range (ns):
            contig=nmapped_source[i][0]
            if nmapped_sink[contig].similar > KS_THRESHOLD:
                nsource+=nmapped_source[i][1]
                nsink+=nmapped_sink[contig].n_mapped_reads
    if nsource + nsink == 0:
        return (tot_mapped, n50, nsink, nsource, ns, ns)
    else:
        p=float(nsink)/(nsink+nsource)  # ratio for binomial test, differs from contam ratio
    for i in range (n50):
        contig=nmapped_source[i][0]
        pbi=stats.binom_test(nmapped_sink[contig].n_mapped_reads, nmapped_sink[contig].n_mapped_reads+nmapped_source[i][1], p)
        pchi2=stats.chisquare([nmapped_sink[contig].n_mapped_reads,], [p*nmapped_source[i][1],])[1]
        if nmapped_source[i][0] in nmapped_sink:
            mssg='contig: %s,\tratio is %d / %d = %g,\tks_profile=%g,\tks_flat=%g,\tp-value=%g, ch2=%g' % (contig, nmapped_sink[contig].n_mapped_reads, nmapped_source[i][1], nmapped_sink[contig].n_mapped_reads/float(nmapped_source[i][1]),\
            nmapped_sink[contig].similar, nmapped_sink[contig].flatness, \
            pbi, pchi2)
            logging.info(mssg)

    # binomial test
    for iter in range (niter):
        nsource=0
        nsink=0
        ns_removed=0
        for i in range (ns):
            contig=nmapped_source[i][0]
            if contig in nmapped_sink and stats.binom_test(nmapped_sink[contig].n_mapped_reads, nmapped_sink[contig].n_mapped_reads+nmapped_source[i][1], p) > P_RARIO:
                nsource+=nmapped_source[i][1]
                nsink+=nmapped_sink[contig].n_mapped_reads
            else:
                ns_removed+=1
        if nsource + nsink == 0:
            return (tot_mapped, n50, nsink, nsource, ns, ns_removed)
        p=float(nsink)/(nsink+nsource)
    return (tot_mapped, n50, nsink, nsource, ns, ns_removed)

def contam_contig(nmapped_source, nmapped_sink, KS_THRESHOLD=0.001, KS_FLAT_THRESHOLD=0.001, P_RARIO=0.05) :
    """
    determine if a contig is contaminated
    """
    # count # of sources
    nsource = 0
    prop=0.
    N=0
    for lib_fr in nmapped_source:
        if lib_fr in source_libs:
            # found a source lib
            nsource += 1
            prop+=source_libs[lib_fr]*nmapped_source[lib_fr].nmap      #proportion
            N+=nmapped_source[lib_fr].nmap
            slib=lib_fr
    if N!=0:
        prop /= N
    if nsource==0 or (nsource==1 and nmapped_source[slib].similar < KS_THRESHOLD):
        # not a contam
        pass
    if nsource==1:
        p_ratio=P_RARIO_STRICT
    else:
        p_ratio=P_RARIO

    if stat.binom_test(n, n+N, prop) < p_ratio:
        # not a contam
        pass
                
    #for lib_fr in nmapped_source:
        
    

# sorting nmapped_source
#sorted(l,key=lambda t:t[1],reverse=True)
import sys
import logging
import collections

if len(sys.argv)==2:
    libs_in_plate=sys.argv[1].split(',')
else:
    sys.stderr.write('usage: %s list of comma-separated library names\n' % sys.argv[0])
    sys.exit(-1)
# setup log file
logging.basicConfig(filename='{0}.log'.format(libs_in_plate[0]),level=logging.INFO)
logging.debug('%s' % libs_in_plate)

Pvalue_strand_balance=0.0001

#collections.defaultdict(dict) #double dict
#myhash = collections.defaultdict(lambda : collections.defaultdict(dict))  #triple dict
nmap_cross=collections.defaultdict(lambda : collections.defaultdict(dict))
nmap_cross_c=collections.defaultdict(lambda : collections.defaultdict(dict))
contig_set=collections.defaultdict(set)
nmap_self={}
for lib in libs_in_plate:
    fn=lib+'-ks.txt'
    with open(fn) as f:
        for line in f:
            fields=line.strip().split('\t')
            #print fields
            lib_fr=fields[0]
            nmap_fw_fr=int(fields[1])
            nmap_rv_fr=int(fields[2])
            # min is more strict than max 
            ks_flatness_fr=min(float(fields[3]),float(fields[4]))
            lib_to=fields[5]
            assert lib_to == lib, 'tolib is incorrect in {0}'.format(fn)
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
            nmap_cross[lib_to][lib_fr][contig_name]=sink_read(nmap_fw_fr+nmap_rv_fr,ks_profile,ks_flatness_fr)
            #print contig_name, len(nmap_cross[lib_to][lib_fr])
            nmap_cross_c[lib_to][contig_name][lib_fr]=sink_read(nmap_fw_fr+nmap_rv_fr,ks_profile,ks_flatness_fr)
            # check if the contig is already stored
            if  contig_name not in contig_set[lib_to]:
                contig_set[lib_to].add(contig_name)
                if lib_to not in nmap_self:
                    nmap_self[lib_to] = list()
                nmap_self[lib_to].append((contig_name, nmap_fw_to+nmap_rv_to, ks_flatness_to))
        #output
        logging.info('In file {0}'.format(fn))

    n_contig_weighted_by_read=0.9
    for lib_to in nmap_self:
        nmap_target=sorted(nmap_self[lib_to],key=lambda t:t[1],reverse=True)
        tot_mapped=0
        for x in nmap_target:
            tot_mapped += x[1]
        sum=0
        for i in range(len(nmap_target)):
            sum+=nmap_target[i][1]
            if sum>n_contig_weighted_by_read*tot_mapped:
                n50=i+1
                break
        logging.info('total mapped reads: {:d}'.format(tot_mapped))
        logging.info('The first {:d} contigs contains {:.1f}% of the reads'.format(n50, 100.*n_contig_weighted_by_read))

        for lib_fr in nmap_cross[lib_to]:
            print '%s %s' % (lib_to, lib_fr)
            logging.info('processing source library {0}, sink library {1}'.format(lib_to, lib_fr))
            print contam_ratio(n50, nmap_target, nmap_cross[lib_to][lib_fr])
