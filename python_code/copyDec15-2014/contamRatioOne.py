#!/usr/bin/env python

#test sorting

# take list of mapped reads

class sink_read:
    def __init__ (self,nmr=0,sim=0.,flat=1.):
        self.n_mapped_reads=nmr
        self.similar=sim
        self.flatness=flat

from scipy import stats
import math

# looping over the list so it only works well when the number of contigs is less than a few thousand
def contam_ratio(n50, nmapped_source, nmapped_sink, KS_THRESHOLD=0.001, KS_FLAT_THRESHOLD=0.001, P_RARIO=0.05, Pextreme=0.001, Nread_cut=10) :
    # n_contig_weighted_by_read: 90% of the reads are mapped to the first n50 contigs of the source
    ns=n50
    nave=ave=ave2=aven=average=0.
    for i in range(n50):
        contig=nmapped_source[i][0]
        if contig not in nmapped_sink or nmapped_sink[contig].n_mapped_reads <= Nread_cut:  # this means n_mapped_reads == 0
            if ns == n50:    # record first time n_mapped_reads == 0
                ns=i    
        else:
            src=nmapped_source[i][1]
            snk=nmapped_sink[contig].n_mapped_reads
            x=float(snk)/src
            average+=x
            ave += x*snk
            ave2 += x*x*snk
            aven += snk
            nave+=1


    if nave != 0:
        average /= nave
        var = (ave2-2.*average*ave+average*average*aven)/nave/(average*average)
    else:
        var=0.

    # compute intial ratio
    tsnk=0
    tsrc=0

    nsource=0
    nsink=0
    first_none_match=ns
    niter=-2
    for i in range (ns):
        contig=nmapped_source[i][0]
        src=nmapped_source[i][1]
        snk=nmapped_sink[contig].n_mapped_reads
        tsnk+=snk
        tsrc+=src
        niter+=1
        if nmapped_sink[contig].similar > KS_THRESHOLD and (nmapped_sink[contig].flatness < KS_FLAT_THRESHOLD or (var != 0 and var < 10.)):
            nsource+=src
            nsink+=snk
        if nmapped_sink[contig].similar < KS_THRESHOLD and first_none_match == ns:
            first_none_match=i
    mssg='before rm outlier, nsink: %d, nsource: %d' % (nsink, nsource)
    logging.info(mssg)

    #remove outliers
    outlier=set()
    for iter in range(niter):
        outlier_pv=Pextreme
        for i in range (ns):
            contig=nmapped_source[i][0]
            src=nmapped_source[i][1]
            snk=nmapped_sink[contig].n_mapped_reads
            if nmapped_sink[contig].similar > KS_THRESHOLD and nmapped_sink[contig].flatness < KS_FLAT_THRESHOLD and i not in outlier:
                pfe=stats.fisher_exact([[nsink, nsource], [snk, src]])[1]
                if pfe < Pextreme and pfe < outlier_pv:
                    ol_src=src
                    ol_snk=snk
                    outlier_pv=pfe
                    outlier.add(i)
        if outlier_pv >= Pextreme:
            break
        if (nsink == ol_snk and nsource == ol_src) and (var != 0 and var > 10):
            nsink -= ol_snk
            nsource -= ol_src

    #print 1, n50, ns, nsink, nsource
    #if nsource + nsink == 0:  # now relax flatness
    #    ns_removed=0
    #    for i in range (ns):
    #        contig=nmapped_source[i][0]
    #        if nmapped_sink[contig].similar > KS_THRESHOLD:
    #            nsource+=nmapped_source[i][1]
    #            nsink+=nmapped_sink[contig].n_mapped_reads
    #        else:
    #            ns_removed+=1
    #print 2, n50, ns, nsink, nsource
    mssg='first pass nsink: %d, nsource: %d' % (nsink, nsource)
    logging.info(mssg)
    if nsource + nsink == 0:
        return (nsink, nsource, var, ns, first_none_match, tsnk, tsrc)

    # fisher exact test
    # TODO: this loop can be made better by iteratively dropping worst outliers
    for i in range (ns):
        contig=nmapped_source[i][0]
        snk=nmapped_sink[contig].n_mapped_reads
        src=nmapped_source[i][1]
        pfe=stats.fisher_exact([[nsink, nsource], [snk, src]])[1]
        mssg='contig: %s,\tratio is %d / %d = %g,\tks_profile=%g,\tks_flat=%g,\tfisher_exact=%g' % (contig, snk, src, snk/float(src),\
        nmapped_sink[contig].similar, nmapped_sink[contig].flatness, pfe)
        logging.info(mssg)
    #    if (not (nmapped_sink[contig].similar > KS_THRESHOLD and nmapped_sink[contig].flatness < KS_FLAT_THRESHOLD)) and pfe > P_RARIO:
    #        nsink+=src
    #        nsource+=snk
    if nsource + nsink == 0:
        return (nsink, nsource, var, ns, first_none_match, tsnk, tsrc)
        #p=float(nsink)/(nsink+nsource)
    #print 4, n50, ns, nsink, nsource


#    for i in range (n50):
#        contig=nmapped_source[i][0]
#        if nmapped_source[i][0] in nmapped_sink:
#            snk=nmapped_sink[contig].n_mapped_reads
#            src=nmapped_source[i][1]
            #pbi=stats.binom_test(snk, snk+src, p)
            #mssg='contig: %s,\tratio is %d / %d = %g,\tks_profile=%g,\tks_flat=%g,\tbinomial=%g, fisher_exact=%g' % (contig, snk, src, snk/float(src),\
            #nmapped_sink[contig].similar, nmapped_sink[contig].flatness, pbi, pfe)
    return (nsink, nsource, var, ns, first_none_match, tsnk, tsrc)

    

# sorting nmapped_source
#sorted(l,key=lambda t:t[1],reverse=True)
import sys
import logging
import collections

if len(sys.argv)==2:
    lib_to=sys.argv[1]
    ksfile=lib_to+'-ks.txt'
    contamfile=lib_to+'-contam.txt'
    flog=lib_to+'.log'
else:
    sys.stderr.write('usage: %s ks-file\n' % sys.argv[0])
    sys.exit(-1)
# setup log file
#lib_to=ksfile[:ksfile.find('-')]
logging.basicConfig(filename=flog,format='%(asctime)s %(message)s',level=logging.INFO)
logging.debug('%s' % ksfile)

Pvalue_strand_balance=0.0001

#collections.defaultdict(dict) #double dict
#myhash = collections.defaultdict(lambda : collections.defaultdict(dict))  #triple dict
nmap_cross=collections.defaultdict(dict)
contig_set=set()
nmap_self=list()

with open(ksfile, 'r') as f:
    for line in f:
        fields=line.strip().split('\t')
        #print fields
        lib_fr=fields[0]
        nmap_fw_fr=int(fields[1])
        nmap_rv_fr=int(fields[2])
        # min is more strict than max 
        ks_flatness_fr=min(float(fields[3]),float(fields[4]))
        lib_2=fields[5]
        assert lib_2 == lib_to, 'column 5 must {0}'.format(lib_to)
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
        nmap_cross[lib_fr][contig_name]=sink_read(nmap_fw_fr+nmap_rv_fr,ks_profile,ks_flatness_fr)
        # check if the contig is already stored
        if  contig_name not in contig_set:
            contig_set.add(contig_name)
            nmap_self.append((contig_name, nmap_fw_to+nmap_rv_to, ks_flatness_to))
    #output

n_contig_weighted_by_read=0.9
nmap_target=sorted(nmap_self,key=lambda t:t[1],reverse=True)
tot_mapped=0
for x in nmap_target:
    tot_mapped += x[1]
sum=0
n50=len(nmap_target)
for i in range(len(nmap_target)):
    sum+=nmap_target[i][1]
    if sum>n_contig_weighted_by_read*tot_mapped:
        n50=i+1
        break
logging.info('total mapped reads: {:d}'.format(tot_mapped))
logging.info('The first {:d} contigs contains {:.1f}% of the reads'.format(n50, 100.*n_contig_weighted_by_read))

with open(contamfile,'w') as outf:
    for lib_fr in nmap_cross:
        logging.info('processing source library {0}, sink library {1}'.format(lib_to, lib_fr))
        nsnk, nsrc, var, ns, first_none_match, tsnk, tsrc = contam_ratio(n50, nmap_target, nmap_cross[lib_fr])
        if tsnk != 0 and tot_mapped != 0:
            outf.write('\t'.join([lib_to, lib_fr, str(nsnk), str(nsrc), str(ns), str(first_none_match), str(tot_mapped), '{:.2g}'.format(var), str(tsnk), str(tsrc)]))
            outf.write('\n')
            outf.flush()

