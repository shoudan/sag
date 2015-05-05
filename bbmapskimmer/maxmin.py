#!/usr/bin/env python

def segment_maxmin(sprof, nwindow=100,resolution=10,len_min_contig=2000):
    '''
    compute maxmin for each segment of length, len_min_contig, of the contig.
    Maxmin is designed to simulate whether there is enough coverage to assemble
    input is the beginning positions of reads from self library mapped to the contig
    nwindow (in bp) is the window size over which to compute the average, which is computed every 'resolution' bp
    return maxmin reported as number of read coverage per bp averaged over nwindow
    '''
    r=resolution
    sum=0
    seg_cov=[]      # number of reads over window of size resolution
    for idx in range(len(sprof)):
        while sprof[idx][0] >= r:
            seg_cov.append(sum)
            sum=0
            r+=resolution
        sum += sprof[idx][1]
    seg_cov.append(sum)
#sum over segments to get averages over windows
    nseg=nwindow/resolution
    wsum=0
#    print seg_cov
    for idx in range(min(nseg, len(seg_cov))):
        wsum+=seg_cov[idx]
    window_cov=[]
    window_cov.append(wsum/nwindow)
    for idx in range(nseg, len(seg_cov)):
        wsum -= seg_cov[idx-nseg]
        wsum += seg_cov[idx]
        window_cov.append(float(wsum)/nwindow)
#    print window_cov
#compute maxmin over min contig length
    ncontig=min_contig/resolution
    contig_max=0
#    print len(window_cov), ncontig
    for idx in range(0,len(window_cov)-ncontig):
        smin=min(window_cov[idx:idx+ncontig])
#        print idx, idx+ncontig, smin
        if smin > contig_max:
            contig_max=smin
    return contig_max


def maxmin_region(prof, reg, ReadLength=150):  # ReadLength is the number of nucleotides in a read
    '''
    input is a profile in short form, return the maxmin coverage in region [from, to]
    '''
    fr=reg[0]
    to=reg[1]
    reg_prof=[r for r in prof if r[0] >= fr-ReadLength and r[0] < to]
    return maxmin(reg_prof)



def maxmin(sprof, nwindow=100,resolution=15,read_length=150,min_contig=150):
    '''
    compute maxmin over contig length. Maxmin is designed to simulate whether there is enough coverage to assemble
    input is the beginning positions of reads from self library mapped to the contig
    nwindow (in bp) is the window size over which to compute the average, which is computed every 'resolution' bp
    return maxmin reported as number of read coverage per bp averaged over nwindow
    '''
    if len(sprof)==0:
        return 0
    init=sprof[0][0]
    size = int(0.5+(sprof[-1][0]+read_length-init)/resolution)
    ilen=int(read_length/resolution)
    seg_cov=[0]*size
    for idx in range(len(sprof)):
        s=(sprof[idx][0] - init)/resolution
        for i in range(ilen):
            seg_cov[s+i]+=sprof[idx][1]
    if (len(seg_cov[ilen:-ilen]) > 0):
        return min(seg_cov[ilen:-ilen])
    return min(seg_cov)
    '''
#sum over segments to get averages over windows
    nseg=nwindow/resolution
    wsum=0
#    print seg_cov
    for idx in range(min(nseg, len(seg_cov))):
        wsum+=seg_cov[idx]
    window_cov=[]
    window_cov.append(wsum/nwindow)
    for idx in range(nseg, len(seg_cov)):
        wsum -= seg_cov[idx-nseg]
        wsum += seg_cov[idx]
        window_cov.append(wsum)
    return min(window_cov)/float(nseg)

#    print window_cov
#compute maxmin over min contig length
    ncontig=min_contig/resolution
    contig_max=0
#    print len(window_cov), ncontig
    for idx in range(0,len(window_cov)-ncontig):
        smin=min(window_cov[idx:idx+ncontig])
#        print idx, idx+ncontig, smin
        if smin > contig_max:
            contig_max=smin
    return contig_max
'''

def test_maxmin():
    prof=[]
    prof.append((1,3))
    prof.append((3,10))
    prof.append((13,10))
    prof.append((23,10))
    prof.append((73,20))
    prof.append((123,10))
    prof.append((223,20))
    prof.append((323,20))
    print maxmin(prof)
    
#test_maxmin()