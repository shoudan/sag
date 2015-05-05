#!/usr/bin/env python

# take source and sink libraries; compare the log of coverage profile

import sys
import coverageProfile 
import string
import collections
import math

from overlap import *
from contig import *
from source_sink import *
from contamCtgHighRatio import contamCtgByHighRatio
from pylab import *

# combine two profiles; compute log of coverage
def get_log_cov(prof_fw, prof_rv):
    profile=comb_2_profiles(prof_fw, prof_rv)
    snk_cov = binned_coverage(profile)
    log_cov=np.zeros(len(snk_cov))
    for i in range(len(snk_cov)):
        if snk_cov[i] == 0:
            continue
        log_cov[i]=math.log(snk_cov[i])
    return log_cov
    
if len(sys.argv) == 2:
    src_lib,snk_lib=string.split(sys.argv[1], ',')
elif len(sys.argv) == 3:
    src_lib=sys.argv[1]
    snk_lib=sys.argv[2]
else:
    print "usage:"
    raise sys.exit(1)  #exit

set_of_contigs=set()
profile_fw=dict()
profile_rv=dict()

coverageProfile.get_profile((src_lib,snk_lib), set_of_contigs, profile_fw, profile_rv)

contigs_by_lib_name={}
fill_contigs_simple(set_of_contigs, profile_fw, profile_rv, contigs_by_lib_name)
for key in contigs_by_lib_name:
    if key[:5]==src_lib[:5]:
        srclib=key
    if key[:5]==snk_lib[:5]:
        snklib=key

print srclib

for lib in contigs_by_lib_name:
    print lib

ctg_list= (contigs_by_lib_name[srclib][ctg_n] for ctg_n in contigs_by_lib_name[srclib])
ctg_by_maxmin = sorted(ctg_list, key=lambda ctg: ctg.maxmin, reverse=True)
x=0;
nw=100
vertical_line=[]
src_x=[]
snk_x=[]
src_y=[]
snk_y=[]
for ctg in ctg_by_maxmin:
    profile=comb_2_profiles(profile_fw[(snklib, ctg.name)], profile_rv[(snklib, ctg.name)])
    snk_cov = binned_coverage(profile, nwindow=nw)
    profile=comb_2_profiles(profile_fw[(srclib, ctg.name)], profile_rv[(srclib, ctg.name)])
    src_cov = binned_coverage(profile, nwindow=nw)
    #print 'len src snk', len(snk_cov), len(src_cov)
    vx_src=[]
    vy_src=[]
    vx_snk=[]
    vy_snk=[]
    for i in range(len(src_cov)):
        if src_cov[i] < 10:
            continue
        x+=nw
        vx_src.append(x)
        vy_src.append(math.log10(src_cov[i]))


        if i >= len (snk_cov) or snk_cov[i] == 0:
            continue
        vx_snk.append(x)
        vy_snk.append(math.log10(snk_cov[i]))

    if len(vx_src) > 20:
        src_x+=vx_src
        src_y+=vy_src
        snk_x+=vx_snk
        snk_y+=vy_snk
        vertical_line.append(x+nw/2)

max_y=max(src_y)
min_y=min(snk_y)
print min_y, max_y
axis([0,x,min_y,max_y])
plot(src_x,src_y)
plot(snk_x,snk_y)
print len(vertical_line)
for x in vertical_line:
    plot((x,x), (min_y, max_y), 'y', linewidth=0.5)
savefig("test.png")


#for ctg_n in contigs_by_lib_name[src_lib]:

