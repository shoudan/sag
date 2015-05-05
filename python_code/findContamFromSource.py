#!/usr/bin/env python


import sys
import coverageProfile 
import kstest 
import string
import collections
import math

from overlap import *
from contig import *
from source_sink import *
from maxmin import *
from contamCtgHighRatio import contamCtgByHighRatio


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
    


# compute partially Pearson coefficient
def partial_Pearson(weight, x, y):
    xy_ave =0.
    x_ave =0.
    x_ave_no_w=0.
    x2_ave_no_w=0.
    y_ave =0.
    y2_ave =0.
    w_ave =0.
    n=min(len(x), len(y))
    for i in range(n):
        if weight[i] == 0:
            continue
        log_y=math.log(y[i])
        xy_ave += weight[i]*x[i]*log_y
        x_ave += weight[i]*x[i]
        x_ave_no_w += x[i]
        x2_ave_no_w += x[i]*x[i]
        y_ave += weight[i]*log_y
        y2_ave += weight[i]*log_y*log_y
        w_ave += weight[i]
    return x_ave_no_w, x2_ave_no_w, x_ave, y_ave, y2_ave, xy_ave, w_ave, n

lst_of_libs=string.split(sys.argv[1], ',')
set_of_contigs=set()
profile_fw=dict()
profile_rv=dict()

coverageProfile.get_profile(lst_of_libs, set_of_contigs, profile_fw, profile_rv)

contigs_by_lib_name={}
fill_contigs_simple(set_of_contigs, profile_fw, profile_rv, contigs_by_lib_name)

#raise sys.exit(1)  #exit
src_log_cov={}
for src_lib in contigs_by_lib_name:
    ctg_list= (contigs_by_lib_name[src_lib][ctg_n] for ctg_n in contigs_by_lib_name[src_lib])
    ctg_by_maxmin = sorted(ctg_list, key=lambda ctg: ctg.maxmin, reverse=True)
    for ctg in ctg_by_maxmin:
        src_log_cov[ctg.name]=get_log_cov(profile_fw[(src_lib, ctg.name)], profile_rv[(src_lib, ctg.name)])

    for snk_lib in contigs_by_lib_name:   # loop through all contaminated sink library
        if src_lib == snk_lib:
            continue

        x_tot_no_w = 0.
        x2_tot_no_w = 0.
        x_tot = 0.
        y_tot = 0.
        y2_tot = 0.
        xy_tot = 0.
        w_tot = 0.
        n = 0.
        for ctg in ctg_by_maxmin:
            profile=comb_2_profiles(profile_fw[(snk_lib, ctg.name)], profile_rv[(snk_lib, ctg.name)])
            snk_cov = binned_coverage(profile)
                
            x_ave_no_w, x2_ave_no_w, x_ave, y_ave, y2_ave, xy_ave, w_ave, nc=partial_Pearson(snk_cov, src_log_cov[ctg.name], snk_cov)
            x_tot_no_w += x_ave_no_w
            x2_tot_no_w += x2_ave_no_w
            x_tot += x_ave
            y_tot += y_ave
            y2_tot += y2_ave
            xy_tot += xy_ave
            w_tot += w_ave
            n += nc

        if 0. == w_tot or 0 == n:
            pearson=0.
        else:
            xbar=x_tot_no_w/n
            ybar=y_tot/w_tot
            sigma_x=x2_tot_no_w/n - xbar*xbar
            sigma_y=y2_tot/w_tot - ybar*ybar
            #print sigma_x, sigma_y, y2_tot, y_tot, w_tot
            pearson = (xy_tot - ybar*x_tot - xbar*y_tot + w_tot*xbar*ybar)/w_tot/math.sqrt(sigma_x*sigma_y)
        print 'contaminant, source = ' , (snk_lib, src_lib), pearson


    

#for lib in contigs_by_lib_name:
#    for ctg in contigs_by_lib_name[lib]:    # loop through all contigs


#def xxx(contigs_by_lib_name):
#    for lib in contigs_by_lib_name:
#        for ctg in contigs_by_lib_name[lib]:    # loop through all contigs
