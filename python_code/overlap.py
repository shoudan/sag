#!/usr/bin/env python

def key(lib1, lib2):
    if lib1 < lib2:
        return (lib1,lib2)
    return (lib2,lib1)

def libname(cn):
    ''' take a contig name cn; and return the lib name
    '''
    return cn.split('_')[0] #first part of contig name is library name

def names_lib_contig(cn):
    fields=cn.split('_')
    return (fields[0], fields[1]) #first part of contig name is library name

    
def contig_len(contigname):
    '''
    take the contig name from spades in the form of PUCC_NODE_5_length_17015_cov_7.47429_ID_9, extract contig length
    '''
    field=contigname.split('_')
    # handle both cases with and without library name
    for i in range(2,len(field)-1):
        if field[i] == 'length':
            return int(field[i+1])


def overlap (lst_all_contigs, lst_contig_overlap):
    '''
    return the inferred genome size from the degree of overlap
    input: set of contig names, included_contigs, that are not judged as contam -> all contigs minus contigs that may be contaminated
    lst_contig_overlap is the output from dedupe which includes all overlaps including containments
    '''
    ctg_size={}  # store total assembly size of a library
    contigs_contam=set()
    for lib in lst_all_contigs:
        # the condition 0 == len(ctg.source) signal that the contig is not contaminated
        ctg_size[lib] = sum([ctg.ctg_len for ctg in lst_all_contigs[lib].itervalues() if not ctg.contam])
        contigs_contam |= set([ctg.name for ctg in lst_all_contigs[lib].itervalues() if ctg.contam])

    overlap_btwn_libs={}
    for co in lst_contig_overlap:
        lib1 = libname(co.contig1name)
        lib2 = libname(co.contig2name)
        if not co.contained and co.contig1name not in contigs_contam and co.contig2name not in contigs_contam:
            libs=key(lib1,lib2)
            if libs not in overlap_btwn_libs:
                overlap_btwn_libs[libs]=0
            overlap_btwn_libs[libs]+=co.overlap_len

    inferred_genome_size={}
    for libs in overlap_btwn_libs:
        if libs[0] not in ctg_size or libs[1] not in ctg_size:  # so we can run for small test set
            continue
        if overlap_btwn_libs[libs] != 0:
            inferred_genome_size[libs] = ctg_size[libs[0]]*ctg_size[libs[1]]/overlap_btwn_libs[libs]
#        print libs, inferred_genome_size[libs]


    while True:
        print 'init overlap_btwn_libs' , inferred_genome_size
        if remove_bad_links(inferred_genome_size):
            break
    print 'final overlap_btwn_libs' , inferred_genome_size

# remove outliers link iteratively: links unlikely to be a part of a clique will be removed
    
    org_id_by_lib = percolating_cluster(inferred_genome_size)
    return org_id_by_lib

def get_sameOrg(org_id_by_lib):
# give org_id_by_lib (org_id_by_lib[lib] indicate organism index (same index means the same organism))
# return a list of libraries in the same organism sameOrg[lib]

    sameOrg=dict()
    for lib1 in org_id_by_lib:
        sameOrg[lib1]=[lib1]
        
    for lib1 in org_id_by_lib:
        for lib2 in org_id_by_lib:
            if lib1 != lib2 and org_id_by_lib[lib1] == org_id_by_lib[lib2]:
                sameOrg[lib1].append(lib2)
                
    return sameOrg



def get_sameOrgByPairLib(org_id_by_lib):
# give org_id_by_lib that, for each lib, returns an index indicating the organism (same index means the same organism)
# return pairs of all libs that are in the same org

    sameOrgByPairLib=set()
    for lib1 in org_id_by_lib:
        for lib2 in org_id_by_lib:
            if lib1 != lib2 and org_id_by_lib[lib1] == org_id_by_lib[lib2]:
                sameOrgByPairLib.add(key(lib1,lib2))
                
    return sameOrgByPairLib


# return org_id; org_id[lib] is an integer indicating the org for lab
def percolating_cluster(overlap_btwn_libs):
    '''
    take library names from key of overlap_btwn_libs
    assign each connect libraries an id from 0 to n_cluster-1
    not ordered by the size of the cluster
    return dict org_id : org_id[lib] is the cluster id for library lib.
    '''

    # first, assign an unique id to each library
    id=0
    org_id={}
    for lib_pair in overlap_btwn_libs:
        if lib_pair[0] not in org_id:
            org_id[lib_pair[0]]=id
            id += 1
        if lib_pair[1] not in org_id:
            org_id[lib_pair[1]]=id
            id += 1

    # then, find percolating cluster iteratively

    while True:
        done=True
        for lib_pair in overlap_btwn_libs:
            if org_id[lib_pair[0]]!=org_id[lib_pair[1]]:
                done=False
                if org_id[lib_pair[0]] > org_id[lib_pair[1]]:
                    org_id[lib_pair[0]] = org_id[lib_pair[1]]
                else:
                    org_id[lib_pair[1]] = org_id[lib_pair[0]]
        if done:
            break
    return org_id


def remove_bad_links(overlap_btwn_libs):
    # find the number of clusters
    org_id_by_lib=percolating_cluster(overlap_btwn_libs)
    #valid_org_id=set([v for v in org_id_by_lib.itervalues()])
    #print 'valid_org_id', valid_org_id
    libs_in_same_org={}
    for lib in org_id_by_lib:
        if org_id_by_lib[lib] not in libs_in_same_org:
            libs_in_same_org[org_id_by_lib[lib]]=set()
        libs_in_same_org[org_id_by_lib[lib]].add(lib)

    links_same_org={}
    for lib_pair in overlap_btwn_libs:
        if org_id_by_lib[lib_pair[0]] not in links_same_org:
            links_same_org[org_id_by_lib[lib_pair[0]]] = set()
        links_same_org[org_id_by_lib[lib_pair[0]]].add(lib_pair)

    #print 'links_same_org', links_same_org

    import numpy as np
    n_neib={}
    for lib_pair in overlap_btwn_libs:
        if lib_pair[0] not in n_neib:
            n_neib[lib_pair[0]] = 1
        else:
            n_neib[lib_pair[0]] += 1
        if lib_pair[1] not in n_neib:
            n_neib[lib_pair[1]] = 1
        else:
            n_neib[lib_pair[1]] += 1

    print 'n_neib', n_neib
    
    # this part is counting the triangle in order to accertain the graph is a clique graph
    done=True
    for org in libs_in_same_org:
        max_n_neib=max(n_neib[lib] for lib in libs_in_same_org[org])
        min_common_neib=max_n_neib
        sum_link=0
        for lib_pair in links_same_org[org]:
            sum_link += overlap_btwn_libs[lib_pair]
            lib1=lib_pair[0]
            lib2=lib_pair[1]
            neib_lib1=set()
            neib_lib2=set()
            for lp in links_same_org[org]:
                if lp[0] == lib1 or lp[1] == lib1:
                    neib_lib1.add(lp[0])
                    neib_lib1.add(lp[1])
                if lp[0] == lib2 or lp[1] == lib2:
                    neib_lib2.add(lp[0])
                    neib_lib2.add(lp[1])

            if lib1 in neib_lib1:
                neib_lib1.remove(lib1)
            if lib2 in neib_lib1:
                neib_lib1.remove(lib2)
            if lib1 in neib_lib2:
                neib_lib2.remove(lib1)
            if lib2 in neib_lib2:
                neib_lib2.remove(lib2)
            in_both=neib_lib1&neib_lib2
            if len(in_both) < min_common_neib:
                min_common_neib=len(in_both)
                min_lib_pair=lib_pair
            elif len(in_both) == min_common_neib and overlap_btwn_libs[min_lib_pair] < overlap_btwn_libs[lib_pair]:
                min_lib_pair=lib_pair
                    
        sum_link -= overlap_btwn_libs[min_lib_pair]
        mean=overlap_btwn_libs[min_lib_pair]+1
        if len(links_same_org[org]) > 1:
            mean=float(sum_link)/(len(links_same_org[org])-1)

        print 'mean', mean
        z_min=100.
        print 'max_n_neib', max_n_neib
        print 'min_common_neib', min_common_neib
        
        if (max_n_neib > 1) and (min_common_neib == 0 or (min_common_neib < max_n_neib-1 and overlap_btwn_libs[min_lib_pair] > mean)): 
            done = False
            del overlap_btwn_libs[min_lib_pair]
    return done


