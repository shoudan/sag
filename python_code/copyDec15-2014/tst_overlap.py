#!/usr/bin/env python

def percolating_cluster(overlap_btwn_libs):
    '''
    take library names from key of overlap_btwn_libs
    assign each connect libraries an id from 0 to n_cluster-1
    not ordered by the size of the cluster
    return dict org_id : org_id[lib] is the cluster id for library lib.
    '''

    # first, find number of neighbors for each library
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
    valid_org_id=set([v for v in org_id_by_lib.itervalues()])
    print valid_org_id
    libs_in_same_org={}
    for id in valid_org_id:
        libs_in_same_org[id]=[]
    for lib in org_id_by_lib:
        libs_in_same_org[org_id_by_lib[lib]].append(lib)

    links_same_org={}
    for id in valid_org_id:
        links_same_org[id]=[]
    for lib_pair in overlap_btwn_libs:
        links_same_org[org_id_by_lib[lib_pair[0]]].append(lib_pair)

    print 'links_same_org', links_same_org

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
    
    
    done=True
    for org in libs_in_same_org:
        min_neib=min(n_neib[lib] for lib in libs_in_same_org[org])
        bad_link=[]
        v_good_link=[]
        for lib_pair in links_same_org[org]:
            if n_neib[lib_pair[0]] == min_neib or n_neib[lib_pair[1]] == min_neib:
                bad_link.append(lib_pair)
            else:
                v_good_link.append(overlap_btwn_libs[lib_pair])

        if len(v_good_link) < len(bad_link):
            continue

        print 'bad_link', bad_link
        print 'v_good_link', v_good_link
        mean=np.mean(v_good_link)
        std=np.std(v_good_link)
        print 'mean', mean
        print 'std', std
        z_min=100.
        
        
        for lib_pair in bad_link:
            z=(overlap_btwn_libs[lib_pair]-mean)/(std+1.)
            if z < z_min:
                z_min=z
        print 'z_min', z_min
        # remove bad links from overlap_btwn_libs
        if z_min > 0:
            done = False
            print 'inside bad_link', bad_link
            for lib_pair in bad_link:
                del overlap_btwn_libs[lib_pair]
    return done



import random

overlap_btwn_libs={}
for i in range(4):
    for j in range(i+1,3):
        overlap_btwn_libs[(i,j)]=random.random()
overlap_btwn_libs[(0,10)]=10

for i in range(4,7):
    for j in range(i+1,8):
        overlap_btwn_libs[(i,j)]=random.random()
overlap_btwn_libs[(4,11)]=0

while True:
    print 'overlap_btwn_libs', overlap_btwn_libs
    if remove_bad_links(overlap_btwn_libs):
        break
print 'final overlap_btwn_libs' , overlap_btwn_libs

org_id_by_lib = percolating_cluster(overlap_btwn_libs)
for lib1 in org_id_by_lib:
    for lib2 in org_id_by_lib:
        if lib1 != lib2 and key(lib1,lib2) not in 
            

clus = percolating_cluster(overlap_btwn_libs)
print clus
