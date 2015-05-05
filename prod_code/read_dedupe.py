#!/usr/bin/env python

from overlap import contig_len

contain_offset=100

class contig_overlap:
    def __init__ (self,cn1,cn2,ori,olen,s1,e1,s2,e2):
        self.contig1name=cn1
        self.contig2name=cn2
        self.strand_orientation=ori
        self.overlap_len=olen
        self.start1=s1
        self.end1=e1
        self.start2=s2
        self.end2=e2
        assert contig_len(cn1)+contain_offset >= contig_len(cn2), "cn1 {0}, cn2 {1}".format(cn1, cn2)
        self.contained = (olen+contain_offset > contig_len(cn2))


def read_dot(dotfile_name, substition_max=2, mutation_max=1):
    '''
    read the .dot file produced by dedupe; fill lst_contig_overlap with contig_overlap class
    '''
    dotfile=open(dotfile_name, 'r')
    line=dotfile.readline()
    assert '{' in line, 'not a dot file'
    lst_contig_overlap=[]
    while True:
        line=dotfile.readline()
        fields=line.strip().split()
        if 1==len(fields):
            if '}' in line:
                break
            cn1=fields[0]
        else:
            assert 4==len(fields), 'wrong format in dotfile'
            assert fields[0]==cn1, 'contig1 name inconsistent'
            cn2=fields[2]
            tmp=fields[3].split('"')[1].split(',')     #unpack
            ori=tmp[0]
            olen,sub,mut,s1,e1,s2,e2=map(int, tmp[1:])
            if sub <= substition_max and mut <= mutation_max:
                lst_contig_overlap.append(contig_overlap(cn1,cn2,ori,olen,s1,e1,s2,e2))
                
    return lst_contig_overlap

