#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

__author__ = 'Serafina Nieves'
__email__ = 'smnieves@ucsc.edu'

import sys
from Bio import Phylo

'''AnnotateTree.py is a program that takes a list of accessions and annotates branches on the tree with those accessions as foreground branches.'''
arth_list = []
with open("/home/smnieves/bin/arth_queries.txt") as prot_list:
    for line in prot_list:
        arth_list.append(line.strip())

treefile = str(sys.argv[1])
outfile = str(sys.argv[2])

newick = Phylo.read(treefile, "newick")
for clade in newick.find_clades(name=True):
    for acc in arth_list:
        try:
            if acc in clade.name:
                clade.name = clade.name + "#1"
        except(TypeError):
             continue
# for clade in newick.find_clades():
    # print(clade.name)

Phylo.write(newick, outfile, format="newick")
