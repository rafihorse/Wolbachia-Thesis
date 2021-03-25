#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

__author__ = 'Serafina Nieves'
__email__ = 'smnieves@ucsc.edu'

from Bio.Phylo.PAML import codeml
import sys

wdir = str(sys.argv[1])
seqfile = str(sys.argv[2])
treefile = str(sys.argv[3])
mod = str(sys.argv[4])
outfile= str(sys.argv[5])

cml = codeml.Codeml(working_dir=wdir, alignment=seqfile, tree=treefile, out_file=outfile)

cml.set_options(noisy=9, verbose=1, runmode=0, seqtype=1, CodonFreq=2,
               ndata=0, clock=0, aaDist=0, model=mod, NSsites=[0], icode=0,
               Mgene=0, fix_kappa=0, kappa=2, fix_omega=0, omega=1,
               fix_alpha=1, alpha=0., Malpha=0, ncatG=8, getSE=0, RateAncestor=0,
               Small_Diff=.5e-6, cleandata=1, fix_blength=1, method=0)

cml.print_options()
cml.run(verbose=True)
