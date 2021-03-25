#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

__author__ = 'Serafina Nieves'
__email__ = 'smnieves@ucsc.edu'

import pandas as pd
import argparse
import numpy as np
import os
import gffutils

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--annotation_input', default="")
args = parser.parse_args()

annotation = args.annotation_input

prot_geneDict = {}
with open("/home/smnieves/bin/arth_queries.txt", "r") as acc_file:
    for line in acc_file:
        prot = line.rstrip()
        prot_geneDict[prot] = None
fn = annotation + ".gff"
db_fn = annotation + ".db"

#gff_db = gffutils.create_db(fn, dbfn=db_fn, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True)

gff_db = gffutils.FeatureDB(db_fn, keep_order=True)

with open("prot_gene_ids.txt", "w") as outfile:
    for acc in prot_geneDict.keys():
        try:
            prot_id = "cds-" + acc
            #rna_id = gff_db[prot_id]["Parent"][0]
            #gene_id = gff_db[rna_id]["Parent"][0]
            gene_id = "gene-" + gff_db[prot_id]["gene"][0]
            print(gene_id)
            prot_geneDict[acc] = gene_id
        except(gffutils.exceptions.FeatureNotFoundError):
            continue
        except(KeyError):
            gene_id = gff_db[prot_id]["Parent"][0]
            prot_geneDict[acc] = gene_id

    for prot, gene in prot_geneDict.items():
        if gene != None:
            outfile.write(prot + "\t" + gene + "\n")
