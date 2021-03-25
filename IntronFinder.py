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


gene_list = []
gene_prot = {}
with open("/scratch/smnieves/all_prot_gene_ids.txt", "r") as acc_file:
    for line in acc_file:
        prot, gene = line.strip().split('\t')
        gene_list.append(gene)
        gene_prot[gene] = prot

fn = annotation + ".gff"
db_fn = annotation + ".db"

#gff_db = gffutils.create_db(fn, dbfn=db_fn, force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True)

gff_db = gffutils.FeatureDB(db_fn, keep_order=True)


outfile_name = annotation + "_introns.txt"
with open(outfile_name, 'w') as outfile:

    acc_worked = []
    contigs = set()
    contains_introns = set()
    intron_check = {}
    for acc in gene_list:
        num_exons = 0
        try:
            for i in gff_db.children(gff_db[acc], featuretype='exon', order_by='start'):
                num_exons += 1
                acc_worked.append(acc)
                contigs.add((acc, gff_db[acc].seqid))

            if num_exons > 1:
                outfile.write("\n\nNon-overlapping exons in feature " + acc + ": " + str(num_exons) + "\n")

                last_exon = False
                intron_count = 0
                for i in gff_db.children(gff_db[acc], featuretype='exon', order_by='start'):
                    if last_exon == False:
                        last_exon = i
                    else:
                        if i.start - last_exon.end > 0:
                            contains_introns.add(acc)
                            outfile.write("\nIntron of length: " + str(i.start - last_exon.end) + "\n")
                            intron_count += 1
                            last_exon = i
                            print(intron_count)
                    outfile.write(str(i))
                intron_check[acc] = intron_count
                print(intron_check[acc])
            elif num_exons <= 1:
                intron_check[acc] = 0
                
        except(gffutils.exceptions.FeatureNotFoundError):
            continue
            

    outfile.write("\n\n" + str(len(acc_worked)) + " HGT genes were found in the GFF file and checked for introns. " + str(len(contains_introns)) + " genes with introns were found:\n")
    for item in contains_introns:
        outfile.write(item + "\n")
    outfile.write("\n\n")
    
    outfile.write("contig\t\tHGT\t\tflanking_gene\n")
    with open("flank_blast2.acc", "a") as outfile2:
        flank_set = set()
        for acc, ident in contigs:
            genes_on_contig = list(gff_db.region(region=(ident), featuretype=['gene']))
            outfile.write(ident + "\t" + acc + "\t")
            for gene in genes_on_contig:
                if gene.id not in gene_list:
                    outfile.write(gene.id + "\t")
                    for i in gff_db.children(gene, featuretype='CDS', order_by='start'):
                        for j in i.attributes["protein_id"]:
                            flank_set.add((ident, gene.id, j))
            outfile.write("\n")
        for cont, gene, prot in flank_set:
            outfile2.write(cont + "\t" + gene + "\t" + prot + "\n")
            
    with open("intron_check.txt", "a") as outfile3:
        for key, value in intron_check.items():
            outfile3.write(key + "\t" + str(value) + "\n")
