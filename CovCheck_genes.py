#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

__author__ = 'Serafina Nieves'
__email__ = 'smnieves@ucsc.edu'

import pandas as pd
import argparse
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import os
import random

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--SRA_input', required=False)
args = parser.parse_args()
# if args.SRA_input != None:
inputHGT = args.SRA_input + "_HGT.cov"
inputAll = args.SRA_input + "_sortedaln.cov"


def calc_coverage():
    try:
        # read samtools bedcov data into dataframe
        HGT = pd.read_csv(inputHGT, sep="\t",
                          names=["chr", "start", "end", "gene_id", "frame", "strand", "source", "type", "idk", "description", "basecov"])
        # calculate length of coding region and coverage over coding region
        HGT["region_len"] = HGT["end"] - HGT["start"]
        HGT["total_cov"] = HGT["basecov"]/HGT["region_len"]
        gene_list = HGT["gene_id"]

        # read in samtools bedcov data for every genome feature
        all_feat = pd.read_csv(inputAll, sep="\t", names=["chr", "start", "end", "gene_id", "frame", "strand", "source", "type", "idk", "description", "basecov"])

        # calculate length of every genome feature and coverage over every feature
        all_feat["region_len"] = all_feat["end"] - all_feat["start"]
        all_feat["total_cov"] = all_feat["basecov"]/all_feat["region_len"]

        # calculate the median coverage of exons within a certain number of basepairs of the coding region length
        def gene_median(row):
            gene_df = all_feat.loc[
                (all_feat["region_len"] < (row["region_len"] + 50)) & (
                            all_feat["region_len"] > (row["region_len"] - 50))]
            gene_df = gene_df[~gene_df["gene_id"].isin(gene_list)]
            return gene_df["total_cov"].median()

        # calculate the percentile of HGT coverage in exon coverage distribution
        def perc_geneDist(row):
            gene_df = all_feat.loc[
                (all_feat["region_len"] < (row["region_len"] + 50)) & (
                            all_feat["region_len"] > (row["region_len"] - 50))]
            gene_df = gene_df[~gene_df["gene_id"].isin(gene_list)]
            return stats.percentileofscore(gene_df["total_cov"], row["total_cov"])

        # calculate exon coverage for all coding regions in HGT samtools bedcov file
        HGT["median_gene_cov"] = HGT.apply(gene_median, axis=1)
        #print(HGT.apply(gene_median, axis=1))

        # calculate percentile of exon distribution that HGT coverage falls in
        HGT["cov_percentile"] = HGT.apply(perc_geneDist, axis=1)
        #print(HGT.apply(gene_median, axis=1))

        # print out cleaned dataframe
        with open(args.SRA_input + ".covcheck", "w") as out:
            print("Unique putative HGT events from these reads: ", len(set(HGT['gene_id'])))
            HGT[['chr', 'gene_id', 'region_len', 'total_cov', 'median_gene_cov', 'cov_percentile']].to_csv(args.SRA_input + ".covcheck", sep='\t')
    except(IndexError):
        print("No data in HGT.cov file for {0}. This may mean that current SRA data does not contain HGT events of interest. Try different SRA read set.".format(args.SRA_input))


def speciesHisto():
    # change dimensions of the entire figure
    figureHeight = 3
    figureWidth = 5

    plt.figure(figsize=(figureWidth, figureHeight))

    mainPanelWidth = 3
    mainPanelHeight = 2

    mainPanel = plt.axes([0.08, 0.2, mainPanelWidth / figureWidth, mainPanelHeight / figureHeight])

    keyPanelWidth = .5
    keyPanelHeight = 2

    keyPanel = plt.axes([0.705, 0.2, keyPanelWidth / figureWidth, keyPanelHeight / figureHeight])

    # change tick labels / distances
    mainPanel.set_xticks(np.arange(0, 110, 10))
    mainPanel.set_xticks(np.arange(0, 100, 2.5), minor=True)
    mainPanel.set_xlim(0, 100)
    mainPanel.set_yticks(np.arange(0, 120, 20))
    mainPanel.set_yticks(np.arange(0, 120, 5), minor=True)
    mainPanel.set_ylim(0, 110)
    mainPanel.set_title("HGT percentile in distribution of host exon coverage", pad=10, loc="left")
    mainPanel.set_xlabel("Percentile")
    mainPanel.set_ylabel("Number of exonic regions")

    keyPanel.set_xlim(0, 1)
    keyPanel.set_ylim(0, 15)
    keyPanel.tick_params(bottom=False, labelbottom=False,
                      left=False, labelleft=False,
                      right=False, labelright=False,
                      top=False, labeltop=False)


    speciesDict = {"SRR3458570.covcheck": "armadillidium vulgare", "ERR034187.covcheck": "acromyrmex echinatior", "SRR3458573.covcheck": "armadillidium vulgare", "ERR3437293.covcheck": "cinara cedri",
                   "SRR2912519.covcheck": "ceratina calcarata", "SRR2135644.covcheck": "drosophila ananassae", "SRR341538.covcheck": "drosophila eugracilis",
                   "DRR042089.covcheck": "formica exsecta", "SRR3680351.covcheck": "rhagoletis zephyria", "SRR7188756.covcheck": "sipha flava", "SRR621118.covcheck": "solenopsis invicta",
                   "SRR5139338.covcheck": "trichonephila clavipes", "DRR037409.covcheck": "vollenhovia emeryi", "DRR029092.covcheck": "wasmannia auropunctata", "SRR5830088.covcheck": "laodelphax striatellus"}

    textBottom = 0.21
    keyBottom = 0
    dataDict = {}
    colorpalette = {"SRR3458570.covcheck": "xkcd:purple", "ERR034187.covcheck": "xkcd:dusty teal", "SRR3458573.covcheck": "xkcd:azure", "ERR3437293.covcheck": "xkcd:orangered", "SRR2912519.covcheck": "xkcd:marine",
                    "SRR2135644.covcheck": "xkcd:blueberry", "SRR341538.covcheck": "xkcd:very pale green", "DRR042089.covcheck": "xkcd:ugly pink", "SRR3680351.covcheck": "xkcd:sage green", "SRR7188756.covcheck": "xkcd:dark mauve",
                    "SRR621118.covcheck": "xkcd:light peach", "SRR5139338.covcheck": "xkcd:yellow", "DRR037409.covcheck": "xkcd:sapphire", "DRR029092.covcheck": "xkcd:cerise", "SRR5830088.covcheck": "xkcd:black"}

    open("genomic_check.txt", "w").close()
    for file in speciesDict.keys():
        covcheck = pd.read_csv(file, sep='\t', header=0, names=['row', 'chr', 'id', 'region_len', 'total_cov', 'median_exon_cov', 'cov_percentile'])
        with open("genomic_check.txt", "a") as all_outfile:
            output_df = covcheck[["chr", "id", "cov_percentile"]].copy()
            output_df["species"] = speciesDict[file]
            all_outfile.write(output_df.to_csv(header=False, sep="\t", index=False) + "\n")
        dataDict[file] = covcheck

    for i in np.arange(0, 102.5, 2.5):
        bottom = 0
        columnDict = {}
        for key, covcheck in dataDict.items():

            height = 0
            color = (float(np.random.random(1)), float(np.random.random(1)), float(np.random.random(1)))
            percList = covcheck['cov_percentile'].tolist()
            for perc in percList:
                if int(perc/2.5) * 2.5 == i:
                    height += 1

            left = i
            width = 2.5
            # rectangle = mplpatches.Rectangle([left, bottom], width, height,
            #                                  facecolor=colorpalette[key],
            #                                  edgecolor='black',
            #                                  linewidth=0.1)
            # mainPanel.add_patch(rectangle)
            columnDict[key] = (width, height, left)

            if i == 0:
                keyRectangle = mplpatches.Rectangle([0, keyBottom], 1, 1,
                                                facecolor=colorpalette[key],
                                                edgecolor='black',
                                                linewidth=0.1)

                keyPanel.add_patch(keyRectangle)
                plt.gcf().text(0.825, textBottom, speciesDict[key], fontsize=5)
                textBottom += 0.045
                keyBottom += 1

        for key, value in sorted(columnDict.items(), key=lambda x: x[1][1]):
            # print(i, key, value[2], bottom, value[0], value[1])
            rectangle = mplpatches.Rectangle([left, bottom], value[0], value[1],
                                             facecolor=colorpalette[key],
                                             edgecolor='black',
                                             linewidth=0.1)
            bottom += value[1]
            mainPanel.add_patch(rectangle)

    # save whole figure
    plt.savefig('pleasefuckingwork.png', dpi=600)

def main():
    calc_coverage()
    
if __name__ == '__main__':
    main()
