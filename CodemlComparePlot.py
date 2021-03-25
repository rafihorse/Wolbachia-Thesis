#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

__author__ = 'Serafina Nieves'
__email__ = 'smnieves@ucsc.edu'

from Bio.Phylo.PAML import codeml
import sys
import numpy as np
from scipy.stats.distributions import chi2
import matplotlib.pyplot as plt
import pandas as pd

''' CodemlCompare.py compares the null and 2 ratio model of codeml to determine which model fits the data better. '''

# set parameters for the program (flags)
point_list = []
string_towrite = []
df = pd.DataFrame(columns = ['directory', 'gene_id', 'dN (model 0)', 'dS (model 0)', 'arth_dNdS (model 2)', 'wol_dNdS (model 2)', 'arth_dN (model 2)', 'arth_dS (model 2)',
    'wol_dN (model 2)', 'wol_dS (model 2)', 'lnL (model 0)', 'lnL (model 2)', 'p_value'])

outfile = "/Scratch/smnieves/alignments/CodemlModelCompare.txt"

with open("/home/smnieves/bin/candidates.txt") as cand_file:
    cand_list = [str(line.rstrip()) for line in cand_file]

print(cand_list)

num_list = range(1, 1673)
for num in num_list:
    try:
        wdir = "prot_" + str(num)
        acc_file = "all_prot.part-" + str(num) + "_filtered.acc"
        acc = str(open(acc_file).readline().rstrip())
        infile_null = wdir + "/paml_results_null.out"
        infile_alt = wdir + "/paml_results_alt.out"
        outfile2 = wdir + "/CodemlCompare_prot" + str(num) + ".txt"

        print(wdir)
        print(acc)

        # read in data and parse it for relevant values
        results_null = codeml.read(infile_null)
        results_alt = codeml.read(infile_alt)

        lnL_null = results_null.get("NSsites").get(0).get("lnL")
        lnL_alt = results_alt.get("NSsites").get(0).get("lnL")

        likelihood_ratio = -2 * (lnL_null - lnL_alt)
        p_value = chi2.sf(likelihood_ratio, 1)

        arth_omega = results_alt.get("NSsites").get(0).get("parameters").get("omega")[1]
        background_omega = results_alt.get("NSsites").get(0).get("parameters").get("omega")[0]
        dS_0 = results_null.get("NSsites").get(0).get("parameters").get("dS")
        dN_0 = results_null.get("NSsites").get(0).get("parameters").get("dN")
        dS_2 = results_alt.get("NSsites").get(0).get("parameters").get("dS")
        dN_2 = results_alt.get("NSsites").get(0).get("parameters").get("dN")
        arth_dS = (dN_2 - dS_2*background_omega)/(arth_omega - background_omega)
        wol_dS = (dS_2*arth_omega - dN_2)/(arth_omega - background_omega)
        wol_dN = background_omega*(dS_2*arth_omega - dN_2) / (arth_omega - background_omega)
        arth_dN = arth_omega*(dN_2 - dS_2*background_omega) / (arth_omega - background_omega)

        # create dataframe to contain all information in table

        df = df.append({'directory': wdir, 'gene_id': acc, 'dN (model 0)': dN_0, 'dS (model 0)': dS_0, 'arth_dNdS (model 2)': arth_omega, 'wol_dNdS (model 2)': background_omega,
        'arth_dS (model 2)': arth_dS, 'wol_dS (model 2)': wol_dS, 'arth_dN (model 2)': arth_dN, 'wol_dN (model 2)': wol_dN,
        'lnL (model 0)': lnL_null, 'lnL (model 2)': lnL_alt, 'p_value': p_value}, ignore_index=True)

        # create formatted strings and save points
        genes_interest = []
        arth_list = open("/home/smnieves/bin/arthqueries.txt").readlines()
        if p_value <= 0.05:
            point_list.append((background_omega, arth_omega, p_value, "red", dS_0, dS_2, wol_dS, arth_dS, acc))
            dir_num = wdir[5:]
            acc_file = "all_prot.part-" + dir_num + ".acc"
            for line in open(acc_file):
                if line in arth_list:
                    genes_interest.append(line.strip())
        elif p_value > 0.05:
            point_list.append((background_omega, arth_omega, p_value, "#1f77b4", dS_0, dS_2, wol_dS, arth_dS, acc))

        # string_towrite.append("directory: {0} \t likelihood ratio: {1} \t p-value: {2} \t arth_genes: {3} \n".format(wdir, likelihood_ratio, p_value, ",".join(genes_interest)))


    except (AttributeError, FileNotFoundError, IndexError, TypeError, ZeroDivisionError):
        continue

# change dimensions of the entire figure
figureHeight = 12
figureWidth = 6

plt.figure(figsize=(figureWidth, figureHeight))

mainPanelHeight = 4
mainPanelWidth = 4

pvaluePanelHeight = 2
pvaluePanelWidth = 4

dSPanelHeight = 2
dSPanelWidth = 4

mainPanel = plt.axes([0.175, 0.6, mainPanelWidth / figureWidth, mainPanelHeight / figureHeight])
pvaluePanel = plt.axes([0.175, 0.35, pvaluePanelWidth / figureWidth, pvaluePanelHeight / figureHeight])
dSPanel = plt.axes([0.175, 0.1, dSPanelWidth / figureWidth, dSPanelHeight / figureHeight])

# mainPanel Omega Plot

filtered_x = []
filtered_y = []
filtered_p = []
filtered_c = []
filtered_arthdS = []
filtered_dS2 = []

for point in point_list:
    if  point[0] < 10 and point[1] < 10 and point[6] <= 3 and point[6] >= 0.01 and point[7] <= 3 and point[7] >= 0.01 and point[8] in cand_list:
    #if  point[0] < 10 and point[1] < 10 and point[5] <= 2 and point[5] >= 0.02 and point[8] in cand_list:
        filtered_x.append(point[0])
        filtered_y.append(point[1])
        filtered_p.append(point[2])
        filtered_c.append(point[3])
        filtered_arthdS.append(point[6])
        filtered_dS2.append(point[5])

print(len(filtered_p))
#print(len(filtered_arthvwol))

mainPanel.plot([0, 2], [0, 2],
                   linestyle='dashdot',
                   color='black',
                   linewidth=1,
                   zorder=1
                   )

mainPanel.scatter(filtered_x, filtered_y, alpha=0.25, c=filtered_c, zorder=2)

mainPanel.set_xlim(0, 2)
mainPanel.set_ylim(0, 2)

mainPanel.set_title("Comparison of Omega Values for Background and Foreground", pad=10, loc="center", fontsize=12, fontweight='bold')
mainPanel.set_xlabel("Background (Wolbachia) Omega")
mainPanel.set_ylabel("Foreground (Arthropod) Omega")

# pvaluePanel p-value Distribution Plot
pvaluePanel.hist(filtered_p, bins=np.arange(0, 1.05, 0.05))

pvaluePanel.set_title("P-value Distribution For All Gene Sets", pad=10, loc="center", fontsize=12, fontweight='bold')
pvaluePanel.set_xlabel("P-value Bins (Bin Size = 0.05)")
pvaluePanel.set_ylabel("Number of Gene Sets")


# comparing dS to p_value to look for trends
#dSPanel.scatter(filtered_p, filtered_arthvwol, alpha=0.25)
#dSPanel.scatter(filtered_p, filtered_dS2, alpha=0.25)
dSPanel.hist(filtered_dS2, bins=np.arange(0, 2.05, 0.05))

dSPanel.set_title("dS Distribution for All Gene Sets", pad=10, loc="center", fontsize=12, fontweight='bold')
dSPanel.set_xlabel("dS Bins")
dSPanel.set_ylabel("Number of Gene Sets")


#dSPanel.set_ylim(0, 10)
#dSPanel.set_xticks(np.arange(0, 1.05, 0.1))

plt.savefig('arthwoldS01_totaldShist.png', dpi=300)

df.to_csv('PAML_dataoverview.tsv', index=True, header=True, sep='\t')

# write formatted strings to files
#with open(outfile, 'w') as file:
#    for string in string_towrite:
#        file.write(string)

#with open(outfile2, 'w') as file2:
#    for string in string_towrite:
#        file2.write(string)
