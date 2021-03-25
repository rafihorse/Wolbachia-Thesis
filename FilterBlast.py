#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

__author__ = 'Serafina Nieves'
__email__ = 'smnieves@ucsc.edu'

import argparse
import os
import pandas as pd


class BlastParser:
    """
    ArgParser handles command line arguments for the BLASTAnalyzer class to allow user specified filtering of BLAST data.  These arguments can be specified before
    running BLAST+, but these modules allow re-filtering of BLAST data without having to re-run lengthy computational programs.
    """

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Handle filtering arguments for the FilterBlast class.',
            prefix_chars='-',
            usage='prog [options] -option1[default] <input >output'
        )

        self.parser.add_argument('-file', action='store', help='input file name')
        self.parser.add_argument('-out', action='store', help='output file name')
        self.parser.add_argument('-percID', type=float, action='store', default=0,
                                 help='specify minimum % identity for a closer match')
        self.parser.add_argument('-len_align', type=int, action='store', default=0,
                                 help='specify minimum alignment length')
        self.parser.add_argument('-evalue', type=float, action='store', default=0,
                                 help='specify threshold e-value for more specific results')
        self.parser.add_argument('-bitscore', type=int, action='store', default=0,
                                 help='specify threshold bit score for more specific results')

        self.parser.add_argument('-dir', action='store',
                                 help='specify a valid path to directory containing BLAST output files')
        self.parser.add_argument('-ext', action='store', default='.out',
                                 help='specify file type extension to correctly locate files')

        self.args = self.parser.parse_args()


class FilterBlast:
    file_count = 0
    filtered_count = 0

    def __init__(self, file, out, evalue, bitscore, len_align, percID):
        """
        Receive the blast output file and parse fields using pandas.  Store command-line arguments from user.
        """

        self.file = file
        self.out = out
        self.evalue = evalue
        self.bitscore = bitscore
        self.len_align = len_align
        self.percID = percID

        try:
            self.blast_table = pd.read_table(self.file, engine='python', header=None,
                                             names=["query acc.ver", "subject acc.ver", "% identity",
                                                    "alignment length",
                                                    "mismatches", "gap opens", "q. start", "q. end", "s. start",
                                                    "s. end",
                                                    "evalue", "bit score"])
            self.blast_table.set_index("query acc.ver", inplace=True)
        except OSError:
            print(self.file + " produces an IO error.")

    def filter(self):
        filtered_blast = self.blast_table.loc[
            (self.blast_table['% identity'] >= self.percID) & (self.blast_table['evalue'] <= float(self.evalue)) & (
                    self.blast_table['alignment length'] >= self.len_align) & (
                    self.blast_table['bit score'] >= self.bitscore)].copy()
        filtered_blast.sort_values(by=['subject acc.ver', 'bit score'], ascending=[True, False], inplace=True)

        if filtered_blast.empty:
            FilterBlast.file_count += 1

        elif not filtered_blast.empty:

            FilterBlast.file_count += 1

            out = open(self.out, 'w')
            out.write(filtered_blast.to_csv(index=True, sep='\t', header=True))
            out.close()
            FilterBlast.filtered_count += 1


def main():
    print("Filtering files...")
    blast_parser = BlastParser()

    directory_in_str = blast_parser.args.dir
    directory = os.fsencode(directory_in_str)

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        outfile = filename + ".filtered"
        if filename.endswith(blast_parser.args.ext):
            blast = FilterBlast(directory_in_str + filename, directory_in_str + outfile, blast_parser.args.evalue,
                                blast_parser.args.bitscore, blast_parser.args.len_align, blast_parser.args.percID)

            blast.filter()
            del blast

    print('All significant hits saved.\n{0} total files.\n{1} filtered files matched user parameters.'.format(
        FilterBlast.file_count, FilterBlast.filtered_count))


if __name__ == "__main__":
    main()
