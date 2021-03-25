#!/usr/bin/env bash

blasta=hgtvnr.out
fasta=prot.faa

ProteinGrouper.py -file $blasta

for file in *.txt
do
 outfile=${file%.txt}.fna
 faSomeRecords $fasta $file $outfile
 blastdbcmd -db /Scratch/smnieves/databases/nt/nt -entry_batch $file >> $outfile 
done

find . -name "*.fna" | parallel -j 20 "java -jar ~/bin/macse_v2.03.jar -prog alignSequences -seq {} > {}.macse" 

