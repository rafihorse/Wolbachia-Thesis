\#!/usr/bin/env bash

while getopts s:r:m option
do
case "${option}"
in
s) SRA=${OPTARG};;
r) ref=${OPTARG};;
m) mem=${OPTARG};;
esac
done

echo "1. download SRA data"
fastq-dump --split-files --gzip $SRA

echo "2. prepare files for BWA"
bwa index $ref.fna

echo "3. map SRA reads to genome"
bwa mem -t 8 $ref.fna ${SRA}_1.fastq.gz ${SRA}_2.fastq.gz | samtools sort -@8 -o ${SRA}_sortedaln.bam -

echo "4. find base coverage of bed regions"
gff2bed <$ref.gff >$ref.bed
#cat $ref.gff | awk '!/^[ \t]*#/&&NF{print $1, $4, $5, $3, $8, $7, $9}' OFS='\t' >> $ref.bed
samtools bedcov $ref.bed ${SRA}_sortedaln.bam | awk '{print $1, $2, $3, $4, $5, $6, $7}' OFS='\t' > ${SRA}_sortedaln.cov

echo "5. generate covcheck files"
grep -Ff ~/bin/gene_ids.txt ${SRA}_sortedaln.cov > ${SRA}_HGT.cov
CovCheck.py -i $SRA

echo "bwaly.sh finished"
