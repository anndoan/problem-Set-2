#! /usr/bin/env bash

datasets='/vol2/home/anndoan/data-sets/'
# Question 1: Use BEDtools intersect to identify the size of the largest
# overlap between CTCF and H3K4me3 locations.

CTCF="$datasets/bed/encode.tfbs.chr22.bed.gz"

H3K4me3="$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz"

answer_1=$(bedtools intersect -a $CTCF -b $H3K4me3 \
    |awk 'BEGIN {OFS ="\t"} ($4 == "CTCF") {print $0, $3 -$2}' \
    |sort -k5nr |head -n1 |cut -f5)

echo "answer-1: $answer_1"

# Question 2: Use BEDtools to calculate the GC content of nucleotides
# 19,000,000 to 19,000,500 on chr22 of hg19 genome build. Report the GC
# content as a fraction (e.g., 0.50).

#first: create a bed file and send it to output file called question2.bed.
#echo -e "chr22\t19000000\t19000500" >> question2.bed 

FASTA="$datasets/fasta/hg19.chr22.fa"
Question2="$datasets/bed/question2.bed"

answer_2=$(bedtools nuc -fi $FASTA -bed $Question2 \
    |cut -f5 \
    |tail -n1)

echo "answer-2: $answer_2"

#Question 3: Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e., interval) that has the largest mean signal in
#ctcf.hela.chr22.bg.gz.

CTCFChIPseq="$datasets/bed/encode.tfbs.chr22.bed.gz"
ctcfhela="$datasets/bedtools/ctcf.hela.chr22.bg.gz"

answer_3=$(bedtools map -c 4 -o mean -a $CTCFChIPseq -b $ctcfhela \
    |awk 'BEGIN {OFS ="\t"} ($4 == "CTCF") {print $0, $3 - $2}' \
    |sort -k5nr |head -n1 |cut -f6)

echo "answer-3: $answer_3"

#Question 4: Use BEDtools to identify the gene promoter (defined as 1000 bp upstream of a TSS) with the highest median signal in
#ctcf.hela.chr22.bg.gz. Report the gene name (e.g., 'ABC123')

TSS="$datasets/bed/tss.hg19.chr22.bed.gz"
genome="$datasets/genome/hg19.genome"

answer_4=$(bedtools slop -l 1000 -r -1 -s -i $TSS -g $genome \
    |bedtools sort -i \
    |bedtools map -a - -b $ctcfhela -c 4 -o median -null 0 \
    |sort -k7nr \
    |head -n1 \
    |cut -f4)

echo "answer-4: $answer_4"

#Question 5: Use BEDtools to identify the longest interval on chr22 that is not covered by genes.hg19.bed.gz. Report the interval like
#chr1:100-500.

genes="$datasets/bed/genes.hg19.bed.gz"

answer_5=$(bedtools complement -i $genes -g $genome \
    |bedtools sort -i \
    |awk 'BEGIN {OFS ="\t"} ($1 == "chr22") {print $0, $3 - $2}' \
    |sort -k4nr \
    |head -n1 \
    |awk '{print $1":"$2"-"$3}')

echo "answer-5: $answer_5"    

#Question 6: Extra credit - Use one or more BEDtools that we haven't covered in class. Be creative.
#Use bedtools flank instead bedtools slop to answer question 4 above.

answer_6=$(bedtools flank -l 1000 -r -1 -s -i $TSS -g $genome \
    |bedtools sort -i \
    |bedtools map -a - -b $ctcfhela -c 4 -o median -null 0 \
    |sort -k7nr \
    |head -n1 \
    |cut -f4)

echo "answer-6: $answer_6"


