#!/bin/sh

#  purecliptofasta.sh
#  
#
#  Created by james burns

## Reqired Programs
    # bedtools v2.30.0

## Required Reference Files For Initialization
    # Reference Genome in Fasta Format - Ex. hg19.fa
        # Index will be created automatically - Ex. hg19.fai
    # Chromosome Siles File in TXT Format - Ex. hg19.chrom.sizes.txt

## User Inputs
    # $1 - Pureclip Bedfile - Accepts Peak Files and Region Files
    # $2 - Number of nt to Extend Up/Downstream
    # $3 - Number of Top Pureclip Peaks to Extract
    # Example Command Line: sh purecliptofasta.sh SP1.bam.pureclip.peaks.bed 50 500


## Command Sequence
# Input Pureclip Output > Extend $2 Nucleotides In Each Direction
bedtools slop -i $1 -b $2 -g hg19.chrom.sizes.txt > $1.$2.bed

# Sort All Peaks by PureClip Score
sort -g -k5,5 -r $1.$2.bed > $1.$2.sorted.bed

# Extract Top User Provided Coordinates ($3) to New Bedfile
head -n $3 $1.$2.sorted.bed > $1.$2.sorted.top$3.bed

# Extract Sorted Multi-fasta Sequence Files
bedtools getfasta -fi hg19.fa -bed $1.$2.sorted.top$3.bed -fo $1.$2.sorted.top$3.fa
bedtools getfasta -fi hg19.fa -bed $1.$2.sorted.bed -fo $1.$2.sorted.fa

#Cleanup - Make Directory For Outputs
mkdir -p $1_fastas

#Cleanup Pt2 - Move Fastas/Bedfiles to Folder
mv $1.$2.* $1_fastas
