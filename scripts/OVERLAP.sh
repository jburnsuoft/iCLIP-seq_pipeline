#!/bin/sh

#  OVERLAP.sh
#  
#
#  Created by james burns on 2021-03-03.

## Requirements
    # bedtools
    # intervene 

#SETUP
    #conda activate untitled
        #OR environment with bedtools and intervene installed
    
## Usage from current folder with .bed peak files
    # for f in *.bed; do sh OVERLAP.sh $f 75 1000; done
        # $1 - Pureclip Bedfiles
        # $2 - Extension #nt
        # $3 - Pureclip Cut-off


# Make Directory for Outputs
mkdir -p $2.slopped

# Input Pureclip Output > Extend $2 Nucleotides In Each Direction
bedtools slop -i $1 -b $2 -g $HOME/hg19.chrom.sizes.txt > $1.$2.bed

# Sort All Peaks by PureClip Score >
sort -g -k5,5 -r $1.$2.bed > $1.$2.sorted.bed

# Extract Top $3 Coordinates to New Bedfile
head -n $3 $1.$2.sorted.bed > $1.$2.sorted.top$3.bed

#Cleanup
rm $1.$2.bed
rm $1.$2.sorted.bed
mv $1.$2.sorted.top$3.bed $2.slopped/$1

##GENERATE MATRICES WITH INTERVENE

##FRACTIONAL OVERLAP
intervene pairwise -i $2.slopped/*.bed --compute frac --htype color -o frac_top$3/

##JACARD INDEX
intervene pairwise -i $2.slopped/*.bed --compute jaccard --htype color --sort -o jaccard_top$3/
