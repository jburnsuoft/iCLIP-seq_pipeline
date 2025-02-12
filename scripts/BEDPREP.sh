#!/bin/sh

#  BEDPREP.sh
#  
#
#  Created by james burns on 2021-03-03.

## Requirements
    # bedtools v2.30.0

## Usage
    # sh BEDPREP.sh $1 $2 $3
        # $1 - Pureclip Bedfile
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

mv $1.$2.bed $2.slopped/
mv $1.$2.sorted.bed $2.slopped/
mv $1.$2.sorted.top$3.bed $2.slopped/
