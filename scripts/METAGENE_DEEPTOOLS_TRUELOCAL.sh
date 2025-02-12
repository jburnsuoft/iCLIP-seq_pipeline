#!/bin/sh
#SBATCH --time=6:00:00
#SBATCH --nodes=2
#SBATCH --account=def-zhaolei

#  METAGENE_DEEPTOOLS.sh
#  
#
#  Created by james burns on 2021-09-21.
#  

cd /scratch/z/zhaolei/jburns/clip/bam_dedup
conda activate condaclip
module load samtools

#INDEX WITH SAMTOOLS
samtools index $1

#CALCULATE bamCoverage
bamCoverage --bam $1 -o $1.bamcoverage.bw \
    --binSize 50 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2864785220 \
    --numberOfProcessors max \
    --verbose
    
#COMPUTEMATRIX
computeMatrix scale-regions -S $1.bamcoverage.bw -R /Volumes/Monolith/clip/hg19.gtf -a 500 -b 500 --verbose -o $1.computematrix.scaled.gz

#computeMatrix reference-point \ # choose the mode
       #--referencePoint TSS \ # alternatives: TES, center
       #-b 3000 -a 10000 \ # define the region you are interested in
       #-R hg19.gtf \
       #-S $1.bamcoverage.bw  \
       #--skipZeros \
       #-o $1.computematrix_TSS.gz \ # to be used with plotHeatmap and plotProfile
       #--outFileSortedRegions $1.computematrix_TSS_genes.bed

plotProfile -m $1.computematrix.scaled.gz \
              -out $1.computematrix.scaled.gz.profile.png \
              --numPlotsPerRow 2 \
              --plotTitle "$1 profile"
