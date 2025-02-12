#!/bin/bash
# Activate environment
conda activate condaclip

# Optionally run the integrated pipeline if not already completed
# bash /workspaces/MSc/integrated_pipeline.sh

##########################################
# Step 1: Run PURECLIP on deduplicated BAM files
##########################################
for bam in /workspaces/MSc/bam_dedup/*.bam; do
    echo "Running PURECLIP on $bam"
    bash /workspaces/MSc/pureclip.sh "$bam"
done

##########################################
# Step 2: Convert PURECLIP peaks BED to FASTA
##########################################
for bed in /workspaces/MSc/pureclip/*.pureclip.peaks.bed; do
    echo "Converting PURECLIP peaks $bed to FASTA"
    # Example: extend 50 nt, extract top 500 peaks
    bash /workspaces/MSc/purecliptofasta.sh "$bed" 50 500
done

##########################################
# Step 3: Run OVERLAP analysis on PURECLIP peaks
##########################################
for bed in /workspaces/MSc/pureclip/*.pureclip.peaks.bed; do
    echo "Running OVERLAP analysis on $bed"
    # Example: use extension of 75 nt and cutoff of 1000
    bash /workspaces/MSc/OVERLAP.sh "$bed" 75 1000
done

##########################################
# Step 4: Run BEDPREP on PURECLIP regions
##########################################
for bed in /workspaces/MSc/pureclip/*.pureclip.regions.bed; do
    echo "Running BEDPREP on $bed"
    # Example: extend 50 nt and cutoff 500
    bash /workspaces/MSc/BEDPREP.sh "$bed" 50 500
done

##########################################
# Step 5: Run METAGENE_DEEPTOOLS analysis
##########################################
for bam in /workspaces/MSc/bam_dedup/*.bam; do
    echo "Running METAGENE_DEEPTOOLS on $bam"
    bash /workspaces/MSc/METAGENE_DEEPTOOLS_TRUELOCAL.sh "$bam"
done

##########################################
# Step 6: Run PEAKFILE_OVERLAPS analysis
##########################################
for bed in /workspaces/MSc/pureclip/*.pureclip.peaks.bed; do
    echo "Running PEAKFILE_OVERLAPS on $bed"
    # Example: use extension of 75 nt and cutoff of 1000
    bash /workspaces/MSc/PEAKFILE_OVERLAPS "$bed" 75 1000
done

echo "ICLIP Post-Processing Pipeline completed successfully with METAGENE_DEEPTOOLS and PEAKFILE_OVERLAPS analysis."
