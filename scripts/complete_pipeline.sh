#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --cpus-per-task=40
#SBATCH --account=def-zhaolei
#SBATCH --mail-type=ALL

##########################################
# PIPELINE DOCUMENTATION
##########################################
# INPUTS:
#  $1 - Raw sequencing fastq file for trimming (e.g., sample_RNA.fastq.gz)
#  $2 - Pureclip BED file for pureclip-to-fasta conversion (e.g., sample.bed)
#  $3 - Extension length in nucleotides for pureclip and BEDPREP steps (e.g., 50)
#  $4 - Top number of pureclip peaks to select (e.g., 500)
#  $5 - Input BAM file for PURECLIP analysis (e.g., sample.bam)
#  $6 - BED file for OVERLAP analysis (e.g., overlap_sample.bed)
#  $7 - BED file for BEDPREP analysis (e.g., pureclip_sample.bed)
#
# OUTPUTS:
#  - Trimmed and processed FASTQ files in "trimmed/" and "trimmed/processed/"
#  - STAR index (if generated) in directory "hs88" and STAR genome build in $SCRATCH/STAR
#  - Mapped reads in "mapped_output" and deduplicated BAMs in "bam_dedup"
#  - PURECLIP outputs in "pureclip/" (peaks and regions BED files)
#  - Fasta files from pureclip-to-fasta in appropriate subdirectories
#  - Overlap analysis results (e.g., directories "frac_top500" and "jaccard_top500")
#
# EXAMPLE USAGE:
#   bash /MSc-pipeline/scripts/complete_pipeline.sh sample_RNA.fastq.gz sample.bed 50 500 sample.bam overlap_sample.bed pureclip_sample.bed
##########################################

# Activate environment
conda activate condaclip

##########################################
# INTEGRATED PIPELINE SECTION
##########################################

echo "=== INTEGRATED PIPELINE ==="

# --- STEP 1: TRIM & UMI EXTRACTION ---
if [ -z "$1" ]; then
    echo "Usage: $0 <raw_fastq.gz> <pureclip_bed> <extension> <top_peaks> <input_bam> <overlap_bed> <bedprep_bed>"
    exit 1
fi

echo "Running cutadapt trimming on $1..."
cutadapt --times 1 -e 0.15 -O 1 --quality-cutoff 10 -m 15 \
    -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG \
    -o results/trimmed/"$1" data/"$1"

echo "Extracting UMIs..."
umi_tools extract --stdin=results/trimmed/"$1" --bc-pattern=NNNNNNNNN --log="$1".log --stdout results/processed/"$1"

# --- STEP 2: STAR REFERENCE INDEX & GENOME BUILD ---
echo "Checking for STAR index hs88..."
if [ ! -d hs88 ]; then
    echo "STAR index hs88 not found. Building index..."
    STAR --runMode genomeGenerate --runThreadN 80 --genomeDir hs88 \
         --genomeFastaFiles homo_sapiens.88.fa --genomeSAindexNbases 13 \
         --genomeSAsparseD 2 --outFileNamePrefix hs88 \
         --alignSJoverhangMin 8 --sjdbGTFfile homo_sapiens.88.gtf --sjdbOverhang 100
else
    echo "STAR index hs88 exists. Skipping index generation."
fi

echo "Building STAR genome with hg19..."
module load star/2.7.1a
STAR --runThreadN 38 --runMode genomeGenerate --genomeDir "$SCRATCH/STAR" \
    --genomeFastaFiles "$SCRATCH/clip/hg19.fa" \
    --sjdbGTFfile "$SCRATCH/clip/hg19.gtf" --sjdbOverhang 83

# --- STEP 3: MAPPING WITH STAR ---
echo "Mapping processed reads with STAR..."
cd results/processed
mkdir -p ../../mapped_output
for f in *.fastq.gz; do
    STAR --runMode alignReads --genomeDir ../STAR/hg19/ \
         --readFilesIn "$f" --runThreadN 80 --readFilesCommand zcat \
         --outFilterMultimapNmax 1 \
         --outFileNamePrefix ../../mapped_output/"${f%.fastq.gz}".output \
         --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done
cd ../../

# --- STEP 4: DEDUPLICATION WITH UMI-TOOLS ---
echo "Deduplicating mapped BAMs..."
cd mapped_output
mkdir -p ../bam_dedup
for f in *.bam; do
    umi_tools dedup -I "$f" -S ../bam_dedup/"$f"
done
cd ../

echo "Integrated pipeline completed."

##########################################
# POST-PROCESSING SECTION
##########################################
echo "=== POST-PROCESSING PIPELINE ==="

# Define common parameters (adjust as needed)
PURECLIP_THREADS=100
PURECLIP_IV="chr1;chr2;chr3;"
EXTENSION_FA=50    # for pureclip-to-fasta
TOP_PEAKS=500
OVERLAP_EXT=75     # for overlap and peakfile overlap steps

# --- STEP 1: PURECLIP (Inline code) ---
echo "Running PURECLIP on deduplicated BAM files..."
mkdir -p results/pureclip
for bam in results/bam_dedup/*.bam; do
    echo "Processing $bam with PURECLIP..."
    # Ensure BAM index exists
    if [ ! -f "${bam}.bai" ]; then
        samtools index "$bam"
    fi
    # Run pureclip (inline command; adjust parameters as needed)
    pureclip -i "$bam" -bai "${bam}.bai" -g "$SCRATCH/clip/hg19.fa" -nt $PURECLIP_THREADS \
        -iv "$PURECLIP_IV" \
        -o results/pureclip/$(basename "$bam").pureclip.peaks.bed \
        -or results/pureclip/$(basename "$bam").pureclip.regions.bed
done

# --- STEP 2: PURECLIP TO FASTA (Inline code) ---
echo "Converting PURECLIP peaks BED to FASTA..."
for bed in results/pureclip/*.pureclip.peaks.bed; do
    echo "Converting $bed to FASTA..."
    # Extend peaks using bedtools slop
    bedtools slop -i "$bed" -b $EXTENSION_FA -g hg19.chrom.sizes.txt > "${bed}.${EXTENSION_FA}.bed"
    # Sort by score (assuming column 5 holds the score)
    sort -g -k5,5 -r "${bed}.${EXTENSION_FA}.bed" > "${bed}.${EXTENSION_FA}.sorted.bed"
    # Extract top peaks and get fasta sequences
    head -n $TOP_PEAKS "${bed}.${EXTENSION_FA}.sorted.bed" > "${bed}.${EXTENSION_FA}.sorted.top${TOP_PEAKS}.bed"
    bedtools getfasta -fi hg19.fa -bed "${bed}.${EXTENSION_FA}.sorted.top${TOP_PEAKS}.bed" -fo "${bed}.top${TOP_PEAKS}.fa"
done

# --- STEP 3: OVERLAP ANALYSIS (Inline code) ---
echo "Performing OVERLAP analysis on PURECLIP peaks..."
for bed in results/pureclip/*.pureclip.peaks.bed; do
    echo "Analyzing overlap for $bed..."
    # Create output directory using the overlap extension
    mkdir -p "${OVERLAP_EXT}.slopped"
    bedtools slop -i "$bed" -b $OVERLAP_EXT -g $HOME/hg19.chrom.sizes.txt > "${bed}.${OVERLAP_EXT}.bed"
    # (Additional overlap analysis steps such as computing intersections can be added here)
done

# --- STEP 4: BEDPREP (Inline code) ---
echo "Running BEDPREP on PURECLIP regions..."
for bed in results/pureclip/*.pureclip.regions.bed; do
    echo "Preparing BED for $bed..."
    mkdir -p "${EXTENSION_FA}.slopped"
    bedtools slop -i "$bed" -b $EXTENSION_FA -g $HOME/hg19.chrom.sizes.txt > "${bed}.${EXTENSION_FA}.bed"
    sort -g -k5,5 -r "${bed}.${EXTENSION_FA}.bed" > "${bed}.${EXTENSION_FA}.sorted.bed"
    head -n $TOP_PEAKS "${bed}.${EXTENSION_FA}.sorted.bed" > "${bed}.${EXTENSION_FA}.sorted.top${TOP_PEAKS}.bed"
    mv "${bed}.${EXTENSION_FA}.bed" "${EXTENSION_FA}.slopped/"
    mv "${bed}.${EXTENSION_FA}.sorted.bed" "${EXTENSION_FA}.slopped/"
    mv "${bed}.${EXTENSION_FA}.sorted.top${TOP_PEAKS}.bed" "${EXTENSION_FA}.slopped/"
done

# --- STEP 5: METAGENE_DEEPTOOLS ANALYSIS (Inline code) ---
echo "Running METAGENE_DEEPTOOLS on deduplicated BAM files..."
for bam in results/bam_dedup/*.bam; do
    echo "Analyzing $bam with METAGENE_DEEPTOOLS..."
    # Generate BAM index if missing
    if [ ! -f "${bam}.bai" ]; then
        samtools index "$bam"
    fi
    # Generate coverage file
    bamCoverage --bam "$bam" -o "${bam}.bamcoverage.bw" \
        --binSize 50 --normalizeUsing RPGC --effectiveGenomeSize 2864785220 \
        --numberOfProcessors max --verbose
    # Compute matrix for scale-regions mode
    computeMatrix scale-regions -S "${bam}.bamcoverage.bw" -R /Volumes/Monolith/clip/hg19.gtf \
        -a 500 -b 500 --verbose -o "${bam}.computematrix.scaled.gz"
    # Plot profile
    plotProfile -m "${bam}.computematrix.scaled.gz" \
        -out "${bam}.computematrix.scaled.profile.png" \
        --numPlotsPerRow 2 --plotTitle "$bam profile"
done

# --- STEP 6: PEAKFILE OVERLAPS (Inline code) ---
echo "Running PEAKFILE OVERLAPS analysis..."
for bed in results/pureclip/*.pureclip.peaks.bed; do
    echo "Running peakfile overlaps for $bed..."
    mkdir -p "${OVERLAP_EXT}.slopped.peakfile"
    bedtools slop -i "$bed" -b $OVERLAP_EXT -g $HOME/hg19.chrom.sizes.txt > "${bed}.${OVERLAP_EXT}.bed"
    # (Additional commands to compute overlaps between peakfiles may be included here)
    mv "${bed}.${OVERLAP_EXT}.bed" "${OVERLAP_EXT}.slopped.peakfile/"
done

echo "Post-processing pipeline completed successfully."