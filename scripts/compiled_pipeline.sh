#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --account=def-zhaolei
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=ALL

# Pipeline Documentation:
# -----------------------
# INPUTS:
#   $1 - Raw sequencing fastq file (e.g., sample_RNA.fastq.gz) used by the TRIM step.
#   $2 - Pureclip bed file (e.g., sample.bed) for pureclip-to-fasta.
#   $3 - Number of nucleotides to extend for pureclip and BEDPREP.
#   $4 - Number of top pureclip peaks to select.
#   $5 - Input BAM file for PURECLIP.
#   $6 - Bed file for OVERLAP step.
#   $7 - Bed file for BEDPREP step.
#
# OUTPUTS:
#   - Trimmed reads in "trimmed" and "trimmed/processed" directories.
#   - STAR index directory "hs88" (if generated) and STAR genome build in $SCRATCH/STAR.
#   - PURECLIP outputs in "../pureclip" (peaks and regions BED files).
#   - FASTA and BED outputs in a "{input}_fastas" directory from pureclip-to-fasta.
#   - Overlap results in directories "frac_top<cutoff>" and "jaccard_top<cutoff>".
#   - Demultiplexed files in "demultiplex" directory.
#   - Mapped reads in "mapped_output" and deduplicated BAMs in "bam_dedup".
#
# EXAMPLE USAGE:
#   Input Example:
#     sample_RNA.fastq.gz, sample.bed, 50, 500, sample.bam, overlap_sample.bed, pureclip_sample.bed
#
#   Output Example:
#     trimmed/sample_RNA.fastq.gz and trimmed/processed/sample_RNA.fastq.gz,
#     A STAR index in directory "hs88", STAR genome build in $SCRATCH/STAR,
#     ../pureclip/sample.bam.pureclip.peaks.bed and ../pureclip/sample.bam.pureclip.regions.bed,
#     Fasta files in sample.bed_fastas, Overlap matrices in "frac_top500" and "jaccard_top500",
#     Demultiplexed fastq files in "demultiplex", Mapped outputs in "mapped_output", 
#     and deduplicated BAM files in "bam_dedup".

echo "Activating environment..."
conda activate condaclip

##########################################
# STEP 1: TRIM (from trim.sh)
##########################################
echo "STEP 1: Running trim.sh section"
# Usage: supply input file as first argument
if [ -z "$1" ]; then
    echo "Usage: $0 <input_fastq.gz>"
    exit 1
fi
# ...existing setup...
cutadapt --times 1 -e 0.15 -O 1 --quality-cutoff 10 -m 15 \
    -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG \
    -o trimmed/$1 $1

umi_tools extract --stdin=trimmed/$1 --bc-pattern=NNNNNNNNN --log=$1.log --stdout trimmed/processed/$1

##########################################
# STEP 2: STAR INDEX GENERATION (from STARINDEX.sh)
##########################################
echo "STEP 2: Running STARINDEX.sh section"
if [ -d hs88 ]; then
    echo "STAR index 'hs88' exists, skipping STAR index generation."
else
    # Generates a genome index using STAR (version set for Homo sapiens.88)
    # ...existing code...
    STAR --runMode genomeGenerate --runThreadN 80 --genomeDir hs88 \
         --genomeFastaFiles homo_sapiens.88.fa --genomeSAindexNbases 13 \
         --genomeSAsparseD 2 --outFileNamePrefix hs88 \
         --alignSJoverhangMin 8 --sjdbGTFfile homo_sapiens.88.gtf --sjdbOverhang 100
fi

##########################################
# STEP 3: STAR GENOME BUILD (from STARBUILD.sh)
##########################################
echo "STEP 3: Running STARBUILD.sh section"
# Builds STAR genome index using reference hg19.
# ...existing code...
module load CCEnv
module load nixpkgs/16.09
module load intel/2018.3
module load star/2.7.1a
module load samtools
STAR --runThreadN 38 --runMode genomeGenerate --genomeDir $SCRATCH/STAR \
    --genomeFastaFiles $SCRATCH/clip/hg19.fa.masked \
    --sjdbGTFfile $SCRATCH/clip/hg19.gtf --sjdbOverhang 83

##########################################
# STEP 4: PURECLIP TO FASTA (from purecliptofasta.sh)
##########################################
echo "STEP 4: Running purecliptofasta.sh section"
# Converts pureclip BED output to fasta sequences.
# Extend, sort, extract top peaks, get fasta, and move outputs.
bedtools slop -i $2 -b $3 -g hg19.chrom.sizes.txt > $2.$3.bed
sort -g -k5,5 -r $2.$3.bed > $2.$3.sorted.bed
head -n $4 $2.$3.sorted.bed > $2.$3.sorted.top$4.bed
bedtools getfasta -fi hg19.fa -bed $2.$3.sorted.top$4.bed -fo $2.$3.sorted.top$4.fa
bedtools getfasta -fi hg19.fa -bed $2.$3.sorted.bed -fo $2.$3.sorted.fa
mkdir -p ${2}_fastas
mv $2.$3.* ${2}_fastas

##########################################
# STEP 5: PURECLIP (from pureclip.sh)
##########################################
echo "STEP 5: Running pureclip.sh section"
# Runs pureclip on the provided input BAM file.
pureclip -i $5 -bai $5.bai -g ../../reference/hg19.fa -nt 100 \
    -iv 'chr1;chr2;chr3;' -o ../pureclip/$5.pureclip.peaks.bed -or ../pureclip/$5.pureclip.regions.bed

##########################################
# STEP 6: OVERLAP (from OVERLAP.sh)
##########################################
echo "STEP 6: Running OVERLAP.sh section"
# Extends bed peaks, sorts, extracts top coordinates, generates overlap matrices.
mkdir -p $3.slopped
bedtools slop -i $6 -b $3 -g $HOME/hg19.chrom.sizes.txt > $6.$3.bed
sort -g -k5,5 -r $6.$3.bed > $6.$3.sorted.bed
head -n $4 $6.$3.sorted.bed > $6.$3.sorted.top$4.bed
rm $6.$3.bed $6.$3.sorted.bed
mv $6.$3.sorted.top$4.bed $3.slopped/$6
intervene pairwise -i $3.slopped/*.bed --compute frac --htype color -o frac_top$4/
intervene pairwise -i $3.slopped/*.bed --compute jaccard --htype color --sort -o jaccard_top$4/

##########################################
# STEP 7: DEMULTIPLEX (from demultiplex.sh)
##########################################
echo "STEP 7: Running demultiplex.sh section"
# Demultiplexes reads based on ADAR barcodes.
conda activate condaclip  # ensure environment if needed
cutadapt --action=none --no-indels -e 0 \
    -g ADAR1_R4_2=^NNNAACCNN -g ADAR1_R3_2=^NNNTTGTNN \
    -o "demultiplex/{name}.fastq.1.gz" ADAR.fastq.gz

##########################################
# STEP 8: BED PREPARATION (from BEDPREP.sh)
##########################################
echo "STEP 8: Running BEDPREP.sh section"
# Extends Pureclip bed files and sorts for downstream analysis.
mkdir -p $2.slopped
bedtools slop -i $7 -b $2 -g $HOME/hg19.chrom.sizes.txt > $7.$2.bed
sort -g -k5,5 -r $7.$2.bed > $7.$2.sorted.bed
head -n $3 $7.$2.sorted.bed > $7.$2.sorted.top$3.bed
mv $7.$2.bed $2.slopped/
mv $7.$2.sorted.bed $2.slopped/
mv $7.$2.sorted.top$3.bed $2.slopped/

##########################################
# STEP 9: INTEGRATED PIPELINE (from integrated_pipeline.sh)
##########################################
echo "STEP 9: Running integrated_pipeline.sh section"
# Set working directory and module loads (if separate from earlier steps)
cd /scratch/z/zhaolei/jburns
module load CCEnv
module load StdEnv
source ~/.virtualenvs/clippipe/bin/activate
module load fastqc/0.11.9
module load fastx-toolkit/0.0.14
module load nixpkgs/16.09
module load samtools/1.10  ##SHOULD BE 1.5
module load bedtools/2.27.1
module load gcc/7.3.0
module load star
module load intel/2018.3
module load seqtk/1.2
module load bowtie2

####STEP 1 - DEMULTIPLEX (Integrated)
echo "Integrated STEP 1: Demultiplexing..."
mkdir -p demultiplex
cutadapt --action=none --no-indels -e 0 \
    -g G3BP1_N_A=^NNNGGCANN -g G3BP1_N_B=^NNNTTAANN -g G3BP1_N_C=^NNNAATANN \
    -g G3BP1_N_INPUT_A=^NNNCCACNN -g G3BP1_N_ARSENITE_A=^NNNCCGGNN \
    -g G3BP1_N_ARSENITE_B=^NNNTGGCNN -g G3BP1_N_ARSENITE_C=^NNNGGTCNN \
    -g G3BP1_N_ARSENITE_INPUT_A=^NNNCGGANN -g ADAR1_A=^NNNAACCNN -g ADAR1_B=^NNNTTGTNN \
    -o "demultiplex/{name}.fastq.1.gz" ADAR.fastq.gz
echo "Integrated STEP 1 completed."

####STEP 2 - ADAPTER TRIM + EXTRACT BARCODES (Integrated)
echo "Integrated STEP 2: Adapter trimming and barcode extraction..."
cd demultiplex
mkdir -p trimmed processed
for f in *.fastq.gz; do
    cutadapt --times 1 -e 0.15 -O 1 --quality-cutoff 10 -m 15 \
        -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG \
        -o trimmed/$f $f
    umi_tools extract --stdin=trimmed/$f --bc-pattern=NNNNNNNNN --log=$f.log --stdout processed/$f
done
echo "Integrated STEP 2 completed."

####STEP 3 - MAPPING WITH STAR (Integrated)
echo "Integrated STEP 3: Mapping processed reads with STAR..."
cd processed
module load CCEnv
module load nixpkgs/16.09
module load intel/2018.3
module load star/2.7.1a
mkdir -p ../mapped_output
for f in *.fastq.gz; do
    STAR --runMode alignReads --genomeDir ../STAR/hg19/ --readFilesIn $f \
         --runThreadN 80 --readFilesCommand zcat --outFilterMultimapNmax 1 \
         --outFileNamePrefix ../mapped_output/$f.output --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes All
done
echo "Integrated STEP 3 completed."

####STEP 4 - DEDUPLICATION WITH UMI-TOOLS (Integrated)
echo "Integrated STEP 4: Deduplicating with UMI-tools..."
cd ../mapped_output
mkdir -p ../bam_dedup
for f in *.bam; do
    umi_tools dedup -I $f -S ../bam_dedup/$f
done
echo "Integrated STEP 4 completed."

echo "Compiled pipeline completed successfully."
