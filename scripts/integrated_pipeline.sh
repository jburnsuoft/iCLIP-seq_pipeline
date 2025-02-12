#!/bin/bash

#  integrated_pipeline.sh
#  
#
#  Created by james burns on 2021-01-25 for use on the Niagara Cluster at SciNet.
#  

# Set working directory and load modules
cd /scratch/z/zhaolei/jburns
module load CCEnv
module load StdEnv
source ~/.virtualenvs/clippipe/bin/activate
module load fastqc/0.11.9
module load fastx-toolkit/0.0.14
module load nixpkgs/16.09
module load samtools/1.10 ##SHOULD BE 1.5
module load bedtools/2.27.1
module load gcc/7.3.0
module load star
module load intel/2018.3
module load seqtk/1.2
module load bowtie2


####STEP 1 - DEMULTIPLEX
echo "Starting Step 1: Demultiplexing..."
# Note: Separates reads into distinct files based on defined barcode sequences.
mkdir -p demultiplex
cutadapt --action=none --no-indels -e 0 -g G3BP1_N_A=^NNNGGCANN -g G3BP1_N_B=^NNNTTAANN -g G3BP1_N_C=^NNNAATANN -g G3BP1_N_INPUT_A=^NNNCCACNN -g G3BP1_N_ARSENITE_A=^NNNCCGGNN -g G3BP1_N_ARSENITE_B=^NNNTGGCNN -g G3BP1_N_ARSENITE_C=^NNNGGTCNN -g G3BP1_N_ARSENITE_INPUT_A=^NNNCGGANN -g ADAR1_A=^NNNAACCNN -g ADAR1_B=^NNNTTGTNN -o "demultiplex/{name}.fastq.1.gz" ADAR.fastq.gz
echo "Step 1 completed."

####STEP 2 - ADAPTER TRIM + EXTRACT BARCODES
echo "Starting Step 2: Adapter trimming and barcode extraction..."
# Note: Trims adapter sequences and extracts UMI barcodes for downstream analysis.
cd demultiplex
mkdir -p trimmed processed
for f in *.fastq.gz; do
    cutadapt --times 1 -e 0.15 -O 1 --quality-cutoff 10 -m 15 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG -o trimmed/$f $f
    umi_tools extract --stdin=trimmed/$f --bc-pattern=NNNNNNNNN --log=$f.log --stdout processed/$f
done
echo "Step 2 completed."

##########################################
# STEP 2.5: STAR REFERENCE SETUP (INDEX & GENOME BUILD)
##########################################
if [ ! -d "../STAR/hg19/" ]; then
    echo "STAR index not found. Building STAR genome index..."
    if [ -f "$SCRATCH/clip/hg19.fa" ] && [ -f "$SCRATCH/clip/hg19.gtf" ]; then
        STAR --runThreadN 80 --runMode genomeGenerate --genomeDir ../STAR/hg19/ \
             --genomeFastaFiles "$SCRATCH/clip/hg19.fa" \
             --sjdbGTFfile "$SCRATCH/clip/hg19.gtf" --sjdbOverhang 83
    else
        echo "ERROR: Missing STAR reference genome or GTF file in \$SCRATCH/clip"
        exit 1
    fi
else
    echo "STAR index found. Skipping STAR genome index generation."
fi

####STEP 3 - MAPPING WITH STAR (Integrated)
echo "Integrated STEP 3: Mapping processed reads with STAR..."
cd processed
module load CCEnv
module load nixpkgs/16.09
module load intel/2018.3
module load star/2.7.1a
mkdir -p ../mapped_output
for f in *.fastq.gz; do
    # Inline STAR mapping code (replacing call to STARMAP.sh)
    STAR --runMode alignReads --genomeDir ../STAR/hg19/ \
         --readFilesIn "$f" \
         --runThreadN 80 --readFilesCommand zcat \
         --outFilterMultimapNmax 1 \
         --outFileNamePrefix ../mapped_output/"${f%.fastq.gz}".output \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes All
done
echo "Integrated STEP 3 completed."

####STEP 4 - DEDUPLICATE WITH UMI-TOOLS
echo "Starting Step 4: Deduplicating with UMI-tools..."
# Note: Removes duplicate reads based on UMI sequences to ensure unique alignments.
cd ../mapped_output
mkdir -p ../bam_dedup
for f in *.bam; do
    umi_tools dedup -I $f -S ../bam_dedup/$f
done
echo "Step 4 completed."

echo "Pipeline completed successfully."