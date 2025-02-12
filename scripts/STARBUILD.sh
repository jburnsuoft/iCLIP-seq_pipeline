#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --account=def-zhaolei
#SBATCH --cpus-per-task=40
conda activate condaclip
module load CCEnv
module load nixpkgs/16.09
module load intel/2018.3
module load star/2.7.1a
module load samtools

STAR --runThreadN 38 --runMode genomeGenerate --genomeDir $SCRATCH/STAR --genomeFastaFiles $SCRATCH/clip/hg19.fa.masked --sjdbGTFfile $SCRATCH/clip/hg19.gtf --sjdbOverhang 83
