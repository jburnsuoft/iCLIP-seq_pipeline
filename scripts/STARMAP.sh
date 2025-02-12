#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --account=def-zhaolei
#SBATCH --mail-type=ALL
module load CCEnv
module load nixpkgs/16.09
module load intel/2018.3
module load star/2.7.1a
STAR --runMode alignReads --genomeDir ../STAR/hg19/  --readFilesIn $1 --runThreadN 80 --readFilesCommand zcat --outFilterMultimapNmax 1 --outFileNamePrefix ../mapped_output/$1.output --outSAMtype BAM SortedByCoordinate --outSAMattributes All
