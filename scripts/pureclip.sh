#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=5
#SBATCH --account=def-zhaolei
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=END
conda activate condaclip
pureclip -i $1 -bai $1.bai -g ../../reference/hg19.fa -nt 100 -iv 'chr1;chr2;chr3;' -o ../pureclip/$1.pureclip.peaks.bed -or ../pureclip/$1.pureclip.regions.bed 

