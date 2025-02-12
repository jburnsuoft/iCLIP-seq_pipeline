#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --nodes=2
#SBATCH --account=def-zhaolei
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=ALL
conda activate condaclip

cutadapt --times 1 -e 0.15 -O 1 --quality-cutoff 10 -m 15 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG -o trimmed/$1 $1

umi_tools extract --stdin=trimmed/$1 --bc-pattern=NNNNNNNNN --log=$1.log --stdout trimmed/processed/$1
