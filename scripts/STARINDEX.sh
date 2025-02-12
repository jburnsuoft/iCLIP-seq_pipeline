#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --account=def-zhaolei
#SBATCH --cpus-per-task=40
conda activate condaclip
STAR   --runMode genomeGenerate   --runThreadN 80   --genomeDir hs88   --genomeFastaFiles homo_sapiens.88.fa      --genomeSAindexNbases 13   --genomeSAsparseD 2   --outFileNamePrefix hs88   --alignSJoverhangMin 8   --sjdbGTFfile homo_sapiens.88.gtf   --sjdbOverhang 100
