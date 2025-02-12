# iCLIP-seq_pipeline Documentation

## Project Overview
This project implements a bioinformatics pipeline designed for processing RNA sequencing data. The pipeline includes steps for trimming raw sequencing data, generating STAR indices, mapping reads, deduplicating BAM files, and performing analyses with PURECLIP and other bioinformatics tools.

## Directory Structure
```
MSc-pipeline
├── scripts
│   └── complete_pipeline.sh       # Shell script that orchestrates the entire pipeline
├── data
│   ├── sample_RNA.fastq.gz        # Raw sequencing data in FASTQ format
│   ├── sample.bed                 # BED file for pureclip-to-fasta conversion
│   ├── sample.bam                 # Aligned sequencing reads in BAM format
│   ├── overlap_sample.bed         # BED file for overlap analysis
│   └── pureclip_sample.bed        # BED file for BEDPREP analysis
├── results
│   ├── trimmed                     # Directory for trimmed FASTQ files
│   ├── processed                   # Directory for processed FASTQ files after UMI extraction
│   ├── mapped_output               # Directory for mapped reads in BAM format
│   ├── bam_dedup                   # Directory for deduplicated BAM files
│   ├── pureclip                    # Directory for PURECLIP analysis outputs
│   ├── frac_top500                 # Directory for overlap analysis results (top 500 peaks)
│   └── jaccard_top500              # Directory for additional overlap analysis results
├── .gitignore                      # Specifies files and directories to ignore in Git
└── README.md                       # Documentation for the project
```

## Usage Instructions
To run the pipeline, execute the following command in your terminal:

```bash
bash scripts/complete_pipeline.sh <raw_fastq.gz> <pureclip_bed> <extension> <top_peaks> <input_bam> <overlap_bed> <bedprep_bed>
```

### Parameters:
- `<raw_fastq.gz>`: Raw sequencing FASTQ file for trimming (e.g., `sample_RNA.fastq.gz`)
- `<pureclip_bed>`: Pureclip BED file for pureclip-to-fasta conversion (e.g., `sample.bed`)
- `<extension>`: Extension length in nucleotides for pureclip and BEDPREP steps (e.g., `50`)
- `<top_peaks>`: Top number of pureclip peaks to select (e.g., `500`)
- `<input_bam>`: Input BAM file for PURECLIP analysis (e.g., `sample.bam`)
- `<overlap_bed>`: BED file for OVERLAP analysis (e.g., `overlap_sample.bed`)
- `<bedprep_bed>`: BED file for BEDPREP analysis (e.g., `pureclip_sample.bed`)

## Outputs
The pipeline generates various outputs, including:
- Trimmed and processed FASTQ files in the `results/trimmed` and `results/processed` directories.
- Mapped reads in the `results/mapped_output` directory and deduplicated BAMs in `results/bam_dedup`.
- PURECLIP outputs in the `results/pureclip` directory (peaks and regions BED files).
- Overlap analysis results in the `results/frac_top500` and `results/jaccard_top500` directories.

## Requirements
- Conda environment with necessary bioinformatics tools installed (e.g., cutadapt, STAR, umi_tools, bedtools).
- Access to the required reference genome files (e.g., `homo_sapiens.88.fa`, `hg19.fa`, etc.).

## License
This project is licensed under the MIT License - see the LICENSE file for details.
