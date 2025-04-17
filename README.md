# Whole-Exome Sequencing (WES) Analysis Pipeline for Osteopetrosis Variant Detection

This repository contains a pipeline for processing whole-exome sequencing (WES) data from a family trio (proband, father, mother) to identify candidate variants associated with osteopetrosis, a rare genetic bone disorder. The pipeline starts with raw FASTQ files, performs alignment, variant calling, annotation, and filtering, and generates a report of potential causative variants under autosomal recessive and compound heterozygous inheritance models.

## Overview

The pipeline consists of 12 steps, implemented as Bash scripts, to process WES data for a trio where the proband is affected by osteopetrosis, and the consanguineous parents are unaffected. The goal is to identify rare, high-impact variants in known osteopetrosis genes.

## Prerequisites

- **System**: Linux (tested on Ubuntu)
- **Conda**: For managing environments and dependencies
- **Docker**: For running GATK tools
- **Reference Files**:
  - Human genome: `hg19.fa` (adjust path in scripts)
  - VEP cache: GRCh37 (`....../wes/ref/vep_cache`)
- **Tools** (installed via Conda or Docker):
  - FastQC
  - MultiQC
  - Trimmomatic
  - BWA
  - Samtools
  - Picard
  - GATK (via Docker: `broadinstitute/gatk:latest`)
  - VEP
  - bgzip/tabix
  - vcf2db.py
  - GEMINI
  - bcftools
- **Conda Environments**:
  - `wes_analysis`: For FastQC, MultiQC, Trimmomatic, BWA, Samtools, Picard, VEP, bcftools
  - `vcf2db_env`: For vcf2db.py and GEMINI

  ```bash
  conda create -n wes_analysis fastqc multiqc trimmomatic bwa samtools picard vep bcftools
  conda create -n vcf2db_env python=3.8 gemini vcf2db
  ```

## Pipeline Steps

The pipeline processes FASTQ files through alignment, variant calling, annotation, and variant filtering. Each step is implemented as a Bash script, with the first step (quality control) performed manually or scripted separately.

### Step 1: Download Data (`1_data.sh`)

- **Description**: Downloads raw FASTQ files from Zenodo.
- **Command**:

  ```bash
  ./1_data.sh
  ```
- **Output**: `father_R1.fq.gz`, `father_R2.fq.gz`, `mother_R1.fq.gz`, `mother_R2.fq.gz`, `proband_R1.fq.gz`, `proband_R2.fq.gz`

### Step 2: Quality Control

- **Tools**: FastQC, MultiQC
- **Description**: Assess raw FASTQ file quality for the trio (`father_R1.fq.gz`, `father_R2.fq.gz`, `mother_R1.fq.gz`, `mother_R2.fq.gz`, `proband_R1.fq.gz`, `proband_R2.fq.gz`).
- **Command**:

  ```bash
  mkdir fastqc_results
  fastqc -o fastqc_results *.fq.gz
  multiqc fastqc_results -o multiqc_report
  ```
- **Output**: `fastqc_results/*.html`, `multiqc_report/multiqc_report.html`
- 
### Step 3: Trim Reads (`2_trim.sh`)

- **Tool**: Trimmomatic
- **Description**: Removes low-quality bases and adapters from FASTQ files.
- **Command**:

  ```bash
  ./2_trim.sh
  ```
- **Output**: `${sample}_R1_trimmed.fq.gz`, `${sample}_R2_trimmed.fq.gz`, `${sample}_R1_unpaired.fq.gz`, `${sample}_R2_unpaired.fq.gz`, `${sample}_trimming.log`

### Step 4: Align Reads (`3_align.sh`)

- **Tool**: BWA
- **Description**: Aligns trimmed reads to the hg19 reference genome, producing SAM files.
- **Command**:

  ```bash
  ./3_align.sh
  ```
- **Output**: `${sample}.sam`

### Step 5: Filter Reads (`4_filter.sh`)

- **Tool**: Samtools
- **Description**: Filters SAM files to keep paired, mapped reads with high mapping quality (`-q 20`), converting to BAM.
- **Command**:

  ```bash
  ./4_filter.sh
  ```
- **Output**: `${sample}_filtered.bam`

### Step 6: Sort BAM Files (`5_sort.sh`)

- **Tool**: Samtools
- **Description**: Sorts and indexes filtered BAM files.
- **Command**:

  ```bash
  ./5_sort.sh
  ```
- **Output**: `${sample}_filtered_sorted.bam`, `${sample}_filtered_sorted.bam.bai`

### Step 7: Mark Duplicates (`6_dedup.sh`)

- **Tool**: Picard
- **Description**: Marks and removes duplicate reads from sorted BAM files.
- **Command**:

  ```bash
  ./6_dedup.sh
  ```
- **Output**: `${sample}_dedup.bam`, `${sample}_dedup_metrics.txt`

### Step 8: Index BAM Files (`7_index.sh`)

- **Tool**: Samtools
- **Description**: Creates index files for deduplicated BAMs.
- **Command**:

  ```bash
  ./7_index.sh
  ```
- **Output**: `${sample}_dedup.bam.bai`

### Step 9: Variant Calling (`8_calling.sh`)

- **Tool**: GATK (HaplotypeCaller, GenomicsDBImport, GenotypeGVCFs)
- **Description**: Performs per-sample variant calling (gVCFs), combines them, and conducts joint genotyping.
- **Command**:

  ```bash
  ./8_calling.sh
  ```
- **Output**: `${sample}.g.vcf.gz`, `trio_db`, `trio.vcf.gz`

### Step 10: Post-Process Variants (`9_post_proc.sh`)

- **Tool**: GATK (LeftAlignAndTrimVariants)
- **Description**: Splits multiallelic variants and normalizes indels in the joint VCF.
- **Command**:

  ```bash
  ./9_post_proc.sh
  ```
- **Output**: `trio_normalized_split.vcf.gz`

### Step 11: Annotate Variants (`10_annotate_VEP.sh`)

- **Tool**: VEP
- **Description**: Annotates variants with functional consequences, ClinVar significance, and population frequencies.
- **Command**:

  ```bash
  ./10_annotate_VEP.sh
  ```
- **Output**: `trio_annotated_vep.vcf.gz`, `trio_annotated_vep.vcf.gz.tbi`

### Step 12: Create GEMINI Database (`11_create_db.sh`)

- **Tools**: vcf2db.py, GEMINI
- **Description**: Converts the normalized VCF to a GEMINI SQLite database with pedigree information.
- **Command**:

  ```bash
  ./11_create_db.sh
  ```
- **Output**: `trio.ped`, `trio.db`

### Step 13: Find Candidate Variants (`12_find.sh`)

- **Tool**: bcftools
- **Description**: Filters variants for autosomal recessive (proband homozygous, parents heterozygous) and compound heterozygous patterns, focusing on high/moderate impact variants.
- **Command**:

  ```bash
  ./12_find.sh
  ```
- **Output**: `autosomal_recessive_candidates.vcf`, `autosomal_recessive_candidates.tsv`, `compound_het_temp.vcf`, `compound_het_candidates.tsv`

## File Structure

```
wes_pipeline/
├── 1_data.sh
├── 2_trim.sh
├── 3_align.sh
├── 4_filter.sh
├── 5_sort.sh
├── 6_dedup.sh
├── 7_index.sh
├── 8_calling.sh
├── 9_post_proc.sh
├── 10_annotate_VEP.sh
├── 11_create_db.sh
├── 12_find.sh
├── fastqc_results/
├── multiqc_report/
├── hg19.fa
├── father_R1.fq.gz
├── father_R2.fq.gz
├── mother_R1.fq.gz
├── mother_R2.fq.gz
├── proband_R1.fq.gz
├── proband_R2.fq.gz
├── *.sam
├── *.bam
├── *.vcf.gz
├── trio.ped
├── trio.db
├── autosomal_recessive_candidates.tsv
├── compound_het_candidates.tsv
```

## Usage

1. Clone the repository:

   ```bash
   git clone <repository-url>
   cd wes_pipeline
   ```
2. Set up Conda environments:

   ```bash
   conda env create -n wes_analysis -c bioconda fastqc multiqc trimmomatic bwa samtools picard vep bcftools
   conda env create -n vcf2db_env python=3.8 gemini vcf2db
   ```
3. Ensure `hg19.fa` and VEP cache are in `...../wes/ref/`.
4. Run scripts sequentially:

   ```bash
   bash 1_data.sh
   bash 2_trim.sh
   # ... continue through 12_find.sh
   ```
5. Perform quality control separately before Step 2.

## Notes

- **Sample Names**: Hardcoded as `father`, `mother`, `proband`. Adjust in scripts if different.
- **Reference Genome**: Uses `hg19.fa`. Update paths in scripts for your setup.
- **VEP Cache**: Assumes GRCh37 cache at `......./wes/ref/vep_cache`.
- **Troubleshooting**:
  - Check sample order in VCF headers for `12_find.sh` (`bcftools view -h trio_annotated_vep.vcf.gz | grep "^#CHROM"`).
  - If no candidates are found, relax filters in `12_find.sh` (e.g., remove gene restrictions).
