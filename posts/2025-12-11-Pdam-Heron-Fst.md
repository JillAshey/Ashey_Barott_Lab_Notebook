---
layout: post
title: Pdam Heron Fst
date: '2025-12-11'
categories: Analysis
tags: [Bioinformatics, Fst]
projects: Pdam Heron
---

## Pdam Heron population connectivity 

Running fst bioinformatics for Angela/Marcelina on Unity using Pdam samples from [Brown et al. 2025](https://onlinelibrary.wiley.com/doi/10.1111/mec.17603). Make scratch directory to work in. 

```
ws_allocate pdam_tagseq_analysis 30
cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis
```

Zoe was working with these fastq files so they are already on Unity! Located here: `/project/pi_hputnam_uri_edu/raw_sequencing_data/20230125_Barott_Pdam`.

Code from Marcelina as reference: 

```
#!/bin/bash
#SBATCH --job-name=Pdam_Fst
#SBATCH --output=Pdam_Fst_%j.log
#SBATCH --error=Pdam_Fst_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
# Load modules (modify according to your HPC environment)
module load bwa
module load samtools
module load angsd
# Set paths
REF="Pocillopora_acuta_HIv2.assembly.fasta"
FASTQ_DIR="/path/to/FASTQ"
BAM_DIR="/path/to/BAM"
POP1="pop1_samples.txt"  # list of FASTQ files or sample names for pop1
POP2="pop2_samples.txt"  # list of FASTQ files or sample names for pop2
OUT_DIR="/path/to/output"
mkdir -p ${BAM_DIR} ${OUT_DIR}
# ------------------------------
# Step 1: Align FASTQ to reference
# ------------------------------
for fq in ${FASTQ_DIR}/*_R1.fastq.gz; do
    sample=$(basename ${fq} _R1.fastq.gz)
    fq2="${FASTQ_DIR}/${sample}_R2.fastq.gz"
    bam="${BAM_DIR}/${sample}.sorted.bam"
    echo "Aligning sample ${sample}"
    bwa mem -t 8 ${REF} ${fq} ${fq2} | samtools view -b - | \
        samtools sort -@ 8 -o ${bam}
    # Index BAM
    samtools index ${bam}
done
# ------------------------------
# Step 2: Create BAM lists for ANGSD
# ------------------------------
ls ${BAM_DIR}/*.bam > ${OUT_DIR}/pop_all.bamlist
# If populations are separate, you can create pop-specific lists:
grep -f ${POP1} ${OUT_DIR}/pop_all.bamlist > ${OUT_DIR}/pop1.bamlist
grep -f ${POP2} ${OUT_DIR}/pop_all.bamlist > ${OUT_DIR}/pop2.bamlist
# ------------------------------
# Step 3: Run ANGSD to calculate SAF
# ------------------------------
# Example for all samples
angsd -b ${OUT_DIR}/pop_all.bamlist \
      -ref ${REF} \
      -anc ${REF} \
      -out ${OUT_DIR}/pop_all \
      -doSaf 1 \
      -GL 2 \
      -doMajorMinor 1 \
      -doMaf 1 \
      -minMapQ 30 \
      -minQ 20 \
      -minInd 20 \
      -minDepth 100 \
      -maxDepth 20000 \
      -nThreads 8
# ------------------------------
# Step 4: Estimate the site frequency spectrum (SFS)
# ------------------------------
realSFS ${OUT_DIR}/pop_all.saf.idx -P 8 > ${OUT_DIR}/pop_all.sfs
# ------------------------------
# Step 5: Fst between two populations
# ------------------------------
# First, calculate SAF for each population
angsd -b ${OUT_DIR}/pop1.bamlist -ref ${REF} -anc ${REF} \
      -out ${OUT_DIR}/pop1 -doSaf 1 -GL 2 -doMajorMinor 1 -doMaf 1 \
      -minMapQ 30 -minQ 20 -minInd 10 -minDepth 50 -maxDepth 20000 -nThreads 8
angsd -b ${OUT_DIR}/pop2.bamlist -ref ${REF} -anc ${REF} \
      -out ${OUT_DIR}/pop2 -doSaf 1 -GL 2 -doMajorMinor 1 -doMaf 1 \
      -minMapQ 30 -minQ 20 -minInd 10 -minDepth 50 -maxDepth 20000 -nThreads 8
# Calculate 2D SFS
realSFS ${OUT_DIR}/pop1.saf.idx ${OUT_DIR}/pop2.saf.idx -P 8 > ${OUT_DIR}/pop1_pop2.sfs
# Calculate Fst per site
realSFS fst index ${OUT_DIR}/pop1.saf.idx ${OUT_DIR}/pop2.saf.idx \
       -sfs ${OUT_DIR}/pop1_pop2.sfs -fstout ${OUT_DIR}/pop1_pop2 -P 8
# Print Fst in windows (optional)
realSFS fst print ${OUT_DIR}/pop1_pop2.fst.idx > ${OUT_DIR}/pop1_pop2.fst.txt 
```

I'm going to run it in several steps, as I anticipate that each step will be pretty RAM/time intensive. 

Samples notes:

- 48 samples total 
- Library prep and sequencing done at UT Austin
- Fastq files are single end, not paired end
	- “Sequencing was completed targeting standard coverage of 3–5 million 100-bp single-end reads per sample (Illumina NovaSeq 600 SR100).” (Brown et al., 2025)

### Align files with BWA

Software used and versions:

- BWA (v0.7.17)
- samtools (v1.19.2)

Align untrimmed files with bwa. `nano align_bwa.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis

module load bwa/0.7.17
module load samtools/1.19.2

echo "Index Pacuta genome for BWA" $(date)

cd /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/
bwa index Pocillopora_acuta_HIv2.assembly.fasta -p Pacuta_bwa

echo "Index complete, aligning sequences" $(date)

READDIR="/project/pi_hputnam_uri_edu/raw_sequencing_data/20230125_Barott_Pdam"
OUTDIR="/scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/bwa_alignments"
REFDIR="/work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/"

mkdir -p "$OUTDIR"

cd "$READDIR"

for R1 in *"_R1_001.fastq.gz"; do
    sample=${R1%%_R1_001.fastq.gz}
    echo "Processing $sample"    
    bwa mem "$REFDIR/Pacuta_bwa" "$R1" \
      | samtools sort -@ 8 -O BAM -o "$OUTDIR/${sample}.sorted.bam" -
    samtools index "$OUTDIR/${sample}.sorted.bam"
done

echo "Sequence alignment complete" $(date)
```

Submitted batch job 50461201

See [bwa manual](https://bio-bwa.sourceforge.net/bwa.shtml) and [samtools manual](https://www.htslib.org/doc/samtools.html) for more info about these programs. 

### Organize bam files by population 

There are two Pdam populations of interest - Reef Slope (RS) and Reef Flat (RF). The bam files have this information in the sample name. Make lists of the bam files for these two populations.

```
cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/bwa_alignments

# RS population only
ls RS*.bam > pop_RS.bamlist

# RF population only  
ls RF*.bam > pop_RF.bamlist
```
