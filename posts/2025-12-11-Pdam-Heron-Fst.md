---
layout: post
title: Pdam Heron Fst
date: '2025-12-11'
categories: Analysis
tags: [Bioinformatics, Fst]
projects: Pdam Heron
---

## Pdam Heron population connectivity 

Running fst bioinformatics for Angela/Marcelina on Unity using Pdam samples from [Brown et al. 2025](https://onlinelibrary.wiley.com/doi/10.1111/mec.17603) (github repo [here](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master)). Make scratch directory to work in. 

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

### Sample information

There are 48 samples total from [Brown et al. 2025](https://onlinelibrary.wiley.com/doi/10.1111/mec.17603). Samples were stored in RNAlater and stored at -80C until extraction. RNA was extracted using the Zymo Quick-RNA/RNA Miniprep Plus kit (Zymo Research #D7003). Samples were bead-beat with 0.5mm glass beads at max speed for 1-2 min. Extractions were performed according to Zymo protocol, with a proteinase K digestion step for 15 mins at room temperature and a DNAse I treatment for 15 mins at room temperature. See [here](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/Project-Summary-Barott-and-Brown-Pdam-RNA-DNA-Extractions.md) for more information about RNA extractions. 

Library prep and sequencing was completed at UT Austin (see [here](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/TagSeq_Submission) for TagSeq submission information). Libraries were sequenced with targeting coverage of 3-5 million reads (pretty low coverage but typical for TagSeq data). Generated reads were 100-bp single-end. See the [QC report](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/data/raw_qc) of the raw samples for more information about the raw reads. 

Important to note: for each sample, there are two files: one includes the prefix "L001" and the other "L002". I'm not sure what the difference is between these or if the samples were just sequenced over multiple lanes. However, before proceeding with analysis, files for the same sample must be concatenated together (see Zoe [TagSeq pipeline](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/Heron-Pdam-gene-expression.md) for the same samples). 

### Concatenate L001 and L002 files for each sample 

`nano cat_files.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=50GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis

module load samtools/1.19.2

cd /project/pi_hputnam_uri_edu/raw_sequencing_data/20230125_Barott_Pdam

echo "Concatenate files" $(date)

\ls *R1_001.fastq.gz | awk -F '[_]' '{print $1"_"$2}' | sort | uniq > /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/ID

while read i; do
    cat ${i}_L001_R1_001.fastq.gz ${i}_L002_R1_001.fastq.gz > /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/${i}_ALL.fastq.gz
done < /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/ID

echo "Mission complete." $(date)
```

Submitted batch job 50497324

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
#SBATCH -t 60:00:00
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

READDIR="/scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/"
OUTDIR="/scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/bwa_alignments"
REFDIR="/work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/"

mkdir -p "$OUTDIR"

cd "$READDIR"

for FASTQ in *"_ALL.fastq.gz"; do
    sample=${FASTQ%%_ALL.fastq.gz}
    echo "Processing $sample"    
    bwa mem "$REFDIR/Pacuta_bwa" "$FASTQ" \
      | samtools sort -O BAM -o "$OUTDIR/${sample}.sorted.bam" -
    samtools index "$OUTDIR/${sample}.sorted.bam"
done

echo "Sequence alignment complete" $(date)
```

Submitted batch job 50502094

See [bwa manual](https://bio-bwa.sourceforge.net/bwa.shtml) and [samtools manual](https://www.htslib.org/doc/samtools.html) for more info about these programs. 

### Organize bam files by population 

There are two Pdam populations of interest - Reef Slope (RS) and Reef Flat (RF). The bam files have this information in the file name. Make lists of the bam files for these two populations.

```
cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/bwa_alignments

# RS population only
ls RS*sorted.bam > popRS.bamlist

# RF population only  
ls RF*sorted.bam > popRF.bamlist
```

Make list of all bam files as well. 

```
ls *sorted.bam > pop_all.bamlist
```

### Run ANGSD to determine depth

Software used and versions:

- ANGSD (v0.935)

Because our data is TagSeq (3-5M reads per sample, 3' bias), coverage will be relatively sparse and 3' biased (given the library prep methods). To choose depth filters, I'm going to run ANGSD's `-doDepth 1 -doCounts 1` before proceeding with the other analyses. This will generate depth histograms from the BAM files: `.depthGlobal` (total reads across all samples per site) and `.depthSample` (per-individual depths). These distributions will help us choose suitable depth parameters. 

`nano depth_angsd.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis

module load angsd/0.935

cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/bwa_alignments

echo "Determine sample depth to set downstream depth parameters" $(date)

angsd -bam pop_all.bamlist \
  -doCounts 1 -doDepth 1 \
  -minMapQ 20 -minQ 20 \
  -out pop_all_depth

echo "Sample depth run complete" $(date)
```

Submitted batch job 50546988. See [ANGSD manual](https://www.popgen.dk/angsd/index.php/ANGSD#Overview) and [github](https://github.com/ANGSD/angsd) for more info about the program. ANGSD will be used downstream to determine SAF and Fst. 

Look at the `pop_all_depth.depthGlobal` in R. This shows total read depth across samples

```{r}
getwd()

library(tidyverse)

# Read in data 
data <- read.table("pop_all_depth.depthGlobal", header=FALSE)
sites <- as.numeric(data[1,])  
depth <- 0:(length(sites)-1)

# Plot first 51 bins (0-50)
barplot(sites[1:51], names.arg=0:50, las=2, cex.names=0.7,
        xlab="Total Depth", ylab="# Sites", 
        main="Coverage Distribution", col="steelblue")

# Cumulative % retained
cum_sites <- rev(cumsum(rev(sites)))
pct <- 100 * cum_sites / sum(sites)
plot(depth[1:51], pct[1:51], type="l", lwd=2, ylim=c(0,100),
     xlab="maxDepth", ylab="% Sites Retained")
abline(h=70, col="red", lty=2)
abline(v=3, col="blue", lty=2)  # minDepth suggestion
```





















GL 1 or 2? 1 is Samtools, 2 is GATK. Leaning towards samtools since it is better with low coverage


Data not high overage

- https://academic.oup.com/g3journal/article/15/10/jkaf172/8219480


Other coral papers 

- https://onlinelibrary.wiley.com/doi/full/10.1111/eva.70115 -- they trimmed data before aligning 
- https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.70771 -- this used GL 1 (https://github.com/sanna2110/SYMBIO_WA/tree/main)
- https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12593?casa_token=Pw3tbUZWAD0AAAAA%3A0P6eZeCZRV6wo3G7VJ6rYB6aLZTsgWnDM3oeHPYVW2Jf03YpWelzOjb9fSQCNuOs-M5MrroWGuLIJzfu -- used GL 1 (https://datadryad.org/dataset/doi:10.5061/dryad.ft596)
- https://www.sciencedirect.com/science/article/pii/S004896972107501X?via%3Dihub#s0125 - used GL 1 )https://github.com/jamesfifer/JapanRE/blob/650f342abcc518f059771e20075c502e7c9f8e84/JapanRE_Temperate_Walkthrough.txt#L9)
- https://onlinelibrary.wiley.com/doi/full/10.1111/rec.70234?casa_token=hHNzHcmUgUEAAAAA%3AHmtnIcLNyQtNn5kzWY4Vd4kSInMu7LAXGy-QhhCeho-SJfX6AhfkAD-F47zo9b1_NbeBLedHz5E-sjyM - used GL 1 (https://github.com/Sydney-Bell/RTTproject.SCTLD)

ASK Jacob, Megan, Amy 