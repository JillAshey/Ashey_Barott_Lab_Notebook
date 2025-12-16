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
	- Cite Li and Durbin [2010](https://academic.oup.com/bioinformatics/article/26/5/589/211735?login=true)
- samtools (v1.19.2)
	- Cite Li et al. [2009](https://academic.oup.com/bioinformatics/article/25/16/2078/204688)

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

### Determine depth with ANGSD

Software used and versions:

- ANGSD (v0.935)
	- Cite Nielsen et al. [2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037558) and Korneliussen et al. [2014](https://link.springer.com/article/10.1186/s12859-014-0356-4)

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

Look at the `pop_all_depth.depthGlobal` in R. This shows total read depth across samples per genomic site.

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

The coverage distribution shows exactly what we would expect--a peak around 1-2 total depth for specific sites (ie low coverage sites), with a rapid tapering off. This is consistent with low coverage TagSeq. 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/depth_coverage_distribution.png)

If we set the minimum depth to 3 (as shown in the plot below), we retain the majority of covered sites and exclude the noise of 1-2 coverage. 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/depth_sites_retained.png)

### Calculate Site Allele Frequency likelihoods (SAF) with ANGSD (for all samples and per population)

Software used and versions:

- ANGSD (v0.935)
	- Cite Nielsen et al. [2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037558) and Korneliussen et al. [2014](https://link.springer.com/article/10.1186/s12859-014-0356-4)

SAF (Site Allele Frequency likelihoods) are probability distributions of allele counts at each genomic site across samples, based on genotype likelihoods from the BAM files. This is needed to estimate 1D SDS (diversity within a population) and 2D SDS (Fst between populations). 

I am going to run one job with all samples and one job with samples separated by population. 

`nano saf_all_angsd.sh`

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

module load angsd/0.935

cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/bwa_alignments

echo "Calculate SAF for all samples" $(date)

angsd -b pop_all.bamlist \
      -ref /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.assembly.fasta \
      -anc /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.assembly.fasta \
      -out ../pop_all \
      -doSaf 1 \
      -doCounts 1 \ 
      -GL 1 \
      -doMajorMinor 1 \
      -doMaf 1 \
      -minMapQ 20 \
      -minQ 20 \
      -minInd 24 \
      -setMinDepth 3 \
      -setMaxDepth 10000 

echo "SAF calculations for all samples complete" $(date)
```

Submitted batch job 50548654. 

`nano saf_per_pop_angsd.sh`

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

module load angsd/0.935

cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/bwa_alignments

echo "Calculate SAF for RF samples only" $(date)

angsd -b popRF.bamlist \
      -ref /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.assembly.fasta \
      -anc /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.assembly.fasta \
      -out ../popRF \
      -doSaf 1 \
      -doCounts 1 \ 
      -GL 1 \
      -doMajorMinor 1 \
      -doMaf 1 \
      -minMapQ 20 \
      -minQ 20 \
      -minInd 12 \
      -setMinDepth 3 \
      -setMaxDepth 10000 

echo "SAF calculations for RF samples complete" $(date)
echo "Calculate SAF for RS samples only" $(date)

angsd -b popRS.bamlist \
      -ref /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.assembly.fasta \
      -anc /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.assembly.fasta \
      -out ../popRS \
      -doSaf 1 \
      -doCounts 1 \ 
      -GL 1 \
      -doMajorMinor 1 \
      -doMaf 1 \
      -minMapQ 20 \
      -minQ 20 \
      -minInd 12 \
      -setMinDepth 3 \
      -setMaxDepth 10000 

echo "SAF calculations for RS samples complete" $(date)
```

Submitted batch job 50548673. 

There are lots of different arguments that can go into ANGSD. Here are explanatation of the arguments that I used (changed some from Marcelina's original code based on samples) and justification: 

- `-b` - list of input BAM files for specific population
- `-ref/anc` - reference + ancestral alleles. Since we don't have that information for our data, the reference genome was used as a proxy. 
- `-out` - prefix for output files.
- `-doSaf 1` - computes site allele frequency likelihoods, enabling calculation of SFS/Fst downstream. 
- `doCounts 1` - enables read counting to calculate site allele frequency likelihoods. 
- `-GL 1` - genotype likelihood framework (1 = Samtools model). 
	- There are several [options](https://www.popgen.dk/angsd/index.php/Genotype_Likelihoods) for the `-GL` argument. I selected 1 (the Samtools model), as this is most appropriate when mapping with BWA and with low coverage data. Option 2 is the GATK model, which is similar to the Samtools model, but works best with data aligned with gatk. According to the [ANGSD manual](https://www.popgen.dk/angsd/index.php/Genotype_Likelihoods), there is little difference in output when using 1 or 2 for new bam files. 
	- Other coral papers that I looked at also used `-GL 1`, so it appears to be the standard in the field (see Eriksson et al. [2025](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.70771) and associated [github](https://github.com/sanna2110/SYMBIO_WA/tree/main), Therkildsen & Palumbi [2016](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12593?casa_token=Pw3tbUZWAD0AAAAA%3A0P6eZeCZRV6wo3G7VJ6rYB6aLZTsgWnDM3oeHPYVW2Jf03YpWelzOjb9fSQCNuOs-M5MrroWGuLIJzfu) and associated [dryad](https://datadryad.org/dataset/doi:10.5061/dryad.ft596), Fifer et al. [2022](https://www.sciencedirect.com/science/article/pii/S004896972107501X?via%3Dihub) and associated [github](https://github.com/jamesfifer/JapanRE/blob/650f342abcc518f059771e20075c502e7c9f8e84/JapanRE_Temperate_Walkthrough.txt#L9), Bell et al. [2025](https://onlinelibrary.wiley.com/doi/full/10.1111/rec.70234?casa_token=hHNzHcmUgUEAAAAA%3AHmtnIcLNyQtNn5kzWY4Vd4kSInMu7LAXGy-QhhCeho-SJfX6AhfkAD-F47zo9b1_NbeBLedHz5E-sjyM) and associated [github](https://github.com/Sydney-Bell/RTTproject.SCTLD)). 
- `-doMajorMinor 1` and `-doMaf 1` - calls major and minor allele frequencies (may be needed downstream).
- `-minMapQ 20` and `-minQ 20` - mapping/base quality filters. I set both to 20 (meaning Phred score of 20 or more for mapping and base quality is needed for a site to proceed). Phred 20 is equivalent to ~1% error, and this is standard in popgen analysis. 
- `-minInd` - minimum number of individuals that have a particular SAF for that SAF to proceed in analysis. I set this number to be 50% of the population (ie 24 individuals across all 48 samples, 12 individuals across 24 samples per population). Again, this is standard for popgen, and I have seen some papers that use 70-80% of their samples (see Fifer et al. [2022](https://www.sciencedirect.com/science/article/pii/S004896972107501X?via%3Dihub) and Bell et al. [2025](https://onlinelibrary.wiley.com/doi/full/10.1111/rec.70234?casa_token=hHNzHcmUgUEAAAAA%3AHmtnIcLNyQtNn5kzWY4Vd4kSInMu7LAXGy-QhhCeho-SJfX6AhfkAD-F47zo9b1_NbeBLedHz5E-sjyM)). 
- `-setMinDepth 3` - total sequencing depth across samples must be at least 3 or higher. This will omit any noise from potential sequencing artifacts or lowly expressed genes. 
- `-setMaxDepth 10000` - total sequencing depth across samples cannot be greater than 10k. This is an absurdly high number, given the low coverage that we have, but I wanted to set it well above the threshold. 

The output files produced are the following: 

- `pop*.saf.pos.gz` - site positions (ie chromosome start positions) for each SAF entry. 
- `pop*.saf.idx` - binary index of SAFs. This is the file needed for our next step! 
- `	pop*.saf.gz` - raw likelihood data of SAFs. 
- `pop*.mafs.gz` - minor allele frequencies.
- `pop*.arg` - log file of parameters and runtime. 

### Compute 1D Site Frequency Spectrum (SFS)

Software used and versions:

- ANGSD (v0.935)
	- Cite Nielsen et al. [2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037558) and Korneliussen et al. [2014](https://link.springer.com/article/10.1186/s12859-014-0356-4)

Using the data with all the samples, regardless of population, I will now compress the SAF likelihood data into a 1D Site Frequency Spectrium (SFS), which shows the proportion of genome sites with derived alleles at each frequency (see ANGSD [manual](https://www.popgen.dk/angsd/index.php/RealSFS
) for more about SFS). This will help us determine the number of rare v. common alleles across all samples. Read more about this method in Han et al. [2015](https://academic.oup.com/bioinformatics/article/31/5/720/318186?login=true). 

`nano sfs_global_angsd.sh`

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

module load angsd/0.935

cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis

echo "Estimate 1D SFS for all samples" $(date)
realSFS pop_all.saf.idx > pop_all.sfs
echo "1D SFS estimation for all samples complete" $(date)

echo "Estimate 1D SFS for RF samples" $(date)
realSFS popRF.saf.idx > popRF.sfs
echo "1D SFS estimation for RF samples complete" $(date)

echo "Estimate 1D SFS for RS samples" $(date)
realSFS popRS.saf.idx > popRS.sfs
echo "1D SFS estimation for RS samples complete" $(date)
```

Submitted batch job 50552649. Plot output in R. 

Global SFS (all samples):

```{r}
sfs <- scan("pop_all.sfs")
barplot(sfs[-c(1,length(sfs))]/sum(sfs[-c(1,length(sfs))]), 
        xlab="Derived allele frequency", ylab="Proportion of polymorphic sites",
        main="Global SFS (all samples)")
```

![](https://github.com/JillAshey/Ashey_Barott_Lab_Notebook/blob/main/images/global_SFS.png?raw=true)

The x-axis is showing derived allele frequency classes across all samples, while the y-axis is showing the proportion of polymorphic sites in each allele class. Derived allele frequency is the fraction of chromosomes in your sample that carry the new (mutated) version of a site, assuming the reference/ancestral state is the “old” allele. A polymorphic site is a position in the genome where more than one allele is present in the population at non‑trivial frequency. In the global plot, we are seeing very tall bars at low derived allele frequencies, which means most polymorphic sites have a low derived allele frequency. So there are more sites where the allele is rate. 

RF population: 

```{r}
sfs_RF <- scan("popRF.sfs")
barplot(sfs_RF[-c(1,length(sfs_RF))]/sum(sfs_RF[-c(1,length(sfs_RF))]), 
        xlab="Derived allele frequency", ylab="Proportion of polymorphic sites",
        main="Global SFS (RF samples)")
```

![](https://github.com/JillAshey/Ashey_Barott_Lab_Notebook/blob/main/images/RF_SFS.png?raw=true)

RS population: 

```{r}
sfs_RS <- scan("popRS.sfs")
barplot(sfs_RS[-c(1,length(sfs_RS))]/sum(sfs_RS[-c(1,length(sfs_RS))]), 
        xlab="Derived allele frequency", ylab="Proportion of polymorphic sites",
        main="Global SFS (RS samples)")
```

![](https://github.com/JillAshey/Ashey_Barott_Lab_Notebook/blob/main/images/RS_SFS.png?raw=true)

Similar patterns for both populations. Data is also consistent with low coverage sequencing. 

### Calculate 2D SFS

Software used and versions:

- ANGSD (v0.935)
	- Cite Nielsen et al. [2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037558) and Korneliussen et al. [2014](https://link.springer.com/article/10.1186/s12859-014-0356-4)

I will also use realSFS to calculate the 2D SFSs for the RF and RS populations. This will allow us to examine the Fst between the two populations. 

`nano sfs_fst_angsd.sh`

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

module load angsd/0.935

cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis

echo "Estimate 2D SFS comparing RS and RF pops" $(date)
realSFS popRF.saf.idx popRS.saf.idx > RFpop_RSpop.sfs
echo "2D SFS estimation complete" $(date)

echo "Index SFS results in preperation for stats" $(date)
realSFS fst index popRF.saf.idx popRS.saf.idx -sfs RFpop_RSpop.sfs -fstout RFpop_RSpop

echo "Calculate Fst of the two populations" $(date)
realSFS fst stats RFpop_RSpop.fst.idx

echo "Print fst in window" $(date)
realSFS fst stats2 RFpop_RSpop.fst.idx -win 50000 -step 10000 > RFpop_RSpop.fst.txt

echo "2D SFS and Fst calculations complete" $(date)
```

Submitted batch job 50553096. From the Fst `stats` calculation: 

```
Calculate Fst of the two populations Tue Dec 16 01:36:29 UTC 2025
0.016514	0.101893
```

The first number is the unweighted Fst (0.016514), and the second number is the weighted Fst across sites (0.101893). The weighted and unweighted math is a bit confusing for me but this github [issue](https://github.com/ANGSD/angsd/issues/16) explains pretty well. The stats is using calculations derived from Reynolds et al. [1983](https://academic.oup.com/genetics/article-abstract/105/3/767/5996242?redirectedFrom=fulltext&login=true). 

I also generated stats using a sliding window of 50000bp with a 10000bp step. Let's look at the output: 

```
wc -l RFpop_RSpop.fst.txt
36067 RFpop_RSpop.fst.txt

head RFpop_RSpop.fst.txt
region  chr     midPos  Nsites
(215,1510)(10000,30979)(10000,60000)    Pocillopora_acuta_HIv2___Sc0000000      35000   1297    0.061345
(1257,3241)(30726,69999)(20000,70000)   Pocillopora_acuta_HIv2___Sc0000000      45000   1986    0.044385
(1257,3321)(30726,70079)(30000,80000)   Pocillopora_acuta_HIv2___Sc0000000      55000   2066    0.043606
(1511,4942)(60356,89492)(40000,90000)   Pocillopora_acuta_HIv2___Sc0000000      65000   3433    0.047008
(1511,5136)(60356,96305)(50000,100000)  Pocillopora_acuta_HIv2___Sc0000000      75000   3627    0.048293
(1511,6847)(60356,109999)(60000,110000) Pocillopora_acuta_HIv2___Sc0000000      85000   5338    0.046875
(3242,9248)(70000,119102)(70000,120000) Pocillopora_acuta_HIv2___Sc0000000      95000   6008    0.048543
(3322,9749)(85759,121796)(80000,130000) Pocillopora_acuta_HIv2___Sc0000000      105000  6429    0.054845
(4943,9749)(90035,121796)(90000,140000) Pocillopora_acuta_HIv2___Sc0000000      115000  4808    0.055371
```

Columns mean the following: 

- `region` - region on the chromosome where sites are found (?)
- `chr` - chromosome name from reference genome 
- `midPos` - midpoint of the sliding window 
- `Nsites` - number of sites (ie SNPs) within that window contribute to Fst
- Last column - mean Fst for that window. These can be compared to the global Fsts calculated above to see if particular regions are more differentiated. 

### PCA of genotype likelihoods 

Software used and versions:

- ANGSD (v0.935)
	- Cite Nielsen et al. [2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037558) and Korneliussen et al. [2014](https://link.springer.com/article/10.1186/s12859-014-0356-4)
- PCAngsd (v1.36.4)
	- Cite Meisner & Albrechtsen [2018](https://academic.oup.com/genetics/article/210/2/719/5931101)

I want to make a PCA from the genotype likelihoods using ANGSD (see [PCA info](https://www.popgen.dk/angsd/index.php/PCA)) and [PCAngsd](https://www.popgen.dk/software/index.php/PCAngsd). Through this, covariance estimates between individuals will be done using genotype likelihoods and then decomposed into principal components. This will help us see if populations are clustering together. 

To make this PCA, I need to get the genotype likelihoods in beagle format with ANGSD. 

`nano gl_beagle_all_angsd.sh`

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

module load angsd/0.935

cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis/bwa_alignments

echo "Calculate SAF for all samples in beagle format" $(date)

angsd -b pop_all.bamlist \
       -ref /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.assembly.fasta \
       -anc /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.assembly.fasta \
       -out ../genoGL_pop_all \
       -doSaf 1 \
      	-doCounts 1 \ 
		-GL 1 \
 		-doGlf 2 \               
 		-doMajorMinor 1 \
  		-doMaf 1 \
  		-SNP_pval 1e-6 \          
  		-minMapQ 20 -minQ 20 \
  		-minInd 24 \
  		-setMinDepth 3 -setMaxDepth 10000

echo "SAF beagle calculations for all samples complete" $(date)
```

Submitted batch job 50556304. I added the `-doGlf 1` and `-SNP_pval 1e-6` arguments here to ensure the output file was in beagle format and that the SNP calls had a pvalue cutoff. 

Install PCAngsd (see instructions and information on [github](https://github.com/Rosemeis/pcangsd) on Unity. 

```
cd /work/pi_hputnam_uri_edu/conda/envs/
module load conda/latest
git clone https://github.com/Rosemeis/pcangsd.git
conda env create \
  --prefix /work/pi_hputnam_uri_edu/conda/envs/pcangsd \
  -f pcangsd/environment.yml
conda activate pcangsd
```

`nano all_pcangsd.sh`

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

module load angsd/0.935
module load conda/latest
conda activate /work/pi_hputnam_uri_edu/conda/envs/pcangsd

cd /scratch4/workspace/jillashey_uri_edu-pdam_tagseq_analysis

echo "Perform PCA selection scan and estimate admixture proportions with pcangsd" $(date)

pcangsd --beagle genoGL_pop_all.beagle.gz --out all_pcangsd --selection --admix

echo "pcangsd complete" $(date)
conda deactivate
```

Submitted batch job 50556788. 

The output files are as follows: 

- `all_pcangsd.cov` - pairwise covariance matrix between all individuals, estimated from genotype likelihoods. This is the raw material for making a PCA, giving PC1, PC2, etc. When plotted, it shows clustering of genetic similarity among individuals. 
- `all_pcangsd.selection` - a matrix of selection statistics per SNP per PC. Higher values mean that SNP loads strongly on that PC. 
- `all_pcangsd.admix.4.Q` - admixture proportions (Q matrix) for K=4 genetic clusters. Here, each row is an individual and columns are clusters. This tells us, for each coral, what fraction of the aligned reads comes from each inferred ancestry (ie genetic clusters). K is the number of ancestral clusters in the admixture model, and is inferred from the number of PCs (or eigenvectors) deemed 'significant' by the program. Here, we had 3 significant PCs. In the PCAngsd framework, K = number of significant PCs + 1 (4 in our case). 
- `all_pcangsd.log` - log file for the run. 

I made some plots in R to double check that things look okay. 

PCA plot with the covariance matrix: 

```{r}
# Read covariance matrix
C <- as.matrix(read.table("all_pcangsd.cov"))

# Eigen decomposition
e <- eigen(C)
pcs <- e$vectors   # columns = PC1, PC2, ...
var_expl <- e$values / sum(e$values)  # proportion per PC
pc1_lab <- paste0("PC1 (", round(100 * var_expl[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(100 * var_expl[2], 1), "%)")

# Read in metadata
meta <- read.table("samples.txt", header=TRUE, sep=",")  # columns: ID, POP

# Make PCA plot 
plot(pcs[,1], pcs[,2],
     xlab = pc1_lab,
     ylab = pc2_lab,
     pch = 19,
     col = as.factor(meta$POP))
legend("topright",
       legend = levels(as.factor(meta$POP)),
       col = 1:length(levels(as.factor(meta$POP))), pch = 19)
# text(pcs[,1], pcs[,2],
#      labels = meta$ID,
#      pos = 1, cex = 0.7)
```

ADD PLOT

PCA plot looks very interesting...separation of the two populations across PC1 but some interesting groupings. It looks like, for both populations, there is some genetic divergence within the population itself. Would be worth looking at which samples are the outliers and if these samples had particularly low counts/mapping rates. 

Admixture (Q) plot with the admix data: 

```{r}
Q <- as.matrix(read.table("all_pcangsd.admix.4.Q"))

# Set colors 
cols <- c("steelblue","tomato","gold","darkgreen")

# index used to order Q 
pop <- meta$POP
mainK <- apply(Q, 1, which.max)
ord <- order(pop, mainK)
Qord <- Q[ord, ]
pop_ord <- pop[ord]

bp <- barplot(t(Qord),
              col = cols,
              border = NA, space = 0,
              xlab = "Individuals (ordered by POP and cluster)",
              ylab = "Ancestry proportion")

# bp = x positions for the center of each bar
# find contiguous blocks per population
runs <- rle(as.character(pop_ord))

start_idx <- cumsum(c(1, head(runs$lengths, -1)))
end_idx   <- cumsum(runs$lengths)
mid_x     <- (bp[start_idx] + bp[end_idx]) / 2

# Add labels under x-axis indicating population blocks
axis(side = 1, at = mid_x, labels = runs$values, tick = FALSE, line = 1)
```

ADD PLOT

This is also very interesting! It is showing the ancestry components across individuals (the colors represent the 4 ancestry components) and populations (RF individuals on left, RS individuals on right). We are seeing 5-ish groupings: almost pure blue individuals, almost pure red individuals, almost pure green individuals, mostly yellow with some green and blue individuals, and mostly yellow with some red and green individuals. These separate according to population as well. Some individuals may have intermediate ancestry? Not sure what to make of this. 

### Questions / next steps 

- Do we need explicit SNPs? Would need to rerun angsd with VCF calling 
- Will gene ext help with identifying SNPs? 
- https://github.com/emmastrand/EmmaStrand_Notebook/blob/7da9d9cdd0250cb68743491bba4f7a132c9fc66c/_posts/2023-03-06-SNP-Calling-with-RNASeq-data.md 
