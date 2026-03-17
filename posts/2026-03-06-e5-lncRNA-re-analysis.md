---
layout: post
title: e5 deep dive lncRNA re-analysis
date: '2026-03-06'
categories: e5
tags: e5, ncRNA
---

I am working from the lncRNA fasta and bed files generated from our original analysis (see [github](https://github.com/urol-e5/deep-dive/tree/main)) and running some more analyses for stringency. This will likely affect the counts and downstream results, so these will have to be updated as well. Updated worflow: 

- Filter based on expression (ie at least 5 counts for 3 replicates)
- Run PLEK -- maybe not if it installs??
- Compare against Pfam 
- Update downstream workflows 
	- Number of unique v shared lncRNAs 
	- Length plot 
	- Gene lncRNA overlap and GO enrichment 

Install needed software for updated workflow

PLEK: 

```
# Download PLEK.1.2.tar.gz from https://sourceforge.net/projects/plek/files/
and decompress it.
tar zvxf PLEK.1.2.tar.gz
# Compile PLEK.
cd PLEK.1.2
python PLEK_setup.py

# Had some issues with the PLEK.py command--had to make some syntax edits and permissions adding 
g++ -O3 PLEK_main.c -o PLEK
g++ -O3 svm-train.c svm.cpp -o svm-train
g++ -O3 svm-predict.c svm.cpp -o svm-predict
chmod +x PLEK svm-train svm-predict svm-scale
sed -i "s/'rU'/'r'/g" PLEK.py
sed -i "s/sys.arg\[0\]/sys.argv[0]/g" PLEK.py
python PLEK.py -fasta PLEK_test.fa -out predicted_results.txt -thread 4
```

Download PFAM database 

```
cd /scratch3/workspace/jillashey_uri_edu-e5_deepdive
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
```

## Apul 

Turn lnc bed into gff 

```
awk 'BEGIN{OFS="\t"} 
{
    count++;
    # Logic: if start ($2) is negative, force GTF start to 1.
    # Otherwise, use standard BED-to-GTF conversion ($2 + 1).
    if ($2 < 0) {
        new_start = 1;
    } else {
        new_start = $2 + 1;
    }

    # Safety: Ensure end ($3) is always at least 1 and greater than start
    if ($3 <= new_start) { $3 = new_start + 1; }

    printf "%s\t.\tlncRNA\t%d\t%d\t.\t+\t.\tgene_id \"lncRNA_%03d\"; transcript_id \"lncRNA_%03d\";\n", $1, new_start, $3, count, count
}' Apul_lncRNA.bed > Apul_lncRNA_for_counts.gtf
```

Run feature counts to quantify lncRNAs

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/apul

module load conda/latest # need to load before making any conda envs
conda activate featurecounts

featureCounts \
  -p \
  -a Apul_lncRNA_for_counts.gtf \
  -t lncRNA \
  -g gene_id \
  -o Apul_lncRNA_counts_matrix.txt \
  *.bam
```

Submitted batch job 53420983

See how many lncRNAs are present with at least 5 counts in three biological reps 

```
awk -F'\t' '
  # 1. Print the header lines (anything starting with # or Geneid)
  /^#/ || /^Geneid/ { print; next } 
  
  # 2. For data lines, count how many samples have >= 5 counts
  {
    count = 0;
    for (i = 7; i <= 11; i++) {
      if ($i >= 5) count++;
    }
    
    # 3. If at least 3 samples met the criteria, print the line
    if (count >= 3) print;
  }' Apul_lncRNA_counts_matrix.txt > filtered_Apul_lncRNA_counts_matrix.txt
  
 wc -l filtered_Apul_lncRNA_counts_matrix.txt
7962 filtered_Apul_lncRNA_counts_matrix.txt
```

Went from ~16000 to 7962 lncRNAs from expression filtering alone. 

Get headers from filtered data and pull seqs from fasta 

```
# 1. Create a list of the exact header strings that passed the filter
# We take the filtered IDs, find them in the GTF, and format them to match the FASTA header
awk 'NR>2 {print $1}' filtered_Apul_lncRNA_counts_matrix.txt > Apul_filtered_ids.txt

grep -Fwf Apul_filtered_ids.txt Apul_lncRNA_for_counts.gtf | \
awk '{print "transcript::"$1":"$4"-"$5}' > Apul_matching_headers.txt

# 2. Use seqtk to pull those sequences
seqtk subseq Apul_lncRNA.fasta Apul_matching_headers.txt > exp_filtered_Apul_lncRNAs.fasta
```

Success. 

Get the new lncRNA names and coordinates

```
# Get the Geneids (lncRNA_XXX) from your filtered counts matrix
awk 'NR>2 {print $1}' filtered_Apul_lncRNA_counts_matrix.txt > Apul_filtered_ids.txt

# Extract the corresponding coordinates from your GTF using those IDs
# This creates a list of coordinates like "NC_058066.1:468618-469943"
grep -Fwf Apul_filtered_ids.txt Apul_lncRNA_for_counts.gtf | \
awk '{print $1":"$4"-"$5}' > Apul_filtered_coords.txt
```

Extract from original fasta using coordinates 

```
module load uri/main all/seqtk/1.4-GCC-12.3.0
seqtk subseq Apul_lncRNA.fasta Apul_filtered_coords.txt > exp_filtered_Apul_lncRNAs.fasta
```

Run PLEK on new fasta file. `nano apul_plek.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/apul

# Load conda
module load conda/latest
# Use the full path to the environment to ensure it finds it
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /work/pi_hputnam_uri_edu/conda/envs/PLEK.1.2

# Define the PLEK directory
PLEK_DIR="/work/pi_hputnam_uri_edu/conda/envs/PLEK.1.2"

# Run PLEK using the full path to the script
# We use -thread 8 to match the cpus-per-task requested above
python ${PLEK_DIR}/PLEK.py \
    -fasta exp_filtered_Apul_lncRNAs.fasta \
    -out apul_plek_results.txt \
    -thread 8
```

Submitted batch job 53544542. In the log file: `Coding: 1157/7957=14.540656026140505%, Non-coding: 6800/7957=85.45934397385949%`. 

Extract the non-coding IDs from the PLEK results 

```
awk -F'\t' '$1 == "Non-coding" {print $3}' apul_plek_results.txt | sed 's/^>//' > apul_plek_non_coding_ids.txt
```

Filter fasta based on non-coding IDs

```
module load uri/main all/seqtk/1.4-GCC-12.3.0
seqtk subseq exp_filtered_Apul_lncRNAs.fasta apul_plek_non_coding_ids.txt > plek_exp_filtered_Apul_lncRNAs.fasta
grep -c ">" plek_exp_filtered_Apul_lncRNAs.fasta
6800
```

Run hmmscan to check for protein coding regions. `nano apul_hmmscan.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/apul

module load uri/main bio/HMMER/3.3.2-gompi-2022a

echo "Hmmscan commencing" $(date)
echo "Unzip pfam file" $(date)

gunzip ../Pfam-A.hmm.gz

echo "press pfam db" $(date)
hmmpress ../Pfam-A.hmm

echo "Run hmmscan" $(date)
hmmscan --domtblout apul_hmmscan_output.domtblout ../Pfam-A.hmm plek_exp_filtered_Apul_lncRNAs.fasta > apul_hmmscan_output.txt

echo "hmmscan complete!" $(date)
```

Submitted batch job 53544913

Extract IDs that have protein domain hits

```
# Only remove hits with E-value < 0.001
grep -v '^#' apul_hmmscan_output.domtblout | awk '$7 < 0.001 {print $4}' | sort | uniq > apul_pfam_ids_to_remove.txt
```

Filter the fasta 

```
python3 -c "
import sys
# Load IDs to exclude into a set
with open('apul_pfam_ids_to_remove.txt') as f:
    exclude = set(line.strip() for line in f if line.strip())
with open('plek_exp_filtered_Apul_lncRNAs.fasta') as fasta:
    keep = True
    for line in fasta:
        if line.startswith('>'):
            # Check if the header ID is in our exclusion list
            header_id = line[1:].strip().split()[0]
            keep = header_id not in exclude
        if keep:
            print(line, end='')
" > apul_final_lncRNAs.fasta
grep -c ">" apul_final_lncRNAs.fasta
6798
```

Convert fasta headers to bed 

```
grep ">" apul_final_lncRNAs.fasta | sed 's/>transcript:://' | awk -F'[: -]' '{print $1 "\t" $2 "\t" $3 "\t" "transcript::" $1 ":" $2 "-" $3}' > apul_final_lncRNAs.bed
```

Sort the lncRNA bed file 

```
sort -k1,1 -k2,2n -k3,3n apul_final_lncRNAs.bed > apul_final_lncRNAs_sorted.bed
```

Sort gff file and filter for mRNA only 

```
sort -k1,1 -k4,4n GCF_013753865.1_Amil_v2.1_genomic.gff > GCF_013753865.1_Amil_v2.1_genomic_sorted.gff 

awk '$3 == "mRNA"' GCF_013753865.1_Amil_v2.1_genomic_sorted.gff > GCF_013753865.1_Amil_v2.1_genomic_sorted_mRNA_only.gff
```

Run bed closest 

```
module load bedtools2/2.31.1
bedtools closest -a apul_final_lncRNAs_sorted.bed -b GCF_013753865.1_Amil_v2.1_genomic_sorted_mRNA_only.gff > Apul_lncRNA_closest_mRNA_only.bed
```

Turn lnc final bed into gff 

```
awk 'BEGIN{OFS="\t"} 
{
    count++;
    # Logic: if start ($2) is negative, force GTF start to 1.
    # Otherwise, use standard BED-to-GTF conversion ($2 + 1).
    if ($2 < 0) {
        new_start = 1;
    } else {
        new_start = $2 + 1;
    }

    # Safety: Ensure end ($3) is always at least 1 and greater than start
    if ($3 <= new_start) { $3 = new_start + 1; }

    printf "%s\t.\tlncRNA\t%d\t%d\t.\t+\t.\tgene_id \"lncRNA_%03d\"; transcript_id \"lncRNA_%03d\";\n", $1, new_start, $3, count, count
}' apul_final_lncRNAs.bed > apul_final_lncRNA_for_counts.gtf
```

Run feature counts to quantify final lncRNAs. `nano apul_final_counts.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/apul

module load conda/latest # need to load before making any conda envs
conda activate featurecounts

featureCounts \
  -p \
  -a apul_final_lncRNA_for_counts.gtf \
  -t lncRNA \
  -g gene_id \
  -o Apul_final_lncRNA_counts_matrix.txt \
  *.bam
```

Submitted batch job 53608797

Calculate number of lncRNAs per sample

```
awk 'NR > 2 { 
    for (i=7; i<=NF; i++) 
        if ($i > 0) count[i]++ 
    } 
    END { 
        for (i=7; i<=NF; i++) 
            print "Sample " (i-6) ": " count[i] " lncRNAs" 
    }' Apul_final_lncRNA_counts_matrix.txt

Sample 1: 6373 lncRNAs
Sample 2: 4962 lncRNAs
Sample 3: 6262 lncRNAs
Sample 4: 5690 lncRNAs
Sample 5: 6596 lncRNAs
```

### Ptuh

Turn lnc bed into gff 

```
awk 'BEGIN{OFS="\t"} 
{
    count++;
    # Logic: if start ($2) is negative, force GTF start to 1.
    # Otherwise, use standard BED-to-GTF conversion ($2 + 1).
    if ($2 < 0) {
        new_start = 1;
    } else {
        new_start = $2 + 1;
    }

    # Safety: Ensure end ($3) is always at least 1 and greater than start
    if ($3 <= new_start) { $3 = new_start + 1; }

    printf "%s\t.\tlncRNA\t%d\t%d\t.\t+\t.\tgene_id \"lncRNA_%03d\"; transcript_id \"lncRNA_%03d\";\n", $1, new_start, $3, count, count
}' Pmea_lncRNA.bed > Pmea_lncRNA_for_counts.gtf
```

Run feature counts to quantify lncRNAs

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/ptua

module load conda/latest # need to load before making any conda envs
conda activate featurecounts

featureCounts \
  -p \
  -a Pmea_lncRNA_for_counts.gtf \
  -t lncRNA \
  -g gene_id \
  -o Pmea_lncRNA_counts_matrix.txt \
  *.bam
```

Submitted batch job 53420992

See how many lncRNAs are present with at least 5 counts in 3 biological reps 

```
awk -F'\t' '
  # 1. Print the header lines (anything starting with # or Geneid)
  /^#/ || /^Geneid/ { print; next } 
  
  # 2. For data lines, count how many samples have >= 5 counts
  {
    count = 0;
    for (i = 7; i <= 11; i++) {
      if ($i >= 5) count++;
    }
    
    # 3. If at least 3 samples met the criteria, print the line
    if (count >= 3) print;
  }' Pmea_lncRNA_counts_matrix.txt > filtered_Pmea_lncRNA_counts_matrix.txt
  
wc -l filtered_Pmea_lncRNA_counts_matrix.txt
6918 filtered_Pmea_lncRNA_counts_matrix.txt
```

Went from ~12000 to 6918 lncRNAs from expression filtering alone.

Get headers from filtered data and pull seqs from fasta

```
# 1. Create a list of the exact header strings that passed the filter
awk 'NR > 2 {print "transcript::" $2 ":" $3 "-" $4}' filtered_Pmea_lncRNA_counts_matrix.txt > pmea_coords_to_keep.txt

# Use seqtk to pull seqs 
seqtk subseq Pmea_lncRNA.fasta pmea_coords_to_keep.txt > exp_filtered_Pmea_lncRNAs.fasta
```

Run PLEK on new fasta file. `nano pmea_plek.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/ptua

# Load conda
module load conda/latest
# Use the full path to the environment to ensure it finds it
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /work/pi_hputnam_uri_edu/conda/envs/PLEK.1.2

# Define the PLEK directory
PLEK_DIR="/work/pi_hputnam_uri_edu/conda/envs/PLEK.1.2"

# Run PLEK using the full path to the script
# We use -thread 8 to match the cpus-per-task requested above
python ${PLEK_DIR}/PLEK.py \
    -fasta exp_filtered_Pmea_lncRNAs.fasta \
    -out pmea_plek_results.txt \
    -thread 8
```

Submitted batch job 53546179. In the log file: `Coding: 1987/6916=28.73048004626952%, Non-coding: 4929/6916=71.26951995373048%`. 

Extract the non-coding IDs from the PLEK results 

```
awk -F'\t' '$1 == "Non-coding" {print $3}' pmea_plek_results.txt | sed 's/^>//' > pmea_plek_non_coding_ids.txt
```

Filter fasta based on non-coding IDs

```
module load uri/main all/seqtk/1.4-GCC-12.3.0
seqtk subseq exp_filtered_Pmea_lncRNAs.fasta pmea_plek_non_coding_ids.txt > plek_exp_filtered_Pmea_lncRNAs.fasta
grep -c ">" plek_exp_filtered_Pmea_lncRNAs.fasta
4929
```

Run hmmscan to check for protein coding regions. `nano pmea_hmmscan.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/ptua

module load uri/main bio/HMMER/3.3.2-gompi-2022a

echo "Hmmscan commencing" $(date)
#echo "Unzip pfam file" $(date)

#gunzip ../Pfam-A.hmm.gz

#echo "press pfam db" $(date)
#hmmpress ../Pfam-A.hmm

echo "Run hmmscan" $(date)
hmmscan --domtblout pmea_hmmscan_output.domtblout ../Pfam-A.hmm plek_exp_filtered_Pmea_lncRNAs.fasta > pmea_hmmscan_output.txt

echo "hmmscan complete!" $(date)
```

Submitted batch job 53547214

Extract IDs that have protein domain hits

```
# Only remove hits with E-value < 0.001
grep -v '^#' pmea_hmmscan_output.domtblout | awk '$7 < 0.001 {print $4}' | sort | uniq > pmea_pfam_ids_to_remove.txt
```

Filter the fasta 

```
python3 -c "
import sys
# Load IDs to exclude into a set
with open('pmea_pfam_ids_to_remove.txt') as f:
    exclude = set(line.strip() for line in f if line.strip())
with open('plek_exp_filtered_Pmea_lncRNAs.fasta') as fasta:
    keep = True
    for line in fasta:
        if line.startswith('>'):
            # Check if the header ID is in our exclusion list
            header_id = line[1:].strip().split()[0]
            keep = header_id not in exclude
        if keep:
            print(line, end='')
" > pmea_final_lncRNAs.fasta
grep -c ">" pmea_final_lncRNAs.fasta
4924
```

Convert fasta headers to bed 

```
grep ">" pmea_final_lncRNAs.fasta | sed 's/>transcript:://' | awk -F'[: -]' '{print $1 "\t" $2 "\t" $3 "\t" "transcript::" $1 ":" $2 "-" $3}' > pmea_final_lncRNAs.bed
```

Sort the lncRNA bed file 

```
sort -k1,1 -k2,2n -k3,3n pmea_final_lncRNAs.bed > pmea_final_lncRNAs_sorted.bed
```

Sort gff file and filter for mRNA only 

```
sort -k1,1 -k4,4n Pocillopora_meandrina_HIv1.genes.gff3 > Pocillopora_meandrina_HIv1.genes_sorted.gff3

awk '$3 == "transcript"' Pocillopora_meandrina_HIv1.genes_sorted.gff3 > Pocillopora_meandrina_HIv1.genes_sorted_mRNA_only.gff3
```

Run bed closest 

```
module load bedtools2/2.31.1
bedtools closest -a pmea_final_lncRNAs_sorted.bed -b Pocillopora_meandrina_HIv1.genes_sorted_mRNA_only.gff3 > Pmea_lncRNA_closest_mRNA_only.bed
```

Turn lnc final bed into gff 

```
awk 'BEGIN{OFS="\t"} 
{
    count++;
    # Logic: if start ($2) is negative, force GTF start to 1.
    # Otherwise, use standard BED-to-GTF conversion ($2 + 1).
    if ($2 < 0) {
        new_start = 1;
    } else {
        new_start = $2 + 1;
    }

    # Safety: Ensure end ($3) is always at least 1 and greater than start
    if ($3 <= new_start) { $3 = new_start + 1; }

    printf "%s\t.\tlncRNA\t%d\t%d\t.\t+\t.\tgene_id \"lncRNA_%03d\"; transcript_id \"lncRNA_%03d\";\n", $1, new_start, $3, count, count
}' pmea_final_lncRNAs.bed > pmea_final_lncRNA_for_counts.gtf
```

Run feature counts to quantify final lncRNAs. `nano pmea_final_counts.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/ptua

module load conda/latest # need to load before making any conda envs
conda activate featurecounts

featureCounts \
  -p \
  -a pmea_final_lncRNA_for_counts.gtf \
  -t lncRNA \
  -g gene_id \
  -o Pmea_final_lncRNA_counts_matrix.txt \
  *.bam
```

Submitted batch job 53608810

Calculate number of lncRNAs per sample

```
awk 'NR > 2 { 
    for (i=7; i<=NF; i++) 
        if ($i > 0) count[i]++ 
    } 
    END { 
        for (i=7; i<=NF; i++) 
            print "Sample " (i-6) ": " count[i] " lncRNAs" 
    }' Pmea_final_lncRNA_counts_matrix.txt

Sample 1: 3612 lncRNAs
Sample 2: 4583 lncRNAs
Sample 3: 4672 lncRNAs
Sample 4: 4672 lncRNAs
Sample 5: 3516 lncRNAs
```

### Peve 

Turn bed into gff 

```
awk 'BEGIN{OFS="\t"} {
  count++; 
  printf "%s\t.\tlncRNA\t%d\t%d\t.\t+\t.\tgene_id \"lncRNA_%03d\"; transcript_id \"lncRNA_%03d\";\n", $1, $2+1, $3, count, count
}' Pevermanni_lncRNA.bed > Peve_lncRNA_for_counts.gtf
```

Run feature counts to quantify lncRNAs
```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/peve

module load conda/latest # need to load before making any conda envs
conda activate featurecounts

featureCounts \
  -p \
  -a Peve_lncRNA_for_counts.gtf \
  -t lncRNA \
  -g gene_id \
  -o lncRNA_counts_matrix.txt \
  *.sorted.bam
```

Submitted batch job 53387948

See how many lncRNAs are present with at least 5 counts in 3 biological reps 

```
awk -F'\t' '
  # 1. Print the header lines (anything starting with # or Geneid)
  /^#/ || /^Geneid/ { print; next } 
  
  # 2. For data lines, count how many samples have >= 5 counts
  {
    count = 0;
    for (i = 7; i <= 11; i++) {
      if ($i >= 5) count++;
    }
    
    # 3. If at least 3 samples met the criteria, print the line
    if (count >= 3) print;
  }' lncRNA_counts_matrix.txt > filtered_Peve_lncRNA_counts_matrix.txt
  
wc -l filtered_Peve_lncRNA_counts_matrix.txt
4294 filtered_Peve_lncRNA_counts_matrix.txt
```

Went from ~7000 to 4294 lncRNAs from expression filtering alone.

Get headers from filtered data and pull seqs from fasta

```
# Create a list of the exact header strings that passed the filter
awk 'NR > 2 {print "transcript::" $2 ":" ($3-1) "-" $4}' filtered_Peve_lncRNA_counts_matrix.txt > peve_coords_to_keep.txt

# Use seqtk to pull seqs 
seqtk subseq Peve_lncRNA.fasta peve_coords_to_keep.txt > exp_filtered_Peve_lncRNAs.fasta
```

Run PLEK on new fasta file. `nano peve_plek.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/peve

# Load conda
module load conda/latest
# Use the full path to the environment to ensure it finds it
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /work/pi_hputnam_uri_edu/conda/envs/PLEK.1.2

# Define the PLEK directory
PLEK_DIR="/work/pi_hputnam_uri_edu/conda/envs/PLEK.1.2"

# Run PLEK using the full path to the script
# We use -thread 8 to match the cpus-per-task requested above
python ${PLEK_DIR}/PLEK.py \
    -fasta exp_filtered_Peve_lncRNAs.fasta \
    -out peve_plek_results.txt \
    -thread 8
```

Submitted batch job 53549159. From the log file: `Coding: 1678/4412=38.03263825929284%, Non-coding: 2734/4412=61.96736174070716%`. 

Extract the non-coding IDs from the PLEK results 

```
awk -F'\t' '$1 == "Non-coding" {print $3}' peve_plek_results.txt | sed 's/^>//' > peve_plek_non_coding_ids.txt
```

Filter fasta based on non-coding IDs

```
module load uri/main all/seqtk/1.4-GCC-12.3.0
seqtk subseq exp_filtered_Peve_lncRNAs.fasta peve_plek_non_coding_ids.txt > plek_exp_filtered_Peve_lncRNAs.fasta
grep -c ">" plek_exp_filtered_Peve_lncRNAs.fasta
2814
```

Run hmmscan to check for protein coding regions. `nano peve_hmmscan.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/peve

module load uri/main bio/HMMER/3.3.2-gompi-2022a

echo "Hmmscan commencing" $(date)
#echo "Unzip pfam file" $(date)

#gunzip ../Pfam-A.hmm.gz

#echo "press pfam db" $(date)
#hmmpress ../Pfam-A.hmm

echo "Run hmmscan" $(date)
hmmscan --domtblout peve_hmmscan_output.domtblout ../Pfam-A.hmm plek_exp_filtered_Peve_lncRNAs.fasta > peve_hmmscan_output.txt

echo "hmmscan complete!" $(date)
```

Submitted batch job 53549217

Extract IDs that have protein domain hits

```
# Only remove hits with E-value < 0.001
grep -v '^#' peve_hmmscan_output.domtblout | awk '$7 < 0.001 {print $4}' | sort | uniq > peve_pfam_ids_to_remove.txt
```

Filter the fasta 

```
python3 -c "
import sys
# Load IDs to exclude into a set
with open('peve_pfam_ids_to_remove.txt') as f:
    exclude = set(line.strip() for line in f if line.strip())
with open('plek_exp_filtered_Peve_lncRNAs.fasta') as fasta:
    keep = True
    for line in fasta:
        if line.startswith('>'):
            # Check if the header ID is in our exclusion list
            header_id = line[1:].strip().split()[0]
            keep = header_id not in exclude
        if keep:
            print(line, end='')
" > peve_final_lncRNAs.fasta
grep -c ">" peve_final_lncRNAs.fasta
2812
```

Convert fasta headers to bed 

```
grep ">" peve_final_lncRNAs.fasta | sed 's/>transcript:://' | awk -F'[: -]' '{print $1 "\t" $2 "\t" $3 "\t" "transcript::" $1 ":" $2 "-" $3}' > peve_final_lncRNAs.bed
```

Sort the lncRNA bed file 

```
sort -k1,1 -k2,2n -k3,3n peve_final_lncRNAs.bed > peve_final_lncRNAs_sorted.bed
```

Sort gff file and filter for mRNA only 

```
sort -k1,1 -k4,4n Porites_evermanni_v1.annot.gff > Porites_evermanni_v1.annot_sorted.gff

awk '$3 == "mRNA"' Porites_evermanni_v1.annot_sorted.gff > Porites_evermanni_v1.annot_sorted_mRNA_only.gff
```

Run bed closest 

```
module load bedtools2/2.31.1
bedtools closest -a peve_final_lncRNAs_sorted.bed -b Porites_evermanni_v1.annot_sorted_mRNA_only.gff > Peve_lncRNA_closest_mRNA_only.bed
```

Turn lnc final bed into gff 

```
awk 'BEGIN{OFS="\t"} 
{
    count++;
    # Logic: if start ($2) is negative, force GTF start to 1.
    # Otherwise, use standard BED-to-GTF conversion ($2 + 1).
    if ($2 < 0) {
        new_start = 1;
    } else {
        new_start = $2 + 1;
    }

    # Safety: Ensure end ($3) is always at least 1 and greater than start
    if ($3 <= new_start) { $3 = new_start + 1; }

    printf "%s\t.\tlncRNA\t%d\t%d\t.\t+\t.\tgene_id \"lncRNA_%03d\"; transcript_id \"lncRNA_%03d\";\n", $1, new_start, $3, count, count
}' peve_final_lncRNAs.bed > peve_final_lncRNA_for_counts.gtf
```

Run feature counts to quantify final lncRNAs. `nano peve_final_counts.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive/peve

module load conda/latest # need to load before making any conda envs
conda activate featurecounts

featureCounts \
  -p \
  -a peve_final_lncRNA_for_counts.gtf \
  -t lncRNA \
  -g gene_id \
  -o Peve_final_lncRNA_counts_matrix.txt \
  *.bam
```

Submitted batch job 53608817.

Calculate number of lncRNAs per sample

```
awk 'NR > 2 { 
    for (i=7; i<=NF; i++) 
        if ($i > 0) count[i]++ 
    } 
    END { 
        for (i=7; i<=NF; i++) 
            print "Sample " (i-6) ": " count[i] " lncRNAs" 
    }' Peve_final_lncRNA_counts_matrix.txt
    
Sample 1: 2528 lncRNAs
Sample 2: 2261 lncRNAs
Sample 3: 2503 lncRNAs
Sample 4: 2353 lncRNAs
Sample 5: 2562 lncRNAs
```

## Compare lncRNA sequences across spp 

via blast 

`nano blash_lnc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-e5_deepdive

module load uri/main all/BLAST+/2.13.0-gompi-2022a

echo "make blast dbs" $(date)
makeblastdb -in apul/apul_final_lncRNAs.fasta -dbtype nucl -out Apul_lncRNA
makeblastdb -in peve/peve_final_lncRNAs.fasta -dbtype nucl -out Peve_lncRNA
makeblastdb -in ptua/pmea_final_lncRNAs.fasta -dbtype nucl -out Pmea_lncRNA

echo "Apul v Peve" $(date)
blastn -query apul/apul_final_lncRNAs.fasta \
       -db Peve_lncRNA \
       -outfmt "6 qseqid sseqid pident length evalue" \
       -evalue 1e-40 \
       -max_target_seqs 1 \
       -out apul_vs_peve.tab
       
echo "Apul v Pmea" $(date)
blastn -query apul/apul_final_lncRNAs.fasta \
       -db Pmea_lncRNA \
       -outfmt "6 qseqid sseqid pident length evalue" \
       -evalue 1e-40 \
       -max_target_seqs 1 \
       -out apul_vs_pmea.tab
       
echo "Peve v Pmea" $(date)
blastn -query peve/peve_final_lncRNAs.fasta \
       -db Pmea_lncRNA \
       -outfmt "6 qseqid sseqid pident length evalue" \
       -evalue 1e-40 \
       -max_target_seqs 1 \
       -out peve_vs_pmea.tab
```

Submitted batch job 53581586

Look at number of lncRNAs 

```
wc -l apul_vs_peve.tab
159 apul_vs_peve.tab

wc -l apul_vs_pmea.tab
111 apul_vs_pmea.tab

wc -l apul_vs_peve.tab
159 apul_vs_peve.tab
```