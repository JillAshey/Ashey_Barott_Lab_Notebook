I am working from the lncRNA fasta and bed files generated from our original analysis (see [github](https://github.com/urol-e5/deep-dive/tree/main)) and running some more analyses for stringency. This will likely affect the counts and downstream results, so these will have to be updated as well. Updated worflow: 

- Filter based on expression (ie at least 5 counts for 3 replicates)
- Run PLEK -- maybe not if it installs??
- Compare against Pfam 
- Update downstream workflows 
	- Number of unique v shared lncRNAs 
	- Length plot 
	- Gene lncRNA overlap 

Install needed software for updated workflow

PLEK: 

```
# Download PLEK.1.2.tar.gz from https://sourceforge.net/projects/plek/files/
and decompress it.
tar zvxf PLEK.1.2.tar.gz
# Compile PLEK.
cd PLEK.1.2
python PLEK_setup.py
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
# (If you don't have seqtk, see the python script below)
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


okay workflow is this: 

- get 





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
  }' lncRNA_counts_matrix.txt > filtered_lncRNA_5counts_3reps.txt
  
 wc -l filtered_lncRNA_counts.txt
4294 filtered_lncRNA_counts.txt
```

See how many lncRNAs are present with at least 5 counts in two biological reps 

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
    if (count >= 2) print;
  }' lncRNA_counts_matrix.txt > filtered_lncRNA_5counts_2reps.txt
  
wc -l filtered_lncRNA_5counts_2reps.txt
5218 filtered_lncRNA_5counts_2reps.txt
```

See how many lncRNAs are present with at least 3 counts in three biological reps 

```
awk -F'\t' '
  # 1. Print the header lines (anything starting with # or Geneid)
  /^#/ || /^Geneid/ { print; next } 
  
  # 2. For data lines, count how many samples have >= 5 counts
  {
    count = 0;
    for (i = 7; i <= 11; i++) {
      if ($i >= 3) count++;
    }
    
    # 3. If at least 3 samples met the criteria, print the line
    if (count >= 3) print;
  }' lncRNA_counts_matrix.txt > filtered_lncRNA_3counts_3reps.txt
 
wc -l filtered_lncRNA_3counts_3reps.txt
4682 filtered_lncRNA_3counts_3reps.txt
```

See how many lncRNAs are present with at least 3 counts in two biological reps 

```
awk -F'\t' '
  # 1. Print the header lines (anything starting with # or Geneid)
  /^#/ || /^Geneid/ { print; next } 
  
  # 2. For data lines, count how many samples have >= 5 counts
  {
    count = 0;
    for (i = 7; i <= 11; i++) {
      if ($i >= 3) count++;
    }
    
    # 3. If at least 3 samples met the criteria, print the line
    if (count >= 3) print;
  }' lncRNA_counts_matrix.txt > filtered_lncRNA_2counts_2reps.txt
  
wc -l filtered_lncRNA_2counts_2reps.txt
4682 filtered_lncRNA_2counts_2reps.txt
```

Look at overall counts

```
awk -F'\t' '
  # 1. Capture the sample names from the header line
  /^Geneid/ {
    for (i = 7; i <= NF; i++) {
        sample_names[i] = $i
    }
    next
  }
  
  # 2. Skip the program info line
  /^#/ { next }

  # 3. For every data line, check each sample column
  {
    for (i = 7; i <= NF; i++) {
        if ($i >= 1) {
            count[i]++
        }
    }
  }

  # 4. Print the results at the end
  END {
    printf "%-30s\t%s\n", "Sample_Name", "lncRNAs_Detected"
    print "------------------------------------------------------------"
    for (i = 7; i <= 11; i++) {
        printf "%-30s\t%d\n", sample_names[i], count[i]
    }
  }' lncRNA_counts_matrix.txt
```