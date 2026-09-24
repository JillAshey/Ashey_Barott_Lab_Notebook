Gene expression analysis on POC Heron samples. Using Zoe's [code](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/Heron-Pdam-gene-expression.md) from Brown et al. 2025 as a reference since that paper used the same Pdam samples and similar sequencing methods. Analysis done on Unity. 

### Make a new scratch workspace and organize files 

```
ws_allocate POC_Heron 30

/scratch4/workspace/jillashey_uri_edu-POC_Heron
remaining extensions  : 10000
remaining time in days: 30
```

In new directory, make the following: 

```
cd /scratch4/workspace/jillashey_uri_edu-POC_Heron
mkdir data output scripts refs
cd output
mkdir QC alignment assembly 
cd QC 
mkdir raw trim
cd ../../data
mkdir raw trim 
```

Transferred fastq files into `raw` directory via globus. 

### QC raw data with fastqc and multiqc

Versions: FastQC v0.12.1; MultiQC v1.12

`nano raw_qc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 48:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules 
module load fastqc/0.12.1
module load uri/main MultiQC/1.12-foss-2021b

# Define directory paths
DATA_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/raw"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/QC/raw"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Initial qc of raw sequencing data" $(date)
fastqc "${DATA_DIR}"/*.fastq.gz -o "${OUTPUT_DIR}"
multiqc "${OUTPUT_DIR}" -o "${OUTPUT_DIR}"

echo "Initial qc of raw sequencing data complete and saved to ${OUTPUT_DIR}" $(date)
echo "Count number of raw reads per sample" $(date)
zgrep -c "^+$" "${DATA_DIR}"/*.fastq.gz > "${DATA_DIR}"/raw_read_counts.txt
echo "Read count complete" $(date)
```

Submitted batch job 64455548. Data looks pretty good overall, but high sequence duplication levels across the board. These samples ranged from 64.6-90.7% sequence duplication levels, which are high compared to Brown et al. 2025 [raw sample QC](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/data/raw_qc/raw_qc_multiqc_report.html).  

### Trim data with fastp and re-QC with fastqc and multiqc 

Versions: FastQC v0.12.1; MultiQC v1.12; fastp 0.23.4

`nano trim_qc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules 
module load fastqc/0.12.1
module load uri/main MultiQC/1.12-foss-2021b
module load fastp/0.23.4

# Define directory paths
RAW_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/raw"
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/QC/trim"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}" "${TRIM_DATA}"

echo "Begin the trim" $(date)
for i in "${RAW_DATA}"/*.fastq.gz; do
    fname=$(basename "$i")    
    trimmed_file="${TRIM_DATA}/trim.${fname}"
    fastp --in1 "${i}" \
          --out1 "${trimmed_file}" \
          --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
          --trim_poly_x 6 \
          -q 30 \
          -y \
          -Y 50 
    fastqc "${trimmed_file}" -o "${OUTPUT_DIR}"
done

multiqc "${OUTPUT_DIR}" -o "${OUTPUT_DIR}"
echo "QC of trimmed sequencing data complete and saved to ${OUTPUT_DIR}" $(date)
echo "Count number of trimmed reads per sample" $(date)
zgrep -c "^+$" "${TRIM_DATA}"/*.fastq.gz > "${TRIM_DATA}"/trim_read_counts.txt
echo "Read count complete" $(date)
```

Submitted batch job 64490004. Still a lot of duplication but moving on to alignment. 

### Download genomic resources 

Because the Pocillopora genus is cryptic, I'm going to align these sequences with all available POC genomes. These samples are likely Pdam but the Pdam genome isn't actually Pdam anyway. 

- P. meandrina 
	- Location: Hawai'i
	- Citation: Stephens et al. 2022 
	- Accession: http://cyanophora.rutgers.edu/Pocillopora_meandrina/
	- Not formally validated
- P. tuahiniensis
	- Location: Mo'orea 
	- Citation: Ashey et al. 2026
	- Accession: PRJNA1496434
	- Validated in Huffmyer et al. 2026
- P. acuta
	- Location: Hawai'i
	- Citation: Stephens et al. 2022
	- Accession: http://cyanophora.rutgers.edu/Pocillopora_acuta/
	- Not formally validated, but individual was a triploid 
- P. damicornis 
	- Location: Panama 
	- Citation: Cunning et al. 2018
	- Accession: PRJNA454489
	- Determined to be P. grandis (Oury et al. 2023)
- P. verrucosa
	- Location: Saudi Arabia
	- Citation: Buitrago-Lopez et al. 2020
	- Accession: http://pver.reefgenomics.org/
	- Determined to be P. favosa (Oury et al. 2025)
- P. effusa 
	- Location: Mo'orea 
	- Citation: Noel et al. 2023
	- Accession: https://www.genoscope.cns.fr/corals/genomes.html
	- Validated in Voolstra et al. 2023, but lacks a formal taxonomic description

```
cd /scratch4/workspace/jillashey_uri_edu-POC_Heron/refs

# Pmea
wget http://cyanophora.rutgers.edu/Pocillopora_meandrina/Pocillopora_meandrina_HIv1.assembly.fasta.gz

# Ptua 
### location on unity 
ln -s /scratch4/workspace/jillashey_uri_edu-Ptua_genome/ptua_softmasked/Pocillopora_tuahiniensis_genome_v1.0.fasta.masked .

# Pacuta
wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.assembly.fasta.gz

# Pdam
wget http://pdam.reefgenomics.org/download/pdam_scaffolds.fasta.gz

# Pdam from NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_genomic.fna.gz

# Pver 
wget http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.fasta.gz

# Peff
wget https://www.genoscope.cns.fr/corals/data/Pocillopora_effusa_v3.fa

# Unzip files 
gunzip *
```

### Align reads to genome with hisat2 

Versions: Hisat2 v2.2.1; 

Align against all published POC genomes. 

`nano align_pmea.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Alignment to Pmea genome" $(date)

# Load modules 
module load uri/main all/HISAT2/2.2.1-gompi-2022a
module load samtools/1.19.2

# Define directory paths
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Pmea"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Index Pmea reference genome" $(date)
if [ ! -f "${REF_DIR}/Pmea_ref.1.ht2" ]; then
    echo "Indexing Pmea reference genome $(date)"
    hisat2-build -f "${REF_DIR}/Pocillopora_meandrina_HIv1.assembly.fasta" "${REF_DIR}/Pmea_ref"
else
    echo "Reference index already exists. Skipping build step."
fi

echo "Reference genome indexed, begin alignment" $(date)
for i in "${TRIM_DATA}"/*.fastq.gz; do
    fname=$(basename "$i") 
    sample_name="${fname%.fastq.gz}"
    echo "Aligning ${sample_name}..."
    hisat2 -p 8 --dta -x "${REF_DIR}/Pmea_ref" -U "${i}" | \
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${sample_name}.bam" -
    samtools index -@ 8 "${OUTPUT_DIR}/${sample_name}.bam"
    echo "${sample_name} aligned, sorted, and indexed!"
done

echo "Alignment complete, calculate mapping percentages" $(date)
STATS_FILE="${OUTPUT_DIR}/alignment_Pmea_summary.txt"

for i in "${OUTPUT_DIR}"/*.bam; do
    sample=$(basename "$i")
    echo "=== Sample: ${sample} ===" >> "${STATS_FILE}"
    samtools flagstat "${i}" | grep "mapped (" >> "${STATS_FILE}"
    echo "" >> "${STATS_FILE}"
done

echo "Summary report saved to ${STATS_FILE}"
```

Submitted batch job 64492344

`nano align_ptua.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Alignment to Ptua genome" $(date)

# Load modules 
module load uri/main all/HISAT2/2.2.1-gompi-2022a
module load samtools/1.19.2

# Define directory paths
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Ptua"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Index Ptua reference genome" $(date)
if [ ! -f "${REF_DIR}/Ptua_ref.1.ht2" ]; then
    echo "Indexing Ptua reference genome $(date)"
    hisat2-build -f "${REF_DIR}/Pocillopora_tuahiniensis_genome_v1.0.fasta.masked" "${REF_DIR}/Ptua_ref"
else
    echo "Reference index already exists. Skipping build step."
fi

echo "Reference genome indexed, begin alignment" $(date)
for i in "${TRIM_DATA}"/*.fastq.gz; do
    fname=$(basename "$i") 
    sample_name="${fname%.fastq.gz}"
    echo "Aligning ${sample_name}..."
    hisat2 -p 8 --dta -x "${REF_DIR}/Ptua_ref" -U "${i}" | \
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${sample_name}.bam" -
    samtools index -@ 8 "${OUTPUT_DIR}/${sample_name}.bam"
    echo "${sample_name} aligned, sorted, and indexed!"
done

echo "Alignment complete, calculate mapping percentages" $(date)
STATS_FILE="${OUTPUT_DIR}/alignment_Ptua_summary.txt"

for i in "${OUTPUT_DIR}"/*.bam; do
    sample=$(basename "$i")
    echo "=== Sample: ${sample} ===" >> "${STATS_FILE}"
    samtools flagstat "${i}" | grep "mapped (" >> "${STATS_FILE}"
    echo "" >> "${STATS_FILE}"
done

echo "Summary report saved to ${STATS_FILE}"
```

Submitted batch job 64492477

`nano align_pacu.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Alignment to Pacu genome" $(date)

# Load modules 
module load uri/main all/HISAT2/2.2.1-gompi-2022a
module load samtools/1.19.2

# Define directory paths
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Pacu"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Index Pacu reference genome" $(date)
if [ ! -f "${REF_DIR}/Pacu_ref.1.ht2" ]; then
    echo "Indexing Pacu reference genome $(date)"
    hisat2-build -f "${REF_DIR}/Pocillopora_acuta_HIv2.assembly.fasta" "${REF_DIR}/Pacu_ref"
else
    echo "Reference index already exists. Skipping build step."
fi

echo "Reference genome indexed, begin alignment" $(date)
for i in "${TRIM_DATA}"/*.fastq.gz; do
    fname=$(basename "$i") 
    sample_name="${fname%.fastq.gz}"
    echo "Aligning ${sample_name}..."
    hisat2 -p 8 --dta -x "${REF_DIR}/Pacu_ref" -U "${i}" | \
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${sample_name}.bam" -
    samtools index -@ 8 "${OUTPUT_DIR}/${sample_name}.bam"
    echo "${sample_name} aligned, sorted, and indexed!"
done

echo "Alignment complete, calculate mapping percentages" $(date)
STATS_FILE="${OUTPUT_DIR}/alignment_Pacu_summary.txt"

for i in "${OUTPUT_DIR}"/*.bam; do
    sample=$(basename "$i")
    echo "=== Sample: ${sample} ===" >> "${STATS_FILE}"
    samtools flagstat "${i}" | grep "mapped (" >> "${STATS_FILE}"
    echo "" >> "${STATS_FILE}"
done

echo "Summary report saved to ${STATS_FILE}"
```

Submitted batch job 64492478

`nano align_pdam.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Alignment to Pdam genome (which is actually Pgrandis)" $(date)

# Load modules 
module load uri/main all/HISAT2/2.2.1-gompi-2022a
module load samtools/1.19.2

# Define directory paths
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Pdam"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Index Pdam reference genome" $(date)
if [ ! -f "${REF_DIR}/Pdam_ref.1.ht2" ]; then
    echo "Indexing Pdam reference genome $(date)"
    hisat2-build -f "${REF_DIR}/pdam_scaffolds.fasta" "${REF_DIR}/Pdam_ref"
else
    echo "Reference index already exists. Skipping build step."
fi

echo "Reference genome indexed, begin alignment" $(date)
for i in "${TRIM_DATA}"/*.fastq.gz; do
    fname=$(basename "$i") 
    sample_name="${fname%.fastq.gz}"
    echo "Aligning ${sample_name}..."
    hisat2 -p 8 --dta -x "${REF_DIR}/Pdam_ref" -U "${i}" | \
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${sample_name}.bam" -
    samtools index -@ 8 "${OUTPUT_DIR}/${sample_name}.bam"
    echo "${sample_name} aligned, sorted, and indexed!"
done

echo "Alignment complete, calculate mapping percentages" $(date)
STATS_FILE="${OUTPUT_DIR}/alignment_Pdam_summary.txt"

for i in "${OUTPUT_DIR}"/*.bam; do
    sample=$(basename "$i")
    echo "=== Sample: ${sample} ===" >> "${STATS_FILE}"
    samtools flagstat "${i}" | grep "mapped (" >> "${STATS_FILE}"
    echo "" >> "${STATS_FILE}"
done

echo "Summary report saved to ${STATS_FILE}"
```

Submitted batch job 64492487

`nano align_pdam_ncbi.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Alignment to Pdam genome (which is actually Pgrandis)--NCBI version" $(date)

# Load modules 
module load uri/main all/HISAT2/2.2.1-gompi-2022a
module load samtools/1.19.2

# Define directory paths
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Pdam_NCBI"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Index Pdam NCBI reference genome" $(date)
if [ ! -f "${REF_DIR}/Pdam_ncbi_ref.1.ht2" ]; then
    echo "Indexing Pdam NCBI reference genome $(date)"
    hisat2-build -f "${REF_DIR}/GCF_003704095.1_ASM370409v1_genomic.fna" "${REF_DIR}/Pdam_ncbi_ref"
else
    echo "Reference index already exists. Skipping build step."
fi

echo "Reference genome indexed, begin alignment" $(date)
for i in "${TRIM_DATA}"/*.fastq.gz; do
    fname=$(basename "$i") 
    sample_name="${fname%.fastq.gz}"
    echo "Aligning ${sample_name}..."
    hisat2 -p 8 --dta -x "${REF_DIR}/Pdam_ncbi_ref" -U "${i}" | \
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${sample_name}.bam" -
    samtools index -@ 8 "${OUTPUT_DIR}/${sample_name}.bam"
    echo "${sample_name} aligned, sorted, and indexed!"
done

echo "Alignment complete, calculate mapping percentages" $(date)
STATS_FILE="${OUTPUT_DIR}/alignment_Pdam_ncbi_summary.txt"

for i in "${OUTPUT_DIR}"/*.bam; do
    sample=$(basename "$i")
    echo "=== Sample: ${sample} ===" >> "${STATS_FILE}"
    samtools flagstat "${i}" | grep "mapped (" >> "${STATS_FILE}"
    echo "" >> "${STATS_FILE}"
done

echo "Summary report saved to ${STATS_FILE}"
```
Submitted batch job 64494392

`nano align_pver.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Alignment to Pver genome (which is actually Pfavosa)" $(date)

# Load modules 
module load uri/main all/HISAT2/2.2.1-gompi-2022a
module load samtools/1.19.2

# Define directory paths
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Pver"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Index Pver reference genome" $(date)
if [ ! -f "${REF_DIR}/Pver_ref.1.ht2" ]; then
    echo "Indexing Pver reference genome $(date)"
    hisat2-build -f "${REF_DIR}/Pver_genome_assembly_v1.0.fasta" "${REF_DIR}/Pver_ref"
else
    echo "Reference index already exists. Skipping build step."
fi

echo "Reference genome indexed, begin alignment" $(date)
for i in "${TRIM_DATA}"/*.fastq.gz; do
    fname=$(basename "$i") 
    sample_name="${fname%.fastq.gz}"
    echo "Aligning ${sample_name}..."
    hisat2 -p 8 --dta -x "${REF_DIR}/Pver_ref" -U "${i}" | \
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${sample_name}.bam" -
    samtools index -@ 8 "${OUTPUT_DIR}/${sample_name}.bam"
    echo "${sample_name} aligned, sorted, and indexed!"
done

echo "Alignment complete, calculate mapping percentages" $(date)
STATS_FILE="${OUTPUT_DIR}/alignment_Pver_summary.txt"

for i in "${OUTPUT_DIR}"/*.bam; do
    sample=$(basename "$i")
    echo "=== Sample: ${sample} ===" >> "${STATS_FILE}"
    samtools flagstat "${i}" | grep "mapped (" >> "${STATS_FILE}"
    echo "" >> "${STATS_FILE}"
done

echo "Summary report saved to ${STATS_FILE}"
```

Submitted batch job 64492489

`nano align_peff.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Alignment to Peff genome" $(date)

# Load modules 
module load uri/main all/HISAT2/2.2.1-gompi-2022a
module load samtools/1.19.2

# Define directory paths
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Peff"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Index Peff reference genome" $(date)
if [ ! -f "${REF_DIR}/Peff_ref.1.ht2" ]; then
    echo "Indexing Peff reference genome $(date)"
    hisat2-build -f "${REF_DIR}/Pocillopora_effusa_v3.fa" "${REF_DIR}/Peff_ref"
else
    echo "Reference index already exists. Skipping build step."
fi

echo "Reference genome indexed, begin alignment" $(date)
for i in "${TRIM_DATA}"/*.fastq.gz; do
    fname=$(basename "$i") 
    sample_name="${fname%.fastq.gz}"
    echo "Aligning ${sample_name}..."
    hisat2 -p 8 --dta -x "${REF_DIR}/Peff_ref" -U "${i}" | \
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${sample_name}.bam" -
    samtools index -@ 8 "${OUTPUT_DIR}/${sample_name}.bam"
    echo "${sample_name} aligned, sorted, and indexed!"
done

echo "Alignment complete, calculate mapping percentages" $(date)
STATS_FILE="${OUTPUT_DIR}/alignment_Peff_summary.txt"

for i in "${OUTPUT_DIR}"/*.bam; do
    sample=$(basename "$i")
    echo "=== Sample: ${sample} ===" >> "${STATS_FILE}"
    samtools flagstat "${i}" | grep "mapped (" >> "${STATS_FILE}"
    echo "" >> "${STATS_FILE}"
done

echo "Summary report saved to ${STATS_FILE}"
```

Submitted batch job 64492497

SUMMARY

|         | Pacuta | Pdamicornis | Pdamicornis (NCBI) | Peffusa | Pmeandrina | Ptuahiniensis | Pverrucosa |
| ------- | ------ | ----------- | ------------------ | ------- | ---------- | ------------- | ---------- |
| AVERAGE | 57.09  | 9.15        | 12.54              | 71.71   | 57.94      | 76.82         | 38.19      |

This is surprising to me! I expected Pacuta to have the highest alignment, given that it is most closely related to Pdam (which these samples putatively are). But the samples aligned best against the Ptua genome and the Peffusa genome.

Marcelina plotted the alignment data: 

![](https://github.com/JillAshey/Ashey_Barott_Lab_Notebook/blob/main/images/alignment_timepoint.png?raw=true)

![](https://github.com/JillAshey/Ashey_Barott_Lab_Notebook/blob/main/images/alignment_site.png?raw=true)

Very interesting as well. For the shallow lagoon (SL), it looks like there are two groupings of samples depending on which genome you are looking at. 

Given the alignments (and the fact that our samples are putatively Pdam), I'm going to make gene count matrices for Pacuta and Ptuahiniensis. Should probably use the Ptuahiniensis one, given the higher mapping. 

### Align to symbiont genome and rRNAs 

High levels of duplication and multi-mapping in hisat2 is sus. Going to align to rRNA database and symbiont genome to see if there is contamination. 

Obtain symbiont genome and rRNA database. There is no genome for Cladocopium latusorum, which is the dominant symbiont in the Heron samples. I'm going to use Cladocopium goreaui, a closely related symbiont species. 

```
cd /scratch4/workspace/jillashey_uri_edu-POC_Heron/refs

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/947/184/155/GCA_947184155.2_Cgoreaui_SCF055-01_v2.1/GCA_947184155.2_Cgoreaui_SCF055-01_v2.1_genomic.fna.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_tax_silva.fasta.gz

gunzip *
```

`nano align_sym.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Alignment to sym genome" $(date)

# Load modules 
module load uri/main all/HISAT2/2.2.1-gompi-2022a
module load samtools/1.19.2

# Define directory paths
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Sym"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Index sym reference genome" $(date)
if [ ! -f "${REF_DIR}/Sym_ref.1.ht2" ]; then
    echo "Indexing Sym reference genome $(date)"
    hisat2-build -f "${REF_DIR}/GCA_947184155.2_Cgoreaui_SCF055-01_v2.1_genomic.fna" "${REF_DIR}/Sym_ref"
else
    echo "Reference index already exists. Skipping build step."
fi

echo "Reference genome indexed, begin alignment" $(date)
for i in "${TRIM_DATA}"/*.fastq.gz; do
    fname=$(basename "$i") 
    sample_name="${fname%.fastq.gz}"
    echo "Aligning ${sample_name}..."
    hisat2 -p 8 --dta -x "${REF_DIR}/Sym_ref" -U "${i}" | \
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${sample_name}.bam" -
    samtools index -@ 8 "${OUTPUT_DIR}/${sample_name}.bam"
    echo "${sample_name} aligned, sorted, and indexed!"
done

echo "Alignment complete, calculate mapping percentages" $(date)
STATS_FILE="${OUTPUT_DIR}/alignment_Sym_summary.txt"

for i in "${OUTPUT_DIR}"/*.bam; do
    sample=$(basename "$i")
    echo "=== Sample: ${sample} ===" >> "${STATS_FILE}"
    samtools flagstat "${i}" | grep "mapped (" >> "${STATS_FILE}"
    echo "" >> "${STATS_FILE}"
done

echo "Summary report saved to ${STATS_FILE}"
```

Submitted batch job 64514194. Sym alignments across all samples <2%, which is great! 

`nano align_rRNA.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB                
#SBATCH -t 72:00:00                
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Alignment to rRNA db" $(date)

# Load modules 
module load uri/main all/HISAT2/2.2.1-gompi-2022a
module load samtools/1.19.2

# Define directory paths
TRIM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/data/trim"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/rRNA"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Index rRNA reference genome" $(date)
if [ ! -f "${REF_DIR}/rRNA_ref.1.ht2" ]; then
    echo "Indexing rRNA db $(date)"
    hisat2-build -f "${REF_DIR}/SILVA_138_SSURef_tax_silva.fasta" "${REF_DIR}/rRNA_ref"
else
    echo "Reference index already exists. Skipping build step."
fi

echo "Reference genome indexed, begin alignment" $(date)
for i in "${TRIM_DATA}"/*.fastq.gz; do
    fname=$(basename "$i") 
    sample_name="${fname%.fastq.gz}"
    echo "Aligning ${sample_name}..."
    hisat2 -p 8 --dta -x "${REF_DIR}/rRNA_ref" -U "${i}" | \
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${sample_name}.bam" -
    samtools index -@ 8 "${OUTPUT_DIR}/${sample_name}.bam"
    echo "${sample_name} aligned, sorted, and indexed!"
done

echo "Alignment complete, calculate mapping percentages" $(date)
STATS_FILE="${OUTPUT_DIR}/alignment_rRNA_summary.txt"

for i in "${OUTPUT_DIR}"/*.bam; do
    sample=$(basename "$i")
    echo "=== Sample: ${sample} ===" >> "${STATS_FILE}"
    samtools flagstat "${i}" | grep "mapped (" >> "${STATS_FILE}"
    echo "" >> "${STATS_FILE}"
done

echo "Summary report saved to ${STATS_FILE}"
```

Submitted batch job 64514219. Less than <1% aligning to rRNAs which is also great! Pretty confident that these are all coral sequences. 

### Fix GFFs so they are compatible with stringtie 

#### Pacuta 

Zoe already fixed the Pacuta one using [this script](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/scripts/fix_gff_format.Rmd). I downloaded the gff [here](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/data/Pocillopora_acuta_HIv2.genes_fixed.gff3.gz). 

#### Ptuahiniensis

Fixed this in R on my local computer. 

```{r}
Load libraries
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(R.utils)

# Set working directory 
setwd("~/Desktop/GFFs")

# Read in gff
gff <- read.delim("Ptua_braker_cleaned.gff3", comment.char = "#", header = FALSE)

# Rename cols 
colnames(gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")

# Create transcript ID
gff$transcript_id <- sub(";.*", "", gff$gene)
gff$transcript_id <- gsub("ID=", "", gff$transcript_id) #remove ID= 
gff$transcript_id <- gsub("Parent=", "", gff$transcript_id) #remove Parent= 

# Create Parent ID 
gff$parent_id <- sub(".*Parent=", "", gff$gene)
gff$parent_id <- sub(";.*", "", gff$parent_id)
gff$parent_id <- gsub("ID=", "", gff$parent_id) #remove ID= 

# Add these back to gene column separated by semicolons
gff <- gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", gff$transcript_id, ";gene_id=", gff$parent_id),  paste0(gene)))

# Remove transcript and parent ID cols
gff<-gff %>%
  select(!transcript_id)%>%
  select(!parent_id)

# Save file and upload to unity to use in stringtie
write.table(gff, file="Ptua_fixed.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
```

### Assemble reads with stringtie

#### Pacu

`nano assemble_pacu.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Assembling samples with stringtie $(date)"

# Load modules
module load uri/main StringTie/2.2.1-GCC-11.2.0

# Define directory paths
BAM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Pacu"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/assembly/Pacu"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

for i in "${BAM_DATA}"/*.bam; do
    [ -e "$i" ] || continue  # Skip loop if no .bam files exist
    
    fname=$(basename "$i")
    sample_name="${fname%.bam}"
    
    echo "Assembling ${sample_name}..."
    stringtie -p 8 -e -B \
        -G "${REF_DIR}/Pocillopora_acuta_HIv2.genes_fixed.gff3" \
        -A "${OUTPUT_DIR}/${sample_name}.gene_abund.tab" \
        -o "${OUTPUT_DIR}/${sample_name}.gtf" \
        "${i}"
        
    echo "${sample_name} assembled!"
done

echo "Assembly complete! $(date)"
```

Submitted batch job 64823750

#### Ptua

`nano assemble_ptua.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

echo "Assembling samples with stringtie $(date)"

# Load modules
module load uri/main StringTie/2.2.1-GCC-11.2.0

# Define directory paths
BAM_DATA="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/alignment/Ptua"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/assembly/Ptua"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

for i in "${BAM_DATA}"/*.bam; do
    [ -e "$i" ] || continue  # Skip loop if no .bam files exist
    
    fname=$(basename "$i")
    sample_name="${fname%.bam}"
    
    echo "Assembling ${sample_name}..."
    stringtie -p 8 -e -B \
        -G "${REF_DIR}/Ptua_fixed.gff3" \
        -A "${OUTPUT_DIR}/${sample_name}.gene_abund.tab" \
        -o "${OUTPUT_DIR}/${sample_name}.gtf" \
        "${i}"
        
    echo "${sample_name} assembled!"
done

echo "Assembly complete! $(date)"
```

Submitted batch job 64833305

### Make gene count matrix 

Download the [prep_DE.py](https://github.com/gpertea/stringtie/blob/master/prepDE.py3) script from the stringtie github repo. I'm downloading the one compatible with python3 since that is what we have on the server. Once downloaded, make the file executable. 

```
chmod +x /scratch4/workspace/jillashey_uri_edu-POC_Heron/prepDE.py3
```

#### Pacu

`nano prepDE_pacu.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules
module load uri/main StringTie/2.2.1-GCC-11.2.0 GffCompare/0.12.6-GCC-11.2.0
module load python/3.9.19

# Define directory paths
GTF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/assembly/Pacu"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Make gtf list for stringtie merge
ls "${GTF_DIR}"/*.gtf > "${GTF_DIR}/gtf_list_pacu.txt"

echo "Merge gtfs with stringtie $(date)"
stringtie --merge -e -p 8 \
    -G "${REF_DIR}/Pocillopora_acuta_HIv2.genes_fixed.gff3" \
    -o "${GTF_DIR}/POC_Heron_Pacu_merged.gtf" \
    "${GTF_DIR}/gtf_list_pacu.txt"

echo "Merge complete, check completeness $(date)"
gffcompare -r "${REF_DIR}/Pocillopora_acuta_HIv2.genes_fixed.gff3" \
    -G \
    -o "${GTF_DIR}/merge" \
    "${GTF_DIR}/POC_Heron_Pacu_merged.gtf"

echo "Check complete, assemble count matrix $(date)"

# Create listGTF.txt mapping sample_id to file path for prepDE.py
> "${GTF_DIR}/listGTF.txt" # Clear file if exists
for f in "${GTF_DIR}"/*.gtf; do
    [ -e "$f" ] || continue
    # Exclude the merged output GTF
[[ "$f" == *"POC_Heron_Pacu_merged.gtf"* || "$f" == *"merge.annotated.gtf"* ]] && continue    
    sample_name=$(basename "$f" .gtf)
    echo "${sample_name} ${f}" >> "${GTF_DIR}/listGTF.txt"
done

# Run prepDE.py
/scratch4/workspace/jillashey_uri_edu-POC_Heron/prepDE.py3 -g "${GTF_DIR}/POC_Heron_Pacu_gene_counts.csv" -t "${GTF_DIR}/POC_Heron_Pacu_transcript_counts.csv" \
          -i "${GTF_DIR}/listGTF.txt"

echo "Counts matrix complete $(date)"
```

Submitted batch job 64831578. Look at how many genes are expressed (ie >1 in at least 1 sample) and how many are not expressed (ie 0 in all samples). 

```
cd /scratch4/workspace/jillashey_uri_edu-POC_Heron/output/assembly/Pacu

awk -F',' 'NR>1 {
    expressed = 0;
    for (i=2; i<=NF; i++) {
        if ($i > 0) { expressed = 1; break }
    }
    if (expressed) pass++; else zero++;
} END {
    print "Expressed (>0 in >=1 sample): " pass;
    print "Unexpressed (0 in all samples): " zero;
    print "Total genes:                   " pass + zero;
}' POC_Heron_Pacu_gene_counts.csv
Expressed (>0 in >=1 sample): 19279
Unexpressed (0 in all samples): 14451
Total genes:                   33730
```

#### Ptua

`nano prepDE_ptua.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules
module load uri/main StringTie/2.2.1-GCC-11.2.0 GffCompare/0.12.6-GCC-11.2.0
module load python/3.9.19

# Define directory paths
GTF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/output/assembly/Ptua"
REF_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs"

# Make gtf list for stringtie merge
ls "${GTF_DIR}"/*.gtf > "${GTF_DIR}/gtf_list_ptua.txt"

echo "Merge gtfs with stringtie $(date)"
stringtie --merge -e -p 8 \
    -G "${REF_DIR}/Ptua_fixed.gff3" \
    -o "${GTF_DIR}/POC_Heron_Ptua_merged.gtf" \
    "${GTF_DIR}/gtf_list_ptua.txt"

echo "Merge complete, check completeness $(date)"
gffcompare -r "${REF_DIR}/Ptua_fixed.gff3" \
    -G \
    -o "${GTF_DIR}/merge" \
    "${GTF_DIR}/POC_Heron_Ptua_merged.gtf"

echo "Check complete, assemble count matrix $(date)"

# Create listGTF.txt mapping sample_id to file path for prepDE.py
> "${GTF_DIR}/listGTF.txt" # Clear file if exists
for f in "${GTF_DIR}"/*.gtf; do
    [ -e "$f" ] || continue
    # Exclude the merged output GTF
[[ "$f" == *"POC_Heron_Ptua_merged.gtf"* || "$f" == *"merge.annotated.gtf"* ]] && continue    
    sample_name=$(basename "$f" .gtf)
    echo "${sample_name} ${f}" >> "${GTF_DIR}/listGTF.txt"
done

# Run prepDE.py
/scratch4/workspace/jillashey_uri_edu-POC_Heron/prepDE.py3 -g "${GTF_DIR}/POC_Heron_Ptua_gene_counts.csv" -t "${GTF_DIR}/POC_Heron_Ptua_transcript_counts.csv" \
          -i "${GTF_DIR}/listGTF.txt"

echo "Counts matrix complete $(date)"
```

Submitted batch job 64834700. Look at how many genes are expressed (ie >1 in at least 1 sample) and how many are not expressed (ie 0 in all samples). 

```
cd /scratch4/workspace/jillashey_uri_edu-POC_Heron/output/assembly/Ptua

awk -F',' 'NR>1 {
    expressed = 0;
    for (i=2; i<=NF; i++) {
        if ($i > 0) { expressed = 1; break }
    }
    if (expressed) pass++; else zero++;
} END {
    print "Expressed (>0 in >=1 sample): " pass;
    print "Unexpressed (0 in all samples): " zero;
    print "Total genes:                   " pass + zero;
}' POC_Heron_Ptua_gene_counts.csv
Expressed (>0 in >=1 sample): 16754
Unexpressed (0 in all samples): 10327
Total genes:                   27081
```


