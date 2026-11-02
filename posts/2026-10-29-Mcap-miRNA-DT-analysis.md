Mcap miRNA developmental time series processing 

```
ws_allocate mcap_mirna 30
Info: could not read email from users config ~/.ws_user.conf.
Info: reminder email will be sent to local user account
Info: creating workspace.
/scratch3/workspace/jillashey_uri_edu-mcap_mirna
```

### QC 

`nano raw_qc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-mcap_mirna

# Load modules 
module load fastqc/0.12.1 uri/main all/MultiQC/1.12-foss-2021b

# Define paths 
FASTQ_DIR="/project/pi_hputnam_uri_edu/raw_sequencing_data/20250624_Mcap_2023_smRNA"
OUTPUT_DIR="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/fastqc_raw"

# Create output directory if doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Run fastqc + multiqc on data" $(date)
for file in "${FASTQ_DIR}"/*.fastq.gz; do 
    fastqc "$file" --outdir "${OUTPUT_DIR}" -t 2
done

multiqc "${OUTPUT_DIR}" -o "${OUTPUT_DIR}"

echo "Raw QC complete!" $(date)
``` 

Submitted batch job 66284673

Will be using [NEB adapter fasta](https://github.com/urol-e5/deep-dive/blob/3e47da33ae0ca52ca25b48c3c8fa70a51a76fde7/E-Peve/output/06-Peve-sRNAseq-trimming/NEB-adapters.fasta#L3) to remove adapters, since these samples were prepped with the NEB small RNA library prep kit. 

### Trim + QC 

`nano trim_qc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-mcap_mirna

# Load modules 
module load fastqc/0.12.1 uri/main all/fastp/0.23.2-GCC-11.2.0 all/MultiQC/1.12-foss-2021b

# Define paths 
FASTQ_DIR="/project/pi_hputnam_uri_edu/raw_sequencing_data/20250624_Mcap_2023_smRNA"
FASTQ_TRIM="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/trim"
FASTQ_TRIM_QC="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/fastqc_trimmed"

# Create output directory if doesn't exist
mkdir -p "${FASTQ_TRIM}"
mkdir -p "${FASTQ_TRIM_QC}"

# Populate R1 and R2 arrays matching specific naming pattern
fastq_array_R1=("${FASTQ_DIR}"/*_R1_001.fastq.gz)
fastq_array_R2=("${FASTQ_DIR}"/*_R2_001.fastq.gz)

echo "Start trimming" $(date)
echo 'export max_read_length=31'
threads=${SLURM_CPUS_ON_NODE:-2}

R1_names_array=()
for f in "${fastq_array_R1[@]}"; do
    R1_names_array+=("$(basename "$f")")
done
for index in "${!fastq_array_R1[@]}"
do
    R1_sample_name="${R1_names_array[index]%%_*}"
    trimmed_fastp_report_title="${R1_sample_name}-fastp-adapters-polyG-${max_read_length}bp-merged"
    merged_output="${FASTQ_TRIM}/${trimmed_fastp_report_title}.fq.gz"

    # Begin fastp trimming
    fastp \
        --in1 "${fastq_array_R1[index]}" \
        --in2 "${fastq_array_R2[index]}" \
        --adapter_fasta "${NEB_adapters_fasta}" \
        --trim_poly_g \
        --overlap_len_require 17 \
        --length_limit 31 \
        --merge \
        --merged_out "${merged_output}" \
        --thread "${threads}" \
        --html "${FASTQ_TRIM}/${trimmed_fastp_report_title}.html" \
        --report_title "${trimmed_fastp_report_title}"
        
    echo "Finished processing sample: ${R1_sample_name}"
done

echo "Run fastqc on trimmed/merged data" $(date)
for file in "${FASTQ_TRIM}"/*.fq.gz; do
    fastqc "$file" --outdir "${FASTQ_TRIM_QC}" -t "${threads}"
done

echo "Run multiqc on trimmed/merged QC reports" $(date)
multiqc "${FASTQ_TRIM_QC}" -o "${FASTQ_TRIM_QC}"

echo "Trimmed QC complete!" $(date)
``` 

Submitted batch job 66297171. Job stalled when it got to M7 because M7 doesn't have an R2 file??? not clue where it went. So rerunning the fastp script but writing it so that it checks for both R1 and R2 first and if it doesnt find R2 it processes that sample as a single end file. 

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-mcap_mirna

# Load modules 
module load fastqc/0.12.1 uri/main all/fastp/0.23.2-GCC-11.2.0 all/MultiQC/1.12-foss-2021b

# Define paths 
FASTQ_DIR="/project/pi_hputnam_uri_edu/raw_sequencing_data/20250624_Mcap_2023_smRNA"
FASTQ_TRIM="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/trim"
FASTQ_TRIM_QC="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/fastqc_trimmed"

# Create output directories if they don't exist
mkdir -p "${FASTQ_TRIM}"
mkdir -p "${FASTQ_TRIM_QC}"

max_read_length=31
threads=${SLURM_CPUS_ON_NODE:-2}

echo "Start trimming" $(date)

# Loop through all R1 files dynamically
for f1 in "${FASTQ_DIR}"/*_R1_001.fastq.gz; do
    basename_f1=$(basename "$f1")
    r1_sample_name="${basename_f1%%_*}"
    
    # Construct the expected R2 filename path
    f2="${FASTQ_DIR}/${r1_sample_name}_R2_001.fastq.gz"
    
    # Check if the R2 file exists
    if [ -f "$f2" ]; then
        # --- PAIRED-END PROCESSING ---
        trimmed_fastp_report_title="${r1_sample_name}-fastp-adapters-polyG-${max_read_length}bp-merged"
        merged_output="${FASTQ_TRIM}/${trimmed_fastp_report_title}.fq.gz"

        fastp \
            --in1 "$f1" \
            --in2 "$f2" \
            --adapter_fasta "${NEB_adapters_fasta}" \
            --trim_poly_g \
            --overlap_len_require 17 \
            --length_limit 31 \
            --merge \
            --merged_out "${merged_output}" \
            --thread "${threads}" \
            --html "${FASTQ_TRIM}/${trimmed_fastp_report_title}.html" \
            --report_title "${trimmed_fastp_report_title}"
            
    else
        # --- SINGLE-END PROCESSING (e.g., M7) ---
        echo "WARNING: R2 file missing for ${r1_sample_name}. Processing as Single-End."
        trimmed_fastp_report_title="${r1_sample_name}-fastp-adapters-polyG-${max_read_length}bp-merged"
        se_output="${FASTQ_TRIM}/${trimmed_fastp_report_title}.fq.gz"

        fastp \
            --in1 "$f1" \
            --out1 "${se_output}" \
            --adapter_fasta "${NEB_adapters_fasta}" \
            --trim_poly_g \
            --length_limit 31 \
            --thread "${threads}" \
            --html "${FASTQ_TRIM}/${trimmed_fastp_report_title}.html" \
            --report_title "${trimmed_fastp_report_title}"
    fi
        
    echo "Finished processing sample: ${r1_sample_name}"
done

echo "Run fastqc on trimmed/merged data" $(date)
for file in "${FASTQ_TRIM}"/*.fq.gz; do
    fastqc "$file" --outdir "${FASTQ_TRIM_QC}" -t "${threads}"
done

echo "Run multiqc on trimmed/merged QC reports" $(date)
multiqc "${FASTQ_TRIM_QC}" -o "${FASTQ_TRIM_QC}"

echo "Trimmed QC complete!" $(date)
```

Submitted batch job 66371738

### Shortstack - miRNA ID

I made a miRNA fasta file of the miRNAs previously IDed in other studies, as well as those from mirbase. There have been a few other studies recently that have came out so I'm going to add those to my fasta. In Table S3 from the preprint [Grinblat et al.](https://www.biorxiv.org/content/10.64898/2026.05.09.724056v1), there are some more I need to add (from Fridrich et al. 2020, Despard et al. 2026, and Grinbalt et al.). I'll use a python script to add them to the fasta. 

```
module load python/3.9.19

python3 - << 'EOF'
import csv
csv_file = "Grinbalt_et_al_2026_TableS3.csv"
fasta_file = "miRNA_reference.fasta"

# Target keywords to match citations flexibly (catches Grinblat/Grinbalt, Despard, and Fridrich)
target_citations = ["grinblat", "grinbalt", "despard", "fridrich"]

added_count = 0

with open(csv_file, mode="r", encoding="utf-8") as infile, open(
    fasta_file, mode="a", encoding="utf-8"
) as outfile:

    reader = csv.DictReader(infile)

    for row in reader:
        citation = row["Citation"].strip()
        cit_lower = citation.lower()

        # Check if the citation matches any of our targets
        if any(target in cit_lower for target in target_citations):
            miRNA_name = row["miRNA"].strip()
            species = row["Species"].strip()
            # Convert RNA 'U' to DNA 'T' for reference fasta compatibility
            seq = row["Mature_miRNA_sequence"].strip().upper().replace("U", "T")

            # Format the header with miRNA name, species, and citation
            header = f">{miRNA_name}::{species}::{citation}"

            # Write to fasta file
            outfile.write(f"{header}\n{seq}\n")
            added_count += 1

print(
    f"Successfully extracted and appended {added_count} miRNA sequences to {fasta_file}!")
EOF
```

Check number of miRNAs now in fasta

```
grep -c ">" miRNA_reference.fasta 
50011
```

Run shortstack to ID novel and known miRNAs. `nano shortstack.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-mcap_mirna

# Load modules
module load uri/main
module load conda/latest 
conda activate /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 
threads=${SLURM_CPUS_ON_NODE:-2}

# Define paths 
OUT_DIR="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/shortstack_out"
GENOME="/work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.assembly.fasta"

# Create output directories if they don't exist
mkdir -p "${OUT_DIR}"

echo "Running shortstack on trimmed reads for mcap" $(date)
ShortStack \
    --genomefile "${GENOME}" \
    --readfile trim/*.fq.gz \
    --known_miRNAs miRNA_reference.fasta \
    --outdir "${OUT_DIR}" \
    --dn_mirna \
    --threads ${threads}

echo "Shortstack complete!" $(date)
```

Submitted batch job 66382579




### mirdeep2

There are a few steps for mirdeep2. First, I need to map the trimmed reads to the genome with the mapper.pl script
I have 32 samples, need make a config file with each line designated a file location and 3 letter code. 

```
# Define the directory path
DIR="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/trim"

# Create or overwrite config.txt
> config.txt

# Loop through the matching files
for filepath in "$DIR"/M*-fastp-adapters-polyG-31bp-merged.fq; do
    # Skip if no matches found
    [ -e "$filepath" ] || continue
    
    # Get just the filename
    filename=$(basename "$filepath")
    
    # Extract the number after 'M'
    num=$(echo "$filename" | sed -E 's/^M([0-9]+)-.*/\1/')
    
    # Format the 3-letter code (e.g., s6 -> s06, s10 -> s10)
    code=$(printf "s%02d" "$num")
    
    # Append the file path and code to config.txt (tab-separated)
    echo -e "$filepath\t$code" >> config.txt
done
```

Okay should be good to go. `nano mapper_mirdeep2.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-mcap_mirna

# Load modules 
module load uri/main Bowtie/1.3.1-GCC-11.3.0

# Define paths 
OUT_DIR="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/mirdeep2"
GENOME="/work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.assembly.fasta"

# Create output directories if they don't exist
mkdir -p "${OUT_DIR}"

echo "Index Mcap genome" $(date)

# Index the reference genome for Mcap 
bowtie-build "${GENOME}" Mcap_ref.btindex

echo "Referece genome indexed!" $(date)
echo "Unload unneeded packages and run mapper script for trimmed reads" $(date)

# Load modules
module load Bowtie/1.3.1-GCC-11.3.0
module load conda/latest
conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 

mapper.pl config.txt -d -e -h -l 18 -v -m -p Mcap_ref.btindex -s "${OUT_DIR}"/mapped_reads.fa -t "${OUT_DIR}"/mapped_reads_vs_genome.arf

echo "Mapping complete!" $(date)
```

Submitted batch job 66510985

Could also add other arguments like -q (map wih one mismach in seed) (look at [docs](https://www.mdc-berlin.de/content/mirdeep2-documentation?mdcbl%5B0%5D=/n-rajewsky%23t-joinus&mdctl=0&mdcou=20738&mdcot=6&mdcbv=5yVgWgdVqnN3ZrjWVN8SYGoN7LptyyVi7YGzkgWdh9E))

CHECK OUTPUT



`nano id_mirdeep2.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-mcap_mirna

# Load modules 
conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 

# Define paths 
OUT_DIR="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/mirdeep2"
GENOME="/work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.assembly.fasta"

echo "ID miRNAs" $(date)

miRDeep2.pl "${OUT_DIR}"/mapped_reads.fa "${GENOME}" "${OUT_DIR}"/mapped_reads_vs_genome.arf miRNA_reference.fasta none none -P 2>report.log

echo "miRNA ID complete" $(date)

# code for cnidarian sperm paper 
#miRDeep2.pl ahya_1_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta ahya_1_mappings.arf /work/#pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P
```

Submitted batch job 66523625



do i need to clean the miRNA fasta reference??





TRY mirdeep2




For mirna ID

- https://www.mirbase.org/mirsubmit/


For target prediction

- also try https://github.com/carringtonlab/TargetFinder
- remove strict flag from miranda?