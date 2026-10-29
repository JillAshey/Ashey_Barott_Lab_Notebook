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

Submitted batch job 66297171

next: shortstack
- get grinbalt and despard miRNA seqs from here: https://docs.google.com/spreadsheets/d/1ykgRIpTRm3Pm6bR3ZxFTkObTqImGOyo3zVECqaMPNQo/edit?gid=1498714596#gid=1498714596


For mirna ID

- https://www.mirbase.org/mirsubmit/
- what fasta to use????


For target prediction

- also try https://github.com/carringtonlab/TargetFinder
- remove strict flag from miranda?