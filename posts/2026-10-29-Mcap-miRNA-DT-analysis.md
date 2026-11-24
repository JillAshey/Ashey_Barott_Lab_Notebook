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

Could also add other arguments like -q (map wih one mismach in seed) (look at [docs](https://www.mdc-berlin.de/content/mirdeep2-documentation?mdcbl%5B0%5D=/n-rajewsky%23t-joinus&mdctl=0&mdcou=20738&mdcot=6&mdcbv=5yVgWgdVqnN3ZrjWVN8SYGoN7LptyyVi7YGzkgWdh9E)).

Output in error file:

```
Log file for this run is in mapper_logs and called mapper.log_575335
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 543775401	251649719	292125682	46.278	53.722
s06: 14355963	6754655	7601308	47.051	52.949
s07: 13941025	6575355	7365670	47.166	52.834
s08: 31413403	14689649	16723754	46.762	53.238
s09: 15980019	7478829	8501190	46.801	53.199
s10: 14386193	6604329	7781864	45.907	54.093
s11: 13530010	6045859	7484151	44.685	55.315
s13: 18471008	9358351	9112657	50.665	49.335
s14: 16957559	8382152	8575407	49.430	50.570
s23: 19709125	9291651	10417474	47.144	52.856
s24: 14854873	6938496	7916377	46.709	53.291
s26: 22162895	10641446	11521449	48.015	51.985
s28: 19186633	8676968	10509665	45.224	54.776
s35: 16992990	7920318	9072672	46.609	53.391
s36: 19895258	9089137	10806121	45.685	54.315
s37: 18830303	8635828	10194475	45.861	54.139
s39: 20527294	9142625	11384669	44.539	55.461
s47: 18502998	9224076	9278922	49.852	50.148
s48: 11851570	5890693	5960877	49.704	50.296
s51: 13073558	6547662	6525896	50.083	49.917
s52: 16468283	7533372	8934911	45.745	54.255
s60: 15871202	7377522	8493680	46.484	53.516
s61: 17789304	8468605	9320699	47.605	52.395
s62: 15614266	7312845	8301421	46.834	53.166
s63: 15178415	7490228	7688187	49.348	50.652
s72: 18973747	8382931	10590816	44.182	55.818
s73: 16038291	7401173	8637118	46.147	53.853
s74: 14910351	6970109	7940242	46.747	53.253
s75: 13471057	5809978	7661079	43.129	56.871
s85: 16702307	6812227	9890080	40.786	59.214
s86: 19626774	8389260	11237514	42.744	57.256
s87: 12169060	4773551	7395509	39.227	60.773
s88: 16339667	7039839	9299828	43.084	56.916
```

Edit miRNA reference fasta so that there is no whitespace in the headers. 

```
awk '{if (/^>/) { gsub(/[[:space:]]+/, "-"); } print}' miRNA_reference.fasta > miRNA_reference_fixed.fasta
```

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

miRDeep2.pl "${OUT_DIR}"/mapped_reads.fa "${GENOME}" "${OUT_DIR}"/mapped_reads_vs_genome.arf miRNA_reference_fixed.fasta none none -P -g -1 2>report.log

echo "miRNA ID complete" $(date)

# code for cnidarian sperm paper 
#miRDeep2.pl ahya_1_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta ahya_1_mappings.arf /work/#pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P
```

Submitted batch job 66549321. This took almost 10 days wtf. 

Now I need to filter for putative miRNAs based on the following qualities from Fromm et al. 2015 and Praher et al. 2021. `nano parse_mcap_miRNAs.py`

```
#!/usr/bin/env python3
"""
parse_mcap_miRNAs.py
Extracts novel and known miRNAs from miRDeep2 CSV and MRD files, calculates
Praher et al. (2021) criteria (5' consistency, duplex complementarity, numeric 3' overhang in bp),
adds sequence length columns (bp), and outputs a master TSV ready for R downstream analysis.
"""

import sys
import re
from collections import Counter

if len(sys.argv) < 4:
    print("Usage: python3 parse_mcap_miRNAs.py <result_csv> <mrd_file> <output_tsv>")
    sys.exit(1)

csv_path = sys.argv[1]
mrd_path = sys.argv[2]
out_path = sys.argv[3]

# -------------------------------------------------------------------------
# 1. Helper Functions
# -------------------------------------------------------------------------

def calc_complementarity(mat, star):
    """Calculates maximum complementary base matches between mature and star strands."""
    if not mat or not star or mat == '-' or star == '-':
        return 0
    comp_map = {'a':'u', 'u':'a', 'c':'g', 'g':'c', 'A':'U', 'U':'A', 'C':'G', 'G':'C'}
    rev_star_comp = "".join([comp_map.get(b, 'N') for b in reversed(star)])
    
    max_matches = 0
    for shift in range(-len(rev_star_comp), len(mat)):
        matches = sum(1 for i in range(len(mat)) 
                      if 0 <= i - shift < len(rev_star_comp) 
                      and mat[i] == rev_star_comp[i - shift])
        if matches > max_matches:
            max_matches = matches
    return max_matches

def calc_3p_overhang(exp_line):
    """
    Calculates the 3' overhang length (bp) between 5' and 3' duplex strands
    from the miRDeep2 exp string (M=Mature, S=Star).
    """
    if not exp_line or 'M' not in exp_line or 'S' not in exp_line:
        return "NA"
        
    mat_start, mat_end = exp_line.find('M'), exp_line.rfind('M')
    star_start, star_end = exp_line.find('S'), exp_line.rfind('S')
    
    # Determine which strand is on 5' arm vs 3' arm
    if mat_start < star_start:
        arm5_start, arm5_end = mat_start, mat_end
        arm3_start, arm3_end = star_start, star_end
    else:
        arm5_start, arm5_end = star_start, star_end
        arm3_start, arm3_end = mat_start, mat_end
        
    # Calculate 3' overhang offset at duplex terminal ends
    # Distance from 3' end of 3' arm to 5' end of 5' arm
    total_len = len(exp_line)
    overhang_3p = (total_len - arm3_end) - arm5_start
    
    # If flank region is not padded, calculate via strand boundary alignment difference
    if overhang_3p < 0 or overhang_3p > 10:
        overhang_3p = abs((arm5_end - arm5_start) - (arm3_end - arm3_start))
        
    return max(0, overhang_3p)

# -------------------------------------------------------------------------
# 2. Parse miRDeep2 Summary CSV File (Novel + Known Sections)
# -------------------------------------------------------------------------

candidates = {}

with open(csv_path, 'r', encoding='utf-8', errors='ignore') as f:
    current_section = None
    
    for line in f:
        line_str = line.strip()
        if not line_str:
            continue
            
        if 'novel miRNAs predicted by miRDeep2' in line_str:
            current_section = 'novel'
            continue
        elif 'mature miRBase miRNAs detected by miRDeep2' in line_str:
            current_section = 'known'
            continue
            
        parts = line_str.split('\t')
        
        if parts[0] in ['provisional id', 'tag id']:
            continue
            
        if current_section in ['novel', 'known'] and len(parts) >= 16:
            mirna_id = parts[0]
            try:
                score = float(parts[1])
            except ValueError:
                continue
                
            rfam = parts[3]
            tot_reads = parts[4]
            mat_reads = parts[5]
            loop_reads = parts[6]
            star_reads = parts[7]
            randfold_p = parts[8]
            known_name = parts[9]
            mat_seq = parts[13]
            star_seq = parts[14]
            prec_seq = parts[15]
            
            candidates[mirna_id] = {
                'mirna_id': mirna_id,
                'category': current_section,
                'mirdeep2_score': score,
                'rfam_alert': rfam,
                'total_read_count': tot_reads,
                'mature_read_count': mat_reads,
                'loop_read_count': loop_reads,
                'star_read_count': star_reads,
                'randfold_pvalue': randfold_p,
                'known_mirbase_name': known_name if known_name != '-' else 'NA',
                'consensus_mature_seq': mat_seq,
                'mature_len_bp': len(mat_seq) if mat_seq and mat_seq != '-' else 0,
                'consensus_star_seq': star_seq,
                'star_len_bp': len(star_seq) if star_seq and star_seq != '-' else 0,
                'consensus_precursor_seq': prec_seq,
                'precursor_len_bp': len(prec_seq) if prec_seq and prec_seq != '-' else 0,
                'comp_bp': calc_complementarity(mat_seq, star_seq)
            }

print(f"Loaded {len(candidates)} total entries from {csv_path} (Novel & Known).")

# -------------------------------------------------------------------------
# 3. Parse MRD File for 5' Consistency and Numeric 3' Overhang (bp)
# -------------------------------------------------------------------------

mrd_stats = {}

with open(mrd_path, 'r', encoding='utf-8', errors='ignore') as f:
    current_id = None
    in_block = False
    reads = []
    exp_line = ""
    
    for line in f:
        line_str = line.rstrip('\r\n')
        
        if line_str.startswith('>'):
            if current_id in candidates:
                pos_counts = Counter()
                total_mapped = 0
                for start_pos, count in reads:
                    pos_counts[start_pos] += count
                    total_mapped += count
                
                pct_5p = (pos_counts.most_common(1)[0][1] / total_mapped * 100) if total_mapped > 0 else 0.0
                overhang_3p_bp = calc_3p_overhang(exp_line)
                
                mrd_stats[current_id] = {
                    'pct_5p_consistency': pct_5p,
                    'mapped_read_count': total_mapped,
                    'overhang_3p_bp': overhang_3p_bp
                }
            
            current_id = line_str[1:].split()[0]
            reads = []
            exp_line = ""
            in_block = current_id in candidates
            continue
            
        if in_block:
            if line_str.startswith('exp'):
                parts = line_str.split()
                exp_line = parts[-1] if len(parts) > 1 else ""
                
            parts = line_str.split()
            if len(parts) >= 2 and parts[0].startswith('s'):
                read_id = parts[0]
                count_match = re.search(r'_x(\d+)$', read_id)
                count = int(count_match.group(1)) if count_match else 1
                
                seq_match = re.search(r'[a-zA-Z\-]', line_str[len(read_id):])
                if seq_match:
                    start_pos = len(read_id) + seq_match.start()
                    reads.append((start_pos, count))

# Process final candidate block
if current_id in candidates:
    pos_counts = Counter()
    total_mapped = 0
    for start_pos, count in reads:
        pos_counts[start_pos] += count
        total_mapped += count
    pct_5p = (pos_counts.most_common(1)[0][1] / total_mapped * 100) if total_mapped > 0 else 0.0
    overhang_3p_bp = calc_3p_overhang(exp_line)
    
    mrd_stats[current_id] = {
        'pct_5p_consistency': pct_5p,
        'mapped_read_count': total_mapped,
        'overhang_3p_bp': overhang_3p_bp
    }

# -------------------------------------------------------------------------
# 4. Output Master Annotated Table
# -------------------------------------------------------------------------

headers = [
    "mirna_id",
    "category",
    "known_mirbase_name",
    "mirdeep2_score",
    "mature_read_count",
    "star_read_count",
    "total_read_count",
    "significant_randfold_pvalue",
    "rfam_alert",
    "consensus_mature_seq",
    "mature_len_bp",
    "consensus_star_seq",
    "star_len_bp",
    "consensus_precursor_seq",
    "precursor_len_bp",
    "comp_bp",
    "pct_5p_consistency",
    "overhang_3p_bp",
    "passes_praher_strict"
]

with open(out_path, 'w') as out:
    out.write("\t".join(headers) + "\n")
    
    for cid, info in candidates.items():
        m_stats = mrd_stats.get(cid, {'pct_5p_consistency': 0.0, 'overhang_3p_bp': 'NA'})
        pct_5p = m_stats['pct_5p_consistency']
        overhang_bp = m_stats['overhang_3p_bp']
        
        passes_strict = "YES" if (info['mirdeep2_score'] >= 10.0 and info['comp_bp'] >= 16 and pct_5p >= 90.0) else "NO"
        
        row = [
            info['mirna_id'],
            info['category'],
            info['known_mirbase_name'],
            str(info['mirdeep2_score']),
            str(info['mature_read_count']),
            str(info['star_read_count']),
            str(info['total_read_count']),
            info['randfold_pvalue'],
            info['rfam_alert'],
            info['consensus_mature_seq'],
            str(info['mature_len_bp']),
            info['consensus_star_seq'],
            str(info['star_len_bp']),
            info['consensus_precursor_seq'],
            str(info['precursor_len_bp']),
            str(info['comp_bp']),
            f"{pct_5p:.2f}",
            str(overhang_bp),
            passes_strict
        ]
        out.write("\t".join(row) + "\n")

print(f"Master dataset containing all {len(candidates)} miRNAs written to: {out_path}")
```

```
python3 parse_mcap_miRNAs.py ../../result_03_11_2026_t_01_56_18.csv output.mrd mcap_miRNA_master_dataset.tsv
```

I assessed putative miRNAs in R in the miRNA discovery [script](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/scripts/miRNA_discovery.Rmd). I filtered based on: 

- mirdeep2 score of 10 or above 
- No rfam alert 
- Significant randfold pvalue
- Consistent 5' end position
- At least 2-bp 3' overhang--had to manually check the PDFs to confirm
- Minumum of 10 mature reads 
- 16 bp complementarity between mature and star strand 


No 2023, 2036, 2037, 2050 identified. 100, 2022, 2030, and 9425 identified. 9425 did not have consistent 5' end so was not included as putative miRNA. Porites_evermani_scaffold_3072 was identified as 2025, peve-mir-novel-21 were also known ones. 

Pull putative miRNAs from novel and known fastas 

```
cd /scratch3/workspace/jillashey_uri_edu-mcap_mirna/mirna_results_03_11_2026_t_01_56_18

# Make list of miRNA IDs
awk -F',' 'NR>1 {gsub(/"/, "", $2); print $2}' final_precursor_miRNAs.csv | sort -u > target_ids.txt

# Pull miRNAs from fastas
for fa in *.fa; do
    awk '
      NR==FNR { ids[$1]; next }
      /^>/ {
        id = substr($1, 2)
        keep = (id in ids)
      }
      keep { print }
    ' target_ids.txt "$fa" > "filtered_${fa}"
done
```

Cat the novel and known sequences together 

```
cat filtered_known_mature_03_11_2026_t_01_56_18_score-50_to_na.fa filtered_novel_mature_03_11_2026_t_01_56_18_score-50_to_na.fa > mcap_miRNA_mature.fa

cat filtered_known_pres_03_11_2026_t_01_56_18_score-50_to_na.fa filtered_novel_pres_03_11_2026_t_01_56_18_score-50_to_na.fa > mcap_miRNA_precursor.fa

cat filtered_known_star_03_11_2026_t_01_56_18_score-50_to_na.fa filtered_novel_star_03_11_2026_t_01_56_18_score-50_to_na.fa > mcap_miRNA_star.fa
```

Now time to run the quantifier module of mirdeep2 which will determine the expression of the miRNAs. `nano quantifier_mirdeep2.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 5-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch3/workspace/jillashey_uri_edu-mcap_mirna

# Load modules 
module load conda/latest
conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 

# Define paths 
MAP_FA="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/mirdeep2/mapped_reads.fa"
PRECURS_FA="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/mirna_results_03_11_2026_t_01_56_18/mcap_miRNA_precursor.fa"
MATURE_FA="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/mirna_results_03_11_2026_t_01_56_18/mcap_miRNA_mature.fa"
STAR_FA="/scratch3/workspace/jillashey_uri_edu-mcap_mirna/mirna_results_03_11_2026_t_01_56_18/mcap_miRNA_star.fa"

echo "Quantifying miRNAs with mirdeep2" $(date)

quantifier.pl -p "${PRECURS_FA}" -m "${MATURE_FA}" -r "${MAP_FA}" -s "${STAR_FA}" -c config.txt 

echo "Quantifying complete!" $(date)

conda deactivate
```

Submitted batch job 1904546.

Low mapping rate. I think mirdeep2 is confused because the fastas all have the same headers. Add mature or star to header names 

```
# Keep precursor file as is (or ensure it uses T bases)
awk '{if ($0 ~ /^>/) print $1; else {gsub(/U/,"T"); gsub(/u/,"t"); print $0}}' \
  mcap_miRNA_precursor.fa > precursor_clean.fa

# Append '_mature' to mature FASTA headers
awk '{if ($0 ~ /^>/) print $1"_mature"; else {gsub(/U/,"T"); gsub(/u/,"t"); print $0}}' \
  mcap_miRNA_mature.fa > mature_clean.fa

# Append '_star' to star FASTA headers
awk '{if ($0 ~ /^>/) print $1"_star"; else {gsub(/U/,"T"); gsub(/u/,"t"); print $0}}' \
  mcap_miRNA_star.fa > star_clean.fa
  
head -n 4 precursor_clean.fa mature_clean.fa star_clean.fa
==> precursor_clean.fa <==
>Montipora_capitata_HIv3___Scaffold_2_271104
CAATGTTTCGGCTTGTTCCCGTTTTCGGGAACAAGTCGAAACATTAGA
>Montipora_capitata_HIv3___Scaffold_8_1538851
CAATGTTTCGGCTTGTTCCCGTTTTCGGAAGCAAGCCGAAACATTGTA
==> mature_clean.fa <==
>Montipora_capitata_HIv3___Scaffold_2_271104_mature
CAATGTTTCGGCTTGTTCCCG
>Montipora_capitata_HIv3___Scaffold_8_1538851_mature
CAATGTTTCGGCTTGTTCCCG
==> star_clean.fa <==
>Montipora_capitata_HIv3___Scaffold_2_271104_star
GGAACAAGTCGAAACATTAGA
>Montipora_capitata_HIv3___Scaffold_8_1538851_star
GAAGCAAGCCGAAACATTGTA
```

Rerun quantifier with the clean fastas as input. Submitted batch job 1906954. Not much difference. 

try gene ext for polyA libraries?



/home/jillashey_uri_edu/.conda/envs/targetfinder

For target prediction

- also try https://github.com/carringtonlab/TargetFinder
- remove strict flag from miranda?
- https://github.com/jtremblay/mirnatarget