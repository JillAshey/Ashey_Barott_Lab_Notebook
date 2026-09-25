---
layout: post
title: POC comparison
date: '2026-09-24'
categories: POC
tags: POC
---

## POC species comparisons

Pocillopora are so weird. I'm curious to do some comparative gene analysis to see how similar their genomes/proteins are. I'm going to run broccoli on all available POC proteins and then annotate them with swissprot so they are all relatively similar. 

I did not include the Pacuta genome from [Vidal-Dupiol et al. 2020](https://www.biorxiv.org/content/10.1101/698688v3.full), as the link to access the genome and genomic resources does not seem to work anymore. 

### Broccoli

Running [broccoli](https://github.com/rderelle/Broccoli/tree/master) to identify orthologous groups across POC species. 

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
	- Accession: PRJNA454489 or http://pdam.reefgenomics.org/ (I used the one from reef genomics)
	- Determined to be P. grandis (Oury et al. 2023) or P. capitata (Connelly et al. 2025)
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


Download POC protein data 

```
cd /scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins

# Pmea
wget http://cyanophora.rutgers.edu/Pocillopora_meandrina/Pocillopora_meandrina_HIv1.genes.pep.faa.gz

# Ptua 
### location on unity 
ln -s /scratch4/workspace/jillashey_uri_edu-Ptua_genome/braker3_output/braker_clean.aa .

# Pacuta HI
wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.genes.pep.faa.gz

# Pdam 
wget http://pdam.reefgenomics.org/download/pdam_proteins.fasta.gz

# Pver 
wget http://pver.reefgenomics.org/download/Pver_proteins_names_v1.0.faa.gz

# Peff
wget https://www.genoscope.cns.fr/corals/data/Pocillopora_effusa_v3.annot.pep.fa

# Unzip files 
gunzip *
```

Clean up protein header names so we know who is who

```
sed 's/^\(>pdam_[0-9]*\)-RA.*/\1-RA/' pdam_proteins.fasta > pdam_proteins_clean.fasta
rm pdam_proteins.fasta

sed '/^>/ s/$/_Ptua/' braker_clean.aa > braker_clean_Ptua.aa
unlink braker_clean.aa
```

Rename all fasta files so its easier to understand output

```
mv braker_clean_Ptua.aa Ptua_protein.faa
mv pdam_proteins_clean.fasta Pdam_proteins.faa
mv Pocillopora_acuta_HIv2.genes.pep.faa Pacu_proteins.faa
mv Pocillopora_effusa_v3.annot.pep.fa Peff_proteins.faa
mv Pocillopora_meandrina_HIv1.genes.pep.faa Pmea_proteins.faa
mv Pver_proteins_names_v1.0.faa Pver_proteins.faa
```

Run Broccoli. `nano POC_broccoli.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load conda and activate conda environment
module load conda/latest
conda activate /work/pi_hputnam_uri_edu/conda/envs/env-broccoli

# Load additional programs needed to run broccoli
module load uri/main
module load diamond/2.1.7
module load all/FastTree/2.1.11-GCCcore-12.3.0

# Define directory paths 
PROT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins/broccoli"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Running broccoli on all POC proteins" $(date)

cd "${OUTPUT_DIR}"

python /work/pi_hputnam_uri_edu/conda/envs/env-broccoli/Broccoli/broccoli.py -dir "${PROT_DIR}" -phylogenies 'ml' -ext '.faa' -path_fasttree FastTree

echo "Broccoli complete!" $(date)

conda deactivate
```

Submitted batch job 64836362

### Proteome annotations 

Using the swissprot db that I built when annotating proteomes for cnidarian sperm project (see [post](https://github.com/JillAshey/Ashey_Barott_Lab_Notebook/blob/main/posts/2026-08-25-Functional-Annotation-Apoc-Ahya-Nvec.md)). 

Blast POC proteins to swissprot db. `nano pmea_blast.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Annotating Pmea proteins with Swissprot db $(date)"

fasta="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins/Pmea_proteins.faa"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins"
GO_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations"
DB_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations/dbs"

# Run BLAST search
blastp -query "${fasta}" \
       -db "${DB_DIR}/sprot_db" \
       -out "${OUTPUT_DIR}/pmea_blastp_out.tab" \
       -evalue 1E-05 \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -outfmt 6

# Merge with GO data
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {go[$1]=$0; next} {split($2, a, "|"); acc=a[2]; if(acc in go) print $0, go[acc]}' \
    "${GO_DIR}/SwissProt-Annot-GO_20260825.tsv" \
    "${OUTPUT_DIR}/pmea_blastp_out.tab" > "${OUTPUT_DIR}/pmea_blast_GO.tsv"

echo "Pmea annotation complete $(date)"
```

Submitted batch job 64838860

`nano ptua_blast.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Annotating Ptua proteins with Swissprot db $(date)"

fasta="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins/Ptua_proteins.faa"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins"
GO_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations"
DB_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations/dbs"

# Run BLAST search
blastp -query "${fasta}" \
       -db "${DB_DIR}/sprot_db" \
       -out "${OUTPUT_DIR}/ptua_blastp_out.tab" \
       -evalue 1E-05 \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -outfmt 6

# Merge with GO data
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {go[$1]=$0; next} {split($2, a, "|"); acc=a[2]; if(acc in go) print $0, go[acc]}' \
    "${GO_DIR}/SwissProt-Annot-GO_20260825.tsv" \
    "${OUTPUT_DIR}/ptua_blastp_out.tab" > "${OUTPUT_DIR}/ptua_blast_GO.tsv"

echo "Ptua annotation complete $(date)"
```

Submitted batch job 64840055

`nano pacu_blast.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=25GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Annotating Pacu proteins with Swissprot db $(date)"

fasta="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins/Pacu_proteins.faa"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins"
GO_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations"
DB_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations/dbs"

# Run BLAST search
blastp -query "${fasta}" \
       -db "${DB_DIR}/sprot_db" \
       -out "${OUTPUT_DIR}/pacu_blastp_out.tab" \
       -evalue 1E-05 \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -outfmt 6

# Merge with GO data
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {go[$1]=$0; next} {split($2, a, "|"); acc=a[2]; if(acc in go) print $0, go[acc]}' \
    "${GO_DIR}/SwissProt-Annot-GO_20260825.tsv" \
    "${OUTPUT_DIR}/pacu_blastp_out.tab" > "${OUTPUT_DIR}/pacu_blast_GO.tsv"

echo "Pacu annotation complete $(date)"
```

Submitted batch job 64839048

`nano pdam_blast.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=25GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Annotating Pdam proteins with Swissprot db $(date)"

fasta="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins/Pdam_proteins.faa"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins"
GO_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations"
DB_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations/dbs"

# Run BLAST search
blastp -query "${fasta}" \
       -db "${DB_DIR}/sprot_db" \
       -out "${OUTPUT_DIR}/pdam_blastp_out.tab" \
       -evalue 1E-05 \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -outfmt 6

# Merge with GO data
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {go[$1]=$0; next} {split($2, a, "|"); acc=a[2]; if(acc in go) print $0, go[acc]}' \
    "${GO_DIR}/SwissProt-Annot-GO_20260825.tsv" \
    "${OUTPUT_DIR}/pdam_blastp_out.tab" > "${OUTPUT_DIR}/pdam_blast_GO.tsv"

echo "Pdam annotation complete $(date)"
```

Submitted batch job 64839957

`nano pver_blast.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=25GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Annotating Pver proteins with Swissprot db $(date)"

fasta="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins/Pver_proteins.faa"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins"
GO_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations"
DB_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations/dbs"

# Run BLAST search
blastp -query "${fasta}" \
       -db "${DB_DIR}/sprot_db" \
       -out "${OUTPUT_DIR}/pver_blastp_out.tab" \
       -evalue 1E-05 \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -outfmt 6

# Merge with GO data
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {go[$1]=$0; next} {split($2, a, "|"); acc=a[2]; if(acc in go) print $0, go[acc]}' \
    "${GO_DIR}/SwissProt-Annot-GO_20260825.tsv" \
    "${OUTPUT_DIR}/pver_blastp_out.tab" > "${OUTPUT_DIR}/pver_blast_GO.tsv"

echo "Pver annotation complete $(date)"
```

Submitted batch job 64839966

`nano peff_blast.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=25GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-POC_Heron/

# Load modules
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Annotating Peff proteins with Swissprot db $(date)"

fasta="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins/Peff_proteins.faa"
OUTPUT_DIR="/scratch4/workspace/jillashey_uri_edu-POC_Heron/refs/proteins"
GO_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations"
DB_DIR="/scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations/dbs"

# Run BLAST search
blastp -query "${fasta}" \
       -db "${DB_DIR}/sprot_db" \
       -out "${OUTPUT_DIR}/peff_blastp_out.tab" \
       -evalue 1E-05 \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -outfmt 6

# Merge with GO data
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {go[$1]=$0; next} {split($2, a, "|"); acc=a[2]; if(acc in go) print $0, go[acc]}' \
    "${GO_DIR}/SwissProt-Annot-GO_20260825.tsv" \
    "${OUTPUT_DIR}/peff_blastp_out.tab" > "${OUTPUT_DIR}/peff_blast_GO.tsv"

echo "Peff annotation complete $(date)"
```

Submitted batch job 64839990

Check how many protein sequences were annotated relative to total

```
grep -c ">" *.faa
Pacu_proteins.faa:33730
Pdam_proteins.faa:26077
Peff_proteins.faa:32095
Pmea_proteins.faa:31840
Ptua_proteins.faa:32520
Pver_proteins.faa:27439


```




failing at the merge w/ GO info step -- may have to do that after it finishes running 