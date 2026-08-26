---
layout: post
title: Annotations for cnidarian sperm analysis
date: '2026-08-25'
categories: sperm
tags: Cnidarian sperm
---


## Perform annotation on proteomes 

Download swissprot annotations 

```
cd /scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations

# Swissprot protein seqs
curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
mv uniprot_sprot.fasta.gz uniprot_sprot_r20260825.fasta.gz
gunzip -k uniprot_sprot_r20260825.fasta.gz

# Swissprot GO 
curl -H "Accept: text/plain; format=tsv" "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cgo_id&format=tsv&query=%28*%29+AND+%28reviewed%3Atrue%29" -o SwissProt-Annot-GO_20260825.tsv
```

Make the blast db for the Swissprot fasta. `nano make_db.sh`

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
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2

# Load modules 
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Making db" $(date)

makeblastdb -in annotations/dbs/uniprot_sprot_r20260825.fasta -dbtype prot -out annotations/dbs/sprot_db

echo "DB construction complete!" $(date)
```

Submitted batch job 63610450


Annotate the Apoc proteins. `nano apoc_blastp.sh`

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
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2

# Load modules 
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Annotating Apoc proteins with Swissprot db" $(date)

fasta="protein_fastas/Apoc_proteins.faa"

blastp -query $fasta -db annotations/dbs/sprot_db -out annotations/apoc_blastp_out.tab -evalue 1E-05 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "Apoc annotation complete" $(date)
```

Submitted batch job 63610633

Merge the blast results with the GO data 

```
cd /scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {go[$1]=$0; next} {split($2, a, "|"); acc=a[2]; if(acc in go) print $0, go[acc]}' dbs/SwissProt-Annot-GO_20260825.tsv apoc_blastp_out.tab > apoc_blast_GO.tsv
```

Annotate the Ahya proteins. `nano ahya_blastp.sh`

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
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2

# Load modules 
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Annotating Ahya proteins with Swissprot db" $(date)

fasta="protein_fastas/Ahya_proteins.faa"

blastp -query $fasta -db annotations/dbs/sprot_db -out annotations/ahya_blastp_out.tab -evalue 1E-05 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "Ahya annotation complete" $(date)
```

Submitted batch job 63610658.

Merge the blast results with the GO data 

```
cd /scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {go[$1]=$0; next} {split($2, a, "|"); acc=a[2]; if(acc in go) print $0, go[acc]}' dbs/SwissProt-Annot-GO_20260825.tsv ahya_blastp_out.tab > ahya_blast_GO.tsv
```

Annotate the Nvec proteins. `nano nvec_blastp.sh`

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
#SBATCH -D /scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2

# Load modules 
module load uri/main all/BLAST+/2.15.0-gompi-2023a

echo "Annotating Nvec proteins with Swissprot db" $(date)

fasta="protein_fastas/Nvec_proteins.faa"

blastp -query $fasta -db annotations/dbs/sprot_db -out annotations/nvec_blastp_out.tab -evalue 1E-05 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "Nvec annotation complete" $(date)
```

Submitted batch job 63610666

Merge the blast results with the GO data 

```
cd /scratch4/workspace/jillashey_uri_edu-cnidarian_sperm_part2/annotations
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {go[$1]=$0; next} {split($2, a, "|"); acc=a[2]; if(acc in go) print $0, go[acc]}' dbs/SwissProt-Annot-GO_20260825.tsv nvec_blastp_out.tab > nvec_blast_GO.tsv
```


