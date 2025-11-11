---
layout: post
title: Sperm smRNA analysis 
date: '2025-08-26'
categories: Analysis
tags: [Bioinformatics, small RNA]
projects: Sperm smRNA 
---

## Sperm smRNA analysis 

Sperm from three cnidarian species (*Nematostella vectensis*, *Astrangia poculata*, and *Acropora hyacinthus*) was collected, and small RNAs were extracted and sequenced from the sperm. This post includes the details on the sperm smRNA analysis. 

The fastq files were initially stored on Box but I have moved them to Unity for analysis (`/project/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA`). My scripts and output will be stored here: `/work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA`. 

Run initial QC. In the scripts folder: `nano raw_qc.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /project/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA

echo "Running md5sum" $(date)
md5sum *fastq.gz > checkmd5.txt
md5sum -c checkmd5.txt 
mv checkmd5.txt /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data

echo "Counting number of reads in each file" $(date)
zgrep -c "@VL" *fastq.gz

# Create an array of fastq files to process
files=($('ls' *.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/fastqc_raw && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/fastqc_raw

echo "Starting multiqc..." $(date)
multiqc *

echo "Initial QC of smRNA data complete." $(date)
```

Submitted batch job 41679618. Ran in about 8 mins. Here are the number of raw reads per sample: 

```
ahya_1_S27_L001_R1_001.fastq.gz:4893728
ahya_2_S28_L001_R1_001.fastq.gz:9123772
ahya_3_S29_L001_R1_001.fastq.gz:6858991
ahya_4_S30_L001_R1_001.fastq.gz:7389533
apoc_2_S31_L001_R1_001.fastq.gz:6896525
apoc_3_S32_L001_R1_001.fastq.gz:6790067
apoc_4_S33_L001_R1_001.fastq.gz:6194770
nvec_1_S25_L001_R1_001.fastq.gz:6380826
nvec_2_S26_L001_R1_001.fastq.gz:8530065
nvec_3_S34_L001_R1_001.fastq.gz:7685020
nvec_4_S35_L001_R1_001.fastq.gz:7667363
```

MultiQC results are [here](https://github.com/JillAshey/Cnidarian_sperm_smRNA/blob/main/data/multiqc_report_raw.html). 

The QC looks good overall. Reads are 75bp long. There is a lot of adapter content from the Illumina 3' smRNA adapter. Before trimming, make a scratch workspace for intermediate files: 

```
ws_allocate cnidarian_sperm 30
Info: could not read email from users config ~/.ws_user.conf.
Info: reminder email will be sent to local user account
Info: creating workspace.
/scratch3/workspace/jillashey_uri_edu-cnidarian_sperm
remaining extensions  : 5
remaining time in days: 30
```

Run fastp to trim files. In the scripts folder: `nano trim_and_qc.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

# load modules needed
module load uri/main
module load parallel/20240822
module load fastqc/0.12.1
module load all/MultiQC/1.12-foss-2021b
module load fastp/0.23.2-GCC-11.2.0

cd /project/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA

echo "Starting trimming with fastp" $(date)

for fq in *.fastq.gz
do
  base=$(basename $fq .fastq.gz)
  fastp \
    --in1 $fq \
    --out1 /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/${base}_trim.fastq.gz \
    --adapter_sequence TGGAATTCTCGGGTGCCAAGG \
    --qualified_quality_phred 30 \
    --length_required 15 \
    --length_limit 50 \
    #--dedup \
    --html ${base}_fastp_report.html \
    --thread 4
done

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm

echo "Counting number of reads in trimmed file" $(date)
zgrep -c "@VL" *fastq.gz

# Create an array of fastq files to process
files=($('ls' *.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc for trimmed files..." $(date)
parallel -j 20 "fastqc {} -o /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/fastqc_trim && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/fastqc_trim

echo "Starting multiqc..." $(date)
multiqc *

echo "QC of trimmed smRNA data complete." $(date)
```

Submitted batch job 41705347. Here are the number of trimmed reads per sample: 

```
ahya_1_S27_L001_R1_001_trim.fastq.gz:4015652
ahya_2_S28_L001_R1_001_trim.fastq.gz:2643869
ahya_3_S29_L001_R1_001_trim.fastq.gz:6131770
ahya_4_S30_L001_R1_001_trim.fastq.gz:6454522
apoc_2_S31_L001_R1_001_trim.fastq.gz:6663244
apoc_3_S32_L001_R1_001_trim.fastq.gz:6141174
apoc_4_S33_L001_R1_001_trim.fastq.gz:5233509
nvec_1_S25_L001_R1_001_trim.fastq.gz:6182908
nvec_2_S26_L001_R1_001_trim.fastq.gz:8258638
nvec_3_S34_L001_R1_001_trim.fastq.gz:7250936
nvec_4_S35_L001_R1_001_trim.fastq.gz:6457956
```

MultiQC results are [here](https://github.com/JillAshey/Cnidarian_sperm_smRNA/blob/main/data/multiqc_report_trim.html). 

Data looks good. Ahya 2 dropped quite a bit in number of reads but that sample also had the most adapter content. 

### miRNAs

Install shortstack on Unity. See [shortstack github](https://github.com/MikeAxtell/ShortStack?tab=readme-ov-file#installation) for installation instructions. 

```
cd /work/pi_hputnam_uri_edu/conda/envs
module load conda/latest # need to load before making any conda envs
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create --prefix /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 shortstack
conda activate /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 
ShortStack --version
ShortStack 4.1.2                                                                         
```

Success! Downloaded these files as references: 

- Cnidarian mirbase [fasta](https://github.com/urol-e5/deep-dive/blob/main/data/cnidarian-mirbase-mature-v22.1-no_U.fa)
- [Apul miRNA fasta](https://github.com/urol-e5/deep-dive-expression/blob/main/D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/Apul_ShortStack_4.1.0_mature.fasta) from deep dive expression
- [Peve miRNA fasta](https://github.com/urol-e5/deep-dive-expression/blob/main/E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/Peve_ShortStack_4.1.0_mature.fasta)
- [Ptuh/Pmea miRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta)
- [Apoc miRNA fasta](https://github.com/JillAshey/Astrangia_repo/blob/main/output/Molecular/smRNA/shortstack/mir_coords_fixed_mature.fasta)

Add spp names to the header IDs of the Apul, Peve, Ptuh, and Apoc fasta files.
```
sed -i 's/^>\([^ ]*\)/>\1::Acropora_pulchra/' Apul_ShortStack_4.1.0_mature.fasta 
sed -i 's/^>\([^ ]*\)/>\1::Porites_evermanni/' Peve_ShortStack_4.1.0_mature.fasta
sed -i 's/^>\([^ ]*\)/>\1::Astrangia_poculata/' mir_coords_fixed_mature.fasta
sed -i 's/^>\([^ ]*\)/>\1::Pocillopora_meandrina/' mir.fasta 
```

Cat files together for one miRNA reference fasta file 

```
cat Apul_ShortStack_4.1.0_mature.fasta cnidarian-mirbase-mature-v22.1-no_U.fa mir_coords_fixed_mature.fasta mir.fasta Peve_ShortStack_4.1.0_mature.fasta > miRNA_reference.fasta
grep -c ">" miRNA_reference.fasta 
49661
```

Organize trimmed data in species-specific folders 

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm
mkdir apoc ahya nvec
mv ahya_* ahya
mv apoc_* apoc
mv nvec_* nvec
```

Okay now time to actually run shortstack on the trimmed files. Only going to do apoc and nvec for now--not sure what genome to use for Ahya. In the scripts folder: `nano shortstack_apoc_nvec.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 

echo "Running short stack on trimmed reads for apoc"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 

gunzip *fastq.gz

ShortStack \
--genomefile /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.assembly.scaffolds_chromosome_level.fasta \
--readfile apoc_2_S31_L001_R1_001_trim.fastq \
apoc_3_S32_L001_R1_001_trim.fastq \
apoc_4_S33_L001_R1_001_trim.fastq \
--known_miRNAs /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta \
--outdir /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/apoc_shortstack \
--dn_mirna

echo "Apoc shortstack complete, running short stack on trimmed reads for nvec"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec 

gunzip *fastq.gz

ShortStack \
--genomefile /work/pi_hputnam_uri_edu/genomes/Nvec/Nvec_genome.fa \
--readfile nvec_1_S25_L001_R1_001_trim.fastq \
nvec_2_S26_L001_R1_001_trim.fastq \
nvec_3_S34_L001_R1_001_trim.fastq \
nvec_4_S35_L001_R1_001_trim.fastq \
--known_miRNAs /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta \
--outdir /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/nvec_shortstack \
--dn_mirna

echo "Nvec shortstack complete"
echo "Shortstack complete for Apoc and Nvec"

conda deactivate
```

Submitted batch job 42213339

Didn't identify many miRNAs...why???? you know why...because maybe the reads are trimmed too short...because the total hairpin miRNA can be quite long (80-90bp)...but didn't we trim the seqs to 35 for e5? And we used shortstack there...maybe its the genomes? I need to find the genomes online and redownload 

- https://zenodo.org/records/14110456 --apoc
- https://www.nature.com/articles/s41467-023-44080-7 / https://simrbase.stowers.org/starletseaanemone -- nvec
- https://marinegenomics.oist.jp/ahya/viewer/download?project_id=91 OR https://drive.google.com/drive/u/1/folders/15dh6eE8q850frTRVtglndXtzxSwj673g - ahya

Redownloaded Apoc genome. I may also need to unzip the genome...

```
cd /work/pi_hputnam_uri_edu/genomes/Apoc
gunzip apoculata.genome.fasta.gz 
```

Run only apoc shortstack. `nano shortstack_apoc.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 

echo "Running short stack on trimmed reads for apoc"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 

#gunzip *fastq.gz

ShortStack \
--genomefile /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta \
--readfile apoc_2_S31_L001_R1_001_trim.fastq \
apoc_3_S32_L001_R1_001_trim.fastq \
apoc_4_S33_L001_R1_001_trim.fastq \
--known_miRNAs /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta \
--outdir /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/apoc_shortstack \
--dn_mirna

echo "Apoc shortstack complete"
conda deactivate
```

Submitted batch job 42222961. Still got the same result...Results [here](https://github.com/JillAshey/Cnidarian_sperm_smRNA/tree/main/output/shortstack/apoc/trimmed_default). I wonder if I can relax the shortstack parameters? I'm going to look more at the shortstack github [wiki](https://github.com/MikeAxtell/ShortStack/wiki/Vignette-%231-%3A-A-%22complete%22-run). Interestingly, the recommended inputs are "one or more untrimmed FASTQ files of small RNA-seq data". Shortstack has an `--autotrim` flag, which (according to the wiki): "tells ShortStack that the reads from --readfile need to be adapter trimmed, and that ShortStack should infer the adapter sequence. This is the recommended way to trim the adapters off of sRNA-seq data." Let's try running the same code but with the raw reads instead. 

In the scripts folder: `nano shortstack_apoc_rawreads.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 

echo "Running short stack on raw reads for apoc"

cd /project/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA

#gunzip apoc*

ShortStack \
--genomefile /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta \
--readfile apoc_2_S31_L001_R1_001.fastq \
apoc_3_S32_L001_R1_001.fastq \
apoc_4_S33_L001_R1_001.fastq \
--autotrim \
--known_miRNAs /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta \
--outdir /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/apoc_raw_shortstack \
--dn_mirna

echo "Apoc shortstack complete"
conda deactivate
```

Submitted batch job 42223404. Still the same problem...Results [here](https://github.com/JillAshey/Cnidarian_sperm_smRNA/tree/main/output/shortstack/apoc/raw_default). Only IDing three miRNAs. It seems to be IDing other known miRNAs but just not marking them as miRNAs? Looking at the `Results.txt` file: 

```
Locus	Name	Chrom	Start	End	Length	Reads	DistinctSequences	FracTop	Strand	MajorRNA	MajorRNAReads	Short	Long	21	22	23	24	DicerCall	MIRNA	known_miRNAs
Ap11:28741-29152	Cluster_1	Ap11	28741	29152	412	31	7	1.0	+	CCCCCCUCCCGGC	10	31	0	0	0	0	0	N	N	NA
Ap11:70294-70715	Cluster_2	Ap11	70294	70715	422	92	16	0.989	+	UAAGACUCUGAGGAUGGAAUGG	33	40	0	1	34	13	4	N	N	NA
Ap11:171278-171695	Cluster_3	Ap11	171278	171695	418	29	5	0.0	-	AUUUGGUUGGUUAUUGGA	19	29	0	0	0	0	0	N	N	NA
Ap11:423611-424025	Cluster_4	Ap11	423611	424025	415	45	12	1.0	+	GUAGAGCGCAACCU	21	45	0	0	0	0	0	N	N	NA
Ap11:583143-583559	Cluster_5	Ap11	583143	583559	417	52	4	0.0	-	GUUCAAAUCCGGCUCCC	39	52	0	0	0	0	0	N	N	NA
Ap11:636984-637400	Cluster_6	Ap11	636984	637400	417	42	11	0.976	+	GUUCGAAUCCUGCUGCCC	11	42	0	0	0	0	0	N	N	NA
Ap11:735982-736399	Cluster_7	Ap11	735982	736399	418	22	8	0.0	-	AAUAAAUACUGAUACAGGGCU	6	16	0	6	0	0	0	N	N	NA
Ap11:818209-818623	Cluster_8	Ap11	818209	818623	415	31	7	0.871	+	CUGGUUGCUGCAGUU	21	27	4	0	0	0	0	N	N	NA
Ap11:849701-850117	Cluster_9	Ap11	849701	850117	417	168	5	0.0	-	AUGUGAGAGUAGGUCGU	141	168	0	0	0	0	0	N	N	NA
Ap11:1173168-1173582	Cluster_10	Ap11	1173168	1173582	415	32	5	1.0	+	CGUGAGGCUUAACCU	25	32	0	0	0	0	0	N	N	NA
Ap11:1291990-1292408	Cluster_11	Ap11	1291990	1292408	419	116	4	1.0	+	GAGAUUGAUAUUAGGCU	87	116	0	0	0	0	0	N	N	NA
Ap11:1507811-1508227	Cluster_12	Ap11	1507811	1508227	417	51	3	0.0	-	CUACGUCCUGAAGUGCC	47	51	0	0	0	0	0	N	N	NA
```

Most clusters have 0 in the 21-22 size categories (columns 15 and 16), which is typical miRNA size. My trimmed data looks like this: 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/trimmed_fastqc_sequence_length_distribution_plot.png)

There is little to no peak around 21-22bp, where I would expect miRNAs to be. There are lots of peaks between 27-35, where I could expect piRNAs to be. Maybe its mostly piRNAs? Maybe I can relax some of the shortstack parameters to see if that helps? 

```
usage: ShortStack [-h] [--version] (--genomefile GENOMEFILE | --autotrim_only) [--known_miRNAs KNOWN_MIRNAS] (--readfile [READFILE ...] | --bamfile [BAMFILE ...]) [--outdir OUTDIR]
                  [--adapter ADAPTER | --autotrim] [--autotrim_key AUTOTRIM_KEY] [--threads THREADS] [--mmap {u,f,r}] [--align_only] [--dicermin DICERMIN] [--dicermax DICERMAX]
                  [--locifile LOCIFILE | --locus LOCUS] [--nohp] [--dn_mirna] [--strand_cutoff STRAND_CUTOFF] [--mincov MINCOV] [--pad PAD] [--make_bigwigs]
options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --genomefile GENOMEFILE
                        FASTA file of the reference genome
  --autotrim_only       If this switch is set, ShortStack quits after performing auto-trimming of input reads.
  --known_miRNAs KNOWN_MIRNAS
                        FASTA file of known/suspected mature microRNAs
  --readfile [READFILE ...]
                        One or more files of reads (fq, fa, gzip OK)
  --bamfile [BAMFILE ...]
                        One or more BAM alignment files
  --outdir OUTDIR       Output directory name. Defaults to ShortStack_time
  --adapter ADAPTER     3-primer adapter sequence to trim off ofreads. If given applies to all input fastq files. Mutually exclusive with --autotrim.
  --autotrim            If this switch is set, automatically discover the 3-prime adapter from each input readfile, and trim it. This uses the sequence from --autotrim_key to discover the
                        adapter sequence. Mutually exlcusive with --adapter.
  --autotrim_key AUTOTRIM_KEY
                         Sequence of an abundant, known small RNA to be used to discover the 3-prime adapter sequence. Has no effect unless --autotrim is specified. Defaults to
                        TCGGACCAGGCTTCATTCCCC (miR166). Can be upper or lower-case, T or U and must be 20-30 bases long.
  --threads THREADS     Number of threads to use (integer) - default: 1
  --mmap {u,f,r}        Protocol for multi-mapped reads: u, f, or r - default: u
  --align_only          If this switch is set, ShortStack quits after performing alignments without any analyses performed.
  --dicermin DICERMIN   Minimum size of a valid Dicer-processed small RNA. Must be integer >= 15 and <= dicermax. Default: 20.
  --dicermax DICERMAX   Maximum size of a valid Dicer-processed small RNA. Must be integer >= 15 and <= dicermax. Default: 24.
  --locifile LOCIFILE   File listing intervals to analyze. Can be simple tab-delimited, .bed, or .gff3. Tab-delimited format is column 1 with coordinates Chr:start-stop, column 2 with names.
                        Input file assumed to be simple tab-delimited unless file name ends in .bed or .gff3. Mutually exclusive with --locus.
  --locus LOCUS         Analyze the specified interval, given in format Chr:start-stop. Mutually exclusive with --locifile.
  --nohp                If this switch is set, RNA folding will not take place, thus MIRNA loci cannot be annotated. This does however save CPU time.
  --dn_mirna            If this switch is set, a de novo search for new MIRNA loci will be performed. By default, de novo MIRNA finding is not performed and MIRNA searches are limited to loci
                        matching RNAs from --known_miRNAs that align to the genome
  --strand_cutoff STRAND_CUTOFF
                        Cutoff for calling the strandedness of a small RNA locus. Must be a floating point > 0.5 and < 1. Default: 0.8.
  --mincov MINCOV       Minimum alignment depth required to nucleate a small RNA cluster during de novo cluster search. In units of reads per million. Must be a floating point number. Default:
                        0.5
  --pad PAD             Initial peaks (continuous regions with depth exceeding argument mincov are merged if they are this distance or less from each other. Must be an integer >= 1. Default: 200
  --make_bigwigs        If this switch is set then bigwigs will be made from alignments.
```

Lets try to decrease the minimum coverage to see if miRNAs are just lowly expressed. Decreasing from 0.5 to 0.1. Also going back to using the trimmed data instead of letting shortstack do the trimming. 

In the scripts folder: `nano shortstack_mincov0.1_apoc.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 

echo "Running short stack on trimmed reads for apoc with mincov set to 0.1"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 

#gunzip *fastq.gz

ShortStack \
--genomefile /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta \
--readfile apoc_2_S31_L001_R1_001_trim.fastq \
apoc_3_S32_L001_R1_001_trim.fastq \
apoc_4_S33_L001_R1_001_trim.fastq \
--mincov 0.1 \
--known_miRNAs /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta \
--outdir /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/apoc_shortstack_mincov0.1 \
--dn_mirna

echo "Apoc shortstack complete"

conda deactivate
```

Submitted batch job 42225432. Still not identifying a lot of miRNAs. See results [here](https://github.com/JillAshey/Cnidarian_sperm_smRNA/tree/main/output/shortstack/apoc/mincov0.1). In the results.txt file, I am seeing though that we are getting hits to the known miRNAs (including the ones from the Apoc adults). 

```
head(mirna_known)
                   Locus         Name Chrom    Start      End Length Reads DistinctSequences FracTop Strand                MajorRNA
1 Ap11:16932436-16932852 Cluster_1246  Ap11 16932436 16932852    417     3                 2       0      -       GUAGAACCUCAAUGAGC
2 Ap11:17370295-17370714 Cluster_1284  Ap11 17370295 17370714    420     3                 1       1      +    GUUAGAGCAACUGUUAAAGG
3   Ap10:2155340-2155761 Cluster_1831  Ap10  2155340  2155761    422     6                 6       1      +   CUGAAAUUAAAGGACUUGCCU
4 Ap10:13290211-13290627 Cluster_2700  Ap10 13290211 13290627    417     2                 1       1      +       GCUGAAUCCUCCAGGUC
5   Ap14:3503170-3503592 Cluster_3563  Ap14  3503170  3503592    423     5                 3       1      + UGAUCUUGUAGGUUCCAUCUUAU
6   Ap14:4191041-4191458 Cluster_3626  Ap14  4191041  4191458    418     3                 3       0      -      AGGAUUGUAGUGAGAACA
  MajorRNAReads Short Long X21 X22 X23 X24 DicerCall MIRNA                                                              known_miRNAs
1             2     2    0   1   0   0   0         N     N Cluster_137.mature::chromosome_1:16932636-16932657(-)::Astrangia_poculata
2             3     3    0   0   0   0   0         N     N                                                           mmu-miR-466i-5p
3             1     2    0   2   2   0   0         N     N   Cluster_247.mature::chromosome_2:2155540-2155561(+)::Astrangia_poculata
4             2     2    0   0   0   0   0         N     N                                           mmu-miR-466i-5p;mmu-miR-466i-5p
5             3     2    0   0   0   3   0         N     N   Cluster_492.mature::chromosome_3:3503370-3503391(+)::Astrangia_poculata
6             1     3    0   0   0   0   0         N     N                                           mmu-miR-466i-5p;mmu-miR-466i-5p
  mature_length
1            17
2            20
3            21
4            17
5            23
6            18
```

Very strange...I am so confused haha. Requested for mirdeep2 to be added as module to Unity. 

Try adjusting the dicer min and max. The default is 21 and 24, respectively. Going to set it to 18 and 26, respectively. 

In the scripts folder: `nano shortstack_mincov0.1_dmin18_dmax26_apoc.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 

echo "Running short stack on trimmed reads for apoc with mincov set to 0.1, dicer min = 18, and dicer max = 26"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 

#gunzip *fastq.gz

ShortStack \
--genomefile /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta \
--readfile apoc_2_S31_L001_R1_001_trim.fastq \
apoc_3_S32_L001_R1_001_trim.fastq \
apoc_4_S33_L001_R1_001_trim.fastq \
--mincov 0.1 \
--dicermin 18 \
--dicermax 26 \
--known_miRNAs /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta \
--outdir /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/apoc_shortstack_mincov0.1_dmin18_dmax26 \
--dn_mirna

echo "Apoc shortstack complete"

conda deactivate
```

Submitted batch job 42239097. Doesn't look much different bleh; see results [here](https://github.com/JillAshey/Cnidarian_sperm_smRNA/tree/main/output/shortstack/apoc/mincov0.1_dmin18_dmax26). 

Let's try to install mirdeep2. 

```
cd /work/pi_hputnam_uri_edu/conda/envs
module load conda/latest # need to load before making any conda envs
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create --prefix /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 mirdeep2
conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 
```

Success! Okay lets try to run mapping/collapsing step first on one sample. The genome was already indexed through shortstack above so I do not need to do it again. If I did, I would use bowtie (NOT bowtie2) to index the genome. 

```
salloc --mem=30g
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 
module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 
mapper.pl apoc_2_S31_L001_R1_001_trim.fastq -e -h -m -s apoc_2_processed.fa -t apoc_2_mappings.arf -p /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta
```

- `-e` input is fastq format 
- `-h` parse to fasta format 
- `-m` collapse reads 
- `-s` prints processed reads to this file 
- `-t` prints read mappings to this file 
- `-p` path to genome 

Success. Look at the log output

```
Log file for this run is in mapper_logs and called mapper.log_1684750
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 6555995  1977598 4578397 30.165  69.835
seq: 6555995    1977598 4578397 30.165  69.835
# reads processed: 1143802
# reads with at least one alignment: 150321 (13.14%)
# reads that failed to align: 993481 (86.86%)
# reads with alignments suppressed due to -m: 14433 (1.26%)
Reported 326273 alignments
```

Identify miRNAs with the mirdeep2 script 

```
miRDeep2.pl apoc_2_processed.fa /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta apoc_2_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta none none -t N.vectensis -P
```

Immediately got this error: 

```
Error: problem with /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta
Error in line 79: The identifier
cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
contains white spaces
Please check your file for the following issues:
I.  Sequences are allowed only to comprise characters [ACGTNacgtn].
II. Identifiers are not allowed to have withespaces.
You could run remove_white_space_in_id.pl inputfile > newfile
This will remove everything from the id line after the first whitespace
```

Remove white spaces

```
remove_white_space_in_id.pl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta > /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta
```

Going to submit it as a job. In the scripts folder: `nano mirdeep2_id_apoc_2.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 

echo "Running mirdeep2 miRNA ID on apoc2"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 

miRDeep2.pl apoc_2_processed.fa /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta apoc_2_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for apoc2"
conda deactivate
```

Submitted batch job 42994685. Ran! Hmm very interesting...Identifying more possible miRNAs than shortstack. [Here](https://github.com/JillAshey/Cnidarian_sperm_smRNA/tree/main/output/mirdeep2/apoc/apoc_2) are the results for Apoc 2. According to the results, mirdeep2 identified 36 novel miRNAs with a mirdeep2 score > 0 and 101 known miRNAs (I think). See `result_11_09_2025_t_15_17_45.csv` and `result_11_09_2025_t_15_17_45.html` for more details. For the known miRNAs, lots of my AST2021 miRNAs got hits but some were not classified as miRNAs by mirdeep2 because they had a negative mirdeep2 score. Its difficult to know what cutoff to use for mirdeep2 scores. In most coral papers, they use a cutoff of 10 but in other taxa, they use 0. 

I like this figure from [McCreight et al. 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176596): 

![](https://journals.plos.org/plosone/article/figure/image?size=large&id=10.1371/journal.pone.0176596.g004)

Its a nice illustration of the mirdeep2 scores and number of miRNAs. They also include this in the methods: "miRDeep2 was used to predict novel and identify previously annotated miRNA [66]. Merged, trimmed reads were mapped to their respective genomes using the miRDeep2 mapper.pl module with the following parameters: -c -j -l 18 -m -p -s -t–v. MiRDeep2 was executed with default parameters. When making novel miRNA predictions, miRDeep2’s algorithm accounts for already known miRNAs of the species being analyzed and of any related species. We retrieved a list of known miRNAs from miRBase (release 21) for any of our primates that were in the database (H. sapiens, P. troglodytes, P. paniscus, G. gorilla, P. pygmaeus, M. mulatta), and used all known metazoan miRNAs as our “related species” reference. MiRDeep2 assigned a score from -10 to 10 to each miRNA, with a higher number corresponding to increased likelihood that the putative miRNA is functional. This score is partially determined by the availability of any known miRNA, which would inherently result in lower scores for our primates with no information in miRBase. Because of this, we chose a relaxed score cut-off of 0 and minimum read depth of 3 to include a miRNA in our analyses, with the expectation that false positives would be removed during paralog clustering and alignment."

Hmm...what to do...let's run mapper and mirdeep2 on the other apoc samples. I also moved all the apoc 2 stuff into its own folder: `/scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/apoc_2`. 

```
salloc --mem=30g
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 
module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 

# Apoc 3
mapper.pl apoc_3_S32_L001_R1_001_trim.fastq -e -h -m -s apoc_3_processed.fa -t apoc_3_mappings.arf -p /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta
Log file for this run is in mapper_logs and called mapper.log_3219686
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 5968483  1025679 4942804 17.185  82.815
seq: 5968483    1025679 4942804 17.185  82.815
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 986952
# reads with at least one alignment: 84692 (8.58%)
# reads that failed to align: 902260 (91.42%)
# reads with alignments suppressed due to -m: 5554 (0.56%)
Reported 222079 alignments

# Apoc 4
mapper.pl apoc_4_S33_L001_R1_001_trim.fastq -e -h -m -s apoc_4_processed.fa -t apoc_4_mappings.arf -p /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta
Log file for this run is in mapper_logs and called mapper.log_3219860
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 5005980  701506  4304474 14.013  85.987
seq: 5005980    701506  4304474 14.013  85.987
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 978632
# reads with at least one alignment: 73845 (7.55%)
# reads that failed to align: 904787 (92.45%)
# reads with alignments suppressed due to -m: 5118 (0.52%)
Reported 193356 alignments
```

Do apoc 3 first. In the scripts folder: `nano mirdeep2_id_apoc_3.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 

echo "Running mirdeep2 miRNA ID on apoc3"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 

miRDeep2.pl apoc_3_processed.fa /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta apoc_3_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for apoc3"
conda deactivate
```

Submitted batch job 43001459. While this script runs, rerun shortstack on Nvec samples. Downloaded the Nvec genome from [here](https://simrbase.stowers.org/starletseaanemone) and saved it here: `/scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome`.  

In the scripts folder: `nano shortstack_nvec.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec 

echo "Running shortstack for Nvec samples"

#gunzip *fastq.gz

ShortStack \
--genomefile /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/Nvec200.fasta \
--readfile nvec_1_S25_L001_R1_001_trim.fastq \
nvec_2_S26_L001_R1_001_trim.fastq \
nvec_3_S34_L001_R1_001_trim.fastq \
nvec_4_S35_L001_R1_001_trim.fastq \
--known_miRNAs /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta \
--outdir /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/nvec_shortstack \
--dn_mirna

echo "Nvec shortstack complete"

conda deactivate
```

Submitted batch job 43003911. In the meantime, apoc 3 finished running. Results are [here](https://github.com/JillAshey/Cnidarian_sperm_smRNA/tree/main/output/mirdeep2/apoc/apoc_3). Still a very limited number of miRNAs identified...annoying that mirdeep2 does not maintain the names of the miRNAs between samples. Nvec shortstack finished running as well. Same issue with the Apoc data...only like 5 miRNAs identified. See results [here](https://github.com/JillAshey/Cnidarian_sperm_smRNA/tree/main/output/shortstack/nvec/trimmed_default). 

Is this an extraction or library prep issue? I need to contact Nick for extraction and library prep details. 

Okay talked with Nick! Very helpful. He said the following: 

XXXXXXX

My plan is to move forward with mirdeep2 for now. Let's finish the Apoc samples. I already ran the mapping part for Apoc 4. `nano mirdeep2_id_apoc_3.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 30:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 

echo "Running mirdeep2 miRNA ID on apoc4"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 

miRDeep2.pl apoc_4_processed.fa /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta apoc_4_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for apoc4"
conda deactivate
```

Submitted batch job 46990553. While that is running, map the Nvec and Ahya samples with mirdeep2.

```
salloc --mem=30g
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec 
module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 
```

Map Nvec. As a note, the mapping stats is what shows up on the screen after the mapping runs. The reads processed and what not is from the bowtie.log file (mirdeep2 software uses bowtie as an aligner).

```
# Nvec 1
mapper.pl nvec_1_S25_L001_R1_001_trim.fastq -e -h -m -s nvec_1_processed.fa -t nvec_1_mappings.arf -p /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/Nvec200.fasta
Log file for this run is in mapper_logs and called mapper.log_1217577
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 6011442  356824  5654618 5.936   94.064
seq: 6011442    356824  5654618 5.936   94.064
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 946124
# reads with at least one alignment: 148870 (15.73%)
# reads that failed to align: 797254 (84.27%)
# reads with alignments suppressed due to -m: 71963 (7.61%)
Reported 129175 alignments

# Nvec 2
mapper.pl nvec_2_S26_L001_R1_001_trim.fastq -e -h -m -s nvec_2_processed.fa -t nvec_2_mappings.arf -p /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/Nvec200.fasta
Log file for this run is in mapper_logs and called mapper.log_1217687
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 8054309  257649  7796660 3.199   96.801
seq: 8054309    257649  7796660 3.199   96.801
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 1353805
# reads with at least one alignment: 138973 (10.27%)
# reads that failed to align: 1214832 (89.73%)
# reads with alignments suppressed due to -m: 73206 (5.41%)
Reported 108506 alignments

# Nvec 3
mapper.pl nvec_3_S34_L001_R1_001_trim.fastq -e -h -m -s nvec_3_processed.fa -t nvec_3_mappings.arf -p /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/Nvec200.fasta
Log file for this run is in mapper_logs and called mapper.log_1634874
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 7134137  227350  6906787 3.187   96.813
seq: 7134137    227350  6906787 3.187   96.813
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 1014934
# reads with at least one alignment: 156894 (15.46%)
# reads that failed to align: 858040 (84.54%)
# reads with alignments suppressed due to -m: 68474 (6.75%)
Reported 142027 alignments

# Nvec 4
mapper.pl nvec_4_S35_L001_R1_001_trim.fastq -e -h -m -s nvec_4_processed.fa -t nvec_4_mappings.arf -p /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/Nvec200.fasta
Log file for this run is in mapper_logs and called mapper.log_1634913
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 6236194  217165  6019029 3.482   96.518
seq: 6236194    217165  6019029 3.482   96.518
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 1032159
# reads with at least one alignment: 131481 (12.74%)
# reads that failed to align: 900678 (87.26%)
# reads with alignments suppressed due to -m: 61901 (6.00%)
Reported 109402 alignments
```

Hooray now we can run mirdeep2 on the Nvec samples! `nano mirdeep2_nvec.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 

echo "Running mirdeep2 miRNA ID on nvec1"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec 

miRDeep2.pl nvec_1_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/Nvec200.fasta nvec_1_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for nvec1"
echo "Running mirdeep2 miRNA ID on nvec2"

miRDeep2.pl nvec_2_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/Nvec200.fasta nvec_2_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for nvec2"
echo "Running mirdeep2 miRNA ID on nvec3"

miRDeep2.pl nvec_3_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/Nvec200.fasta nvec_3_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for nvec3"
echo "Running mirdeep2 miRNA ID on nvec4"

miRDeep2.pl nvec_4_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/Nvec200.fasta nvec_4_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for nvec4"

conda deactivate
```

Submitted batch job 47009181. This finished running after about 12 hours. 

Before mapping Ahya, I need to generate the bowtie index for the genome (didn't run shortstack with this genome). Need to run as a job bleh. `nano ahya_bowtie_index.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 10:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Bowtie/1.3.1-GCC-11.3.0

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya
bowtie-build Ahyacinthus.chrsV1.fasta Ahyacinthus.chrsV1
```

Submitted batch job 46996202. Success! Mapping Ahya with mirdeep2.

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya

# Ahya 1
mapper.pl ahya_1_S27_L001_R1_001_trim.fastq -e -h -m -s ahya_1_processed.fa -t ahya_1_mappings.arf -p Ahyacinthus.chrsV1.fasta
mapper.pl ahya_1_S27_L001_R1_001_trim.fastq -e -h -m -s ahya_1_processed.fa -t ahya_1_mappings.arf -p Ahyacinthus.chrsV1
Log file for this run is in mapper_logs and called mapper.log_1635617
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 3854254  0       3854254 0.000   100.000
seq: 3854254    0       3854254 0.000   100.000
```

Weird...nothing mapped at all. Tried mapping ahya 2 as well and didn't work. The bowtie log file says this: 

```
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
Illegal instruction
```

Why do you hate me Ahya??? Should I try a different genome? Is there an indexing issue? IDK!!!!! Come back to this...Maybe let's try running shortstack on Ahya and see what we get, as I did not run any shortstack on this species, just apoc and nvec. `nano shortstack_ahya.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/ShortStack4 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya 

echo "Running shortstack for Ahya samples"

#gunzip *fastq.gz

ShortStack \
--genomefile Ahyacinthus.chrsV1.fasta \
--readfile ahya_1_S27_L001_R1_001_trim.fastq \
ahya_2_S28_L001_R1_001_trim.fastq \
ahya_3_S29_L001_R1_001_trim.fastq \
ahya_4_S30_L001_R1_001_trim.fastq \
--known_miRNAs /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference.fasta \
--outdir /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/output/ahya_shortstack \
--dn_mirna

echo "Ahya shortstack complete"

conda deactivate
```

Submitted batch job 47066266. Interestingly, the Ahya shortstack [results](https://github.com/JillAshey/Cnidarian_sperm_smRNA/tree/main/output/shortstack/ahya/trimmed_default) are much more in line with what I was expecting from corals...it looks like 30ish miRNAs were identified, many of them being prior Acropora miRNAs. Why don't the nvec or apoc have more miRNAs present?? Maybe related to how the samples were preserved? I believe the Ahya sample is from Moorea and was preserved in shield, but I'm not sure how the apoc or nvec samples were preserved (need to ask Ben). Or maybe its just low sequencing depth?? 

Let's try running mirdeep2 on Ahya and compare the results. I know I tried to run it above and it didn't work but maybe the script will use the files created by shortstack and that will work. 

Run the mapper script first. 

```
mapper.pl ahya_1_S27_L001_R1_001_trim.fastq -e -h -m -s ahya_1_processed.fa -t ahya_1_mappings.arf -p /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta
Log file for this run is in mapper_logs and called mapper.log_1419116
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 3854254  1450157 2404097 37.625  62.375
seq: 3854254    1450157 2404097 37.625  62.375
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 1059849
# reads with at least one alignment: 414285 (39.09%)
# reads that failed to align: 645564 (60.91%)
# reads with alignments suppressed due to -m: 134761 (12.72%)
Reported 562591 alignments

mapper.pl ahya_2_S28_L001_R1_001_trim.fastq -e -h -m -s ahya_2_processed.fa -t ahya_2_mappings.arf -p /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta
Log file for this run is in mapper_logs and called mapper.log_1419189
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 2358949  367553  1991396 15.581  84.419
seq: 2358949    367553  1991396 15.581  84.419
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 443141
# reads with at least one alignment: 68451 (15.45%)
# reads that failed to align: 374690 (84.55%)
# reads with alignments suppressed due to -m: 21205 (4.79%)
Reported 92142 alignments

mapper.pl ahya_3_S29_L001_R1_001_trim.fastq -e -h -m -s ahya_3_processed.fa -t ahya_3_mappings.arf -p /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta
Log file for this run is in mapper_logs and called mapper.log_1419249
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 5907632  2059409 3848223 34.860  65.140
seq: 5907632    2059409 3848223 34.860  65.140
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 1344290
# reads with at least one alignment: 562244 (41.82%)
# reads that failed to align: 782046 (58.18%)
# reads with alignments suppressed due to -m: 176750 (13.15%)
Reported 743523 alignments

mapper.pl ahya_4_S30_L001_R1_001_trim.fastq -e -h -m -s ahya_4_processed.fa -t ahya_4_mappings.arf -p /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta
Log file for this run is in mapper_logs and called mapper.log_1419312
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 6235983  2122036 4113947 34.029  65.971
seq: 6235983    2122036 4113947 34.029  65.971
Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
# reads processed: 1130403
# reads with at least one alignment: 418508 (37.02%)
# reads that failed to align: 711895 (62.98%)
# reads with alignments suppressed due to -m: 132251 (11.70%)
Reported 559265 alignments
```

WOW!!!! These alignments are so high!! Way more comparable to data that I've analyzed in the past. Very interesting...the plot thickens or at least gives me some more info on the differences between species. `nano mirdeep2_ahya.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load conda/latest # need to load before making any conda envs

conda activate /work/pi_hputnam_uri_edu/conda/envs/mirdeep2 

echo "Running mirdeep2 miRNA ID on ahya 1"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya 

miRDeep2.pl ahya_1_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta ahya_1_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for ahya 1"
echo "Running mirdeep2 miRNA ID on ahya 2"

miRDeep2.pl ahya_2_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta ahya_2_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for ahya 2"
echo "Running mirdeep2 miRNA ID on ahya 3"

miRDeep2.pl ahya_3_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta ahya_3_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for ahya 3"
echo "Running mirdeep2 miRNA ID on ahya 4"

miRDeep2.pl ahya_4_processed.fa /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta ahya_4_mappings.arf /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/data/miRNA_reference_cleaned.fasta none none -t N.vectensis -P

echo "Mirdeep2 complete for ahya 4"

conda deactivate
```

Submitted batch job 47248281.

Based on the miRNA ID step, I can now make a fasta with the expressed miRNAs for target prediction. For Apoc, only 7 miRNAs were expressed. `nano expressed_miRNAs_apoc.fasta`

```
>Cluster_598
UGGCUUGGAGUUCUAGUGAUUGA
>Cluster_597
ACUGAACUCCAAACAUGAUU
>Cluster_595
UGGCUUGGAGUUCUAGUGAUUGA
>Cluster_2651
CGUUCGCGCUCGCCCCCGA
>Cluster_2637
UCAGAAAGAGAAUGAACACCAGA
>Cluster_2306
ACUUUGAAUCCAGCGAUCCG
>Cluster_1200
ACCCGUAGAUCCGAACUUGUGG
```

Install miranda via conda. 

```
cd /work/pi_hputnam_uri_edu/conda/envs
module load conda/latest # need to load before making any conda envs
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create --prefix /work/pi_hputnam_uri_edu/conda/envs/miranda miranda
conda activate /work/pi_hputnam_uri_edu/conda/envs/miranda 
```

Okay now I ran run miranda for Apoc! `nano miranda_strict_mRNA_apoc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

echo "Apoc target prediction with miranda - targeting mRNA seqs"$(date)

module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/miranda 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc

miranda expressed_miRNAs_apoc.fasta /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata_mrna.fasta -en -20 -strict -out miranda_strict_mRNA_apoc.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" miranda_strict_mRNA_apoc.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" miranda_strict_mRNA_apoc.tab | sort | grep '>' > miranda_strict_mRNA_parsed_apoc.txt

echo "counting number of putative interactions predicted" $(date)
wc -l miranda_strict_mRNA_parsed_apoc.txt

echo "Apoc miranda script complete" $(date)
```

Submitted batch job 48591145. 4174 putative interactions predicted. 

For Nvec, 14 miRNAs were expressed. `nano expressed_miRNAs_nvec.fasta`

```
>Cluster_739
UUCGUCUCGUAACUGCCCUGA
>Cluster_1073
UAAUGUUCCUGCUUGUUCCU
>Cluster_1341
UUGCACAUUACCAAUAUUCUG
>Cluster_1343
UUGCACAUCACCAAUGUUCUGAU
>Cluster_1376
CGAGGUAACUGUUGCAGCA
>Cluster_1430
ACCCGUAGAUCCGAACUUGUGG
>Cluster_1640
UAAAAGUGAUCAGUGGUAAGGU
>Cluster_1954
AACUGAUUUUUGGAAACUGGCU
>Cluster_1959
AACUGAUUUUUGGAAACUGGCU
>Cluster_1961
AACUGAUUUUUGGAAACUGGCU
>Cluster_2470
UGACAGAAGAUAGAAGCGCUG
>Cluster_2518
UCUGACAGAUAAUGAUCCACCGCCU
>Cluster_2678
GAAUUCCUCGUUACUGAAGCUC
>Cluster_2764
AAUAUCAGGAGGACGCAUUCG
```

Run miranda for Nvec! `nano miranda_strict_mRNA_nvec.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

echo "Nvec target prediction with miranda - targeting mRNA seqs"$(date)

module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/miranda 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec

miranda expressed_miRNAs_nvec.fasta /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/NV2g.20240221.transcripts.fa -en -20 -strict -out miranda_strict_mRNA_nvec.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" miranda_strict_mRNA_nvec.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" miranda_strict_mRNA_nvec.tab | sort | grep '>' > miranda_strict_mRNA_parsed_nvec.txt

echo "counting number of putative interactions predicted" $(date)
wc -l miranda_strict_mRNA_parsed_nvec.txt

echo "Nvec miranda script complete" $(date)
```

Submitted batch job 48591620

I also want to predict miRNA interactions with the 3'UTR but these gffs do not have these annotated. Since these are not annotated in my genomes, I will have to predict them myself. Like previous analyses, I will say that 3'UTR is 1000bp to the right of the genes. Look at features:

```
cd /work/pi_hputnam_uri_edu/genomes/Apoc
zgrep -v '^#' apoculata.gff3.gz | cut -s -f 3 | sort | uniq -c | sort -rn > all_features.txt

239708 exon
 236997 CDS
  47156 gene
  45867 mRNA
  44823 stop_codon
  44312 start_codon
   2889 five_prime_UTR
   2317 tRNA
    687 three_prime_UTR
```

Make gff for genes only 

```
zgrep $'\tgene\t' apoculata.gff3.gz > apoc_GFFannotation.gene.gff
```

Extract scaffold names and lengths

```
cat is apoculata.genome.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > apoc.Chromosome_lengths.txt

# Sort so that chromosomes are in order 
sort -k1,1 -V apoc.Chromosome_lengths.txt > apoc.Chromosome_lengths_sorted.txt
```

Extract scaffold names 

```
awk -F" " '{print $1}' apoc.Chromosome_lengths_sorted.txt > apoc.Chromosome_names.txt
```

Sort gene GFF by chromosome names 

```
module load bedtools2/2.31.1
sortBed -faidx apoc.Chromosome_names.txt -i apoc_GFFannotation.gene.gff > apoc_GFFannotation.gene_sorted.gff
```

Extract 1000bp flank around genes and subtract any gene overlap

```
# Flank
flankBed -i apoc_GFFannotation.gene_sorted.gff -g apoc.Chromosome_lengths_sorted.txt -l 0 -r 1000 -s | awk '{gsub("gene","3prime_UTR",$3); print $0 }' | awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | tr ' ' '\t' > apoc.GFFannotation.3UTR_1kb.gff

# Subtract 
subtractBed -a apoc.GFFannotation.3UTR_1kb.gff -b apoc_GFFannotation.gene_sorted.gff > apoc.GFFannotation.3UTR_1kb_corrected.gff 

# Sort 
sortBed -faidx apoc.Chromosome_names.txt -i apoc.GFFannotation.3UTR_1kb_corrected.gff > apoc.GFFannotation.3UTR_1kb_corrected.sorted.gff 
```

Look for closest 3UTR to each gene

```
bedtools closest -a apoc_GFFannotation.gene_sorted.gff \
                   -b apoc.GFFannotation.3UTR_1kb_corrected.sorted.gff \
                   -s -D a \
                   -t first > apoc_closest_UTRs.txt
```

Make 3UTR gff 

```
awk '{OFS="\t"; print $10, $11, $12, $13, $14, $15, $16, $17, $18 }' apoc_closest_UTRs.txt > apoc_3UTRs_closest.gff
```

Make fasta of 3UTRs

```
bedtools getfasta -fi apoculata.genome.fasta -bed apoc_3UTRs_closest.gff -fo apoc_3UTRs_closest.fasta -name -s
```

Hooray! Now run miranda with 3UTRs. `nano miranda_strict_3UTR_apoc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

echo "Apoc target prediction with miranda - targeting 3UTR seqs"$(date)

module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/miranda 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc

miranda expressed_miRNAs_apoc.fasta /work/pi_hputnam_uri_edu/genomes/Apoc/apoc_3UTRs_closest.fasta -en -20 -strict -out miranda_strict_3UTR_apoc.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" miranda_strict_3UTR_apoc.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" miranda_strict_3UTR_apoc.tab | sort | grep '>' > miranda_strict_3UTR_parsed_apoc.txt

echo "counting number of putative interactions predicted" $(date)
wc -l miranda_strict_3UTR_parsed_apoc.txt

echo "Apoc miranda script complete" $(date)
```

Submitted batch job 48632520

Do the same 3UTR identification with Nvec. 

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome
grep -v '^#' NV2g.20240221.gff | cut -s -f 3 | sort | uniq -c | sort -rn > all_features.txt

287079 exon
 241911 CDS
  35881 mRNA
  24526 gene
```

Make gff for mRNA only 

```
grep $'\tmRNA\t' NV2g.20240221.gff > nvec_GFFannotation.mRNA.gff
```

Extract scaffold names and lengths

```
cat is Nvec200.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > nvec.Chromosome_lengths.txt
```

Extract scaffold names 

```
awk -F" " '{print $1}' nvec.Chromosome_lengths.txt > nvec.Chromosome_names.txt
```

Sort mRNA GFF by chromosome names 

```
module load bedtools2/2.31.1
sortBed -faidx nvec.Chromosome_names.txt -i nvec_GFFannotation.mRNA.gff > nvec_GFFannotation.mRNA_sorted.gff
```

Extract 1000bp flank around genes and subtract any mRNA overlap

```
# Flank
flankBed -i nvec_GFFannotation.mRNA_sorted.gff -g nvec.Chromosome_lengths.txt -l 0 -r 1000 -s | awk '{gsub("mRNA","3prime_UTR",$3); print $0 }' | awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | tr ' ' '\t' > nvec.GFFannotation.3UTR_1kb.gff

# Subtract 
subtractBed -a nvec.GFFannotation.3UTR_1kb.gff -b nvec_GFFannotation.mRNA_sorted.gff > nvec.GFFannotation.3UTR_1kb_corrected.gff 

# Sort 
sortBed -faidx nvec.Chromosome_names.txt -i nvec.GFFannotation.3UTR_1kb_corrected.gff > nvec.GFFannotation.3UTR_1kb_corrected.sorted.gff 
```

Make fasta of 3UTRs

```
bedtools getfasta -fi Nvec200.fasta -bed nvec.GFFannotation.3UTR_1kb_corrected.sorted.gff -s -name > nvec_3UTRs.fasta
```

Run miranda with 3UTRs. `nano miranda_strict_3UTR_nvec.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

echo "NVec target prediction with miranda - targeting 3UTR seqs"$(date)

module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/miranda 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec

miranda expressed_miRNAs_nvec.fasta /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/nvec_3UTRs.fasta -en -20 -strict -out miranda_strict_3UTR_nvec.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" miranda_strict_3UTR_nvec.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" miranda_strict_3UTR_nvec.tab | sort | grep '>' > miranda_strict_3UTR_parsed_nvec.txt

echo "counting number of putative interactions predicted" $(date)
wc -l miranda_strict_3UTR_parsed_nvec.txt

echo "Apoc miranda script complete" $(date)
```

Submitted batch job 48634665

For Ahya, 27 miRNAs were expressed. `nano expressed_miRNAs_ahya.fasta`

```
>Cluster_1332
AUUGAUUGUAGACAAGCCU
>Cluster_1356
UCCCGUAGAUCCGAACUUGUGG
>Cluster_1518
UCUGCGUUAUCGGUGAAAUUGU
>Cluster_1641
CAAGUGAGAGAAGGUUAGUGUGG
>Cluster_1975
UAUGGGUUGACAGUCGACGGUC
>Cluster_3756
AGUGCACUUUUCUCAGGAUG
>Cluster_3818
UAAUGUUCGCAACUGCCUUGU
>Cluster_7759
UAUGAGACCGUUUGCCAAUGUU
>Cluster_7774
UCUGGCAGUAUGUUAUUUUUCC
>Cluster_8901
UUAACGAGUAGAUAAAUGAAGAGU
>Cluster_9540
UAGGCGUAUUUCCGAUUGUCCU
>Cluster_10318
ACUGCAGCUAAAUACUCCGCUG
>Cluster_10319
AAAAAUGUCGGUUGCUUAAGCU
>Cluster_10493
CAAUGUUUCGGCUUGUUCCCG
>Cluster_11616
AAGAACACCCAAAAUAGCUGAGGA
>Cluster_11958
UUGGUAGAGAGAAAUGACAAAAU
>Cluster_12579
UCUGCCAAUCGUCAGACAAACU
>Cluster_12617
UUUUUGUGAUGUUCGUCAAUA
>Cluster_13735
CAAGGAGGAAGCAUGAUACGUA
>Cluster_16517
UUCAAUGCAAAACAAGUCAGGUAU
>Cluster_16534
UCAUAACAGUGAGGACCAUUCU
>Cluster_16886
AUUCGUCUCUUGUAUUUUCCUGG
>Cluster_18389
ACUUGAAUGCGAGGCUCAUUAGU
>Cluster_18872
UCUCUUACAUUGUUGUGCUGCC
>Cluster_21323
AAAAAUUUCGUUUCAGGGCCA
```

Run miranda with mRNA prediction. `nano miranda_strict_mRNA_ahya.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

echo "Ahya target prediction with miranda - targeting mRNA seqs"$(date)

module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/miranda 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya

miranda expressed_miRNAs_ahya.fasta Ahyacinthus.transcripts.fasta -en -20 -strict -out miranda_strict_mRNA_ahya.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" miranda_strict_mRNA_ahya.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" miranda_strict_mRNA_ahya.tab | sort | grep '>' > miranda_strict_mRNA_parsed_ahya.txt

echo "counting number of putative interactions predicted" $(date)
wc -l miranda_strict_mRNA_parsed_ahya.txt

echo "Ahya miranda script complete" $(date)
```

Submitted batch job 48680891





### piRNAs

Moving to identifying piRNAs from the data, given that there seems to be a high prevalence of piRNAs in the data. Following methods outlined in Ashey et al. 2025. Using code from [deep dive repo](https://github.com/urol-e5/deep-dive/tree/main). 

First, remove redundant reads from trimmed fastqs, keep reads between 25 and 35 nt, remove low complexity reads and change the format for SeqMap alignment. 

In the scripts folder: `nano piRNA_prep.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Apoc piRNA read prep"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/

# Apoc 
for f in /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/*.fastq
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_collapse.pl -i ${f} -o ${f}.collapsed
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_length-filter.pl -i ${f}.collapsed -o ${f}.collapsed.filt -min 25 -max 35
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_duster.pl -i ${f}.collapsed.filt
done

echo "Nvec piRNA read prep"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/

# Nvec 
for f in /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/*.fastq
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_collapse.pl -i ${f} -o ${f}.collapsed
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_length-filter.pl -i ${f}.collapsed -o ${f}.collapsed.filt -min 25 -max 35
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_duster.pl -i ${f}.collapsed.filt
done

echo "Ahya piRNA read prep"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/

# Nvec 
for f in /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/*.fastq
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_collapse.pl -i ${f} -o ${f}.collapsed
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_length-filter.pl -i ${f}.collapsed -o ${f}.collapsed.filt -min 25 -max 35
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_duster.pl -i ${f}.collapsed.filt
done
```

Submitted batch job 43327763

While this is running, install [tRNAscan](https://github.com/UCSC-LoweLab/tRNAscan-SE) and [barnap](https://github.com/tseemann/barrnap) via conda. 

```
cd /work/pi_hputnam_uri_edu/conda/envs
module load conda/latest # need to load before making any conda envs
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# trnascan-se
conda create --prefix /work/pi_hputnam_uri_edu/conda/envs/trnascan trnascan-se
conda activate /work/pi_hputnam_uri_edu/conda/envs/trnascan
tRNAscan-SE -h
conda deactivate 

# barrnap
conda create --prefix /work/pi_hputnam_uri_edu/conda/envs/barrnap barrnap
conda activate /work/pi_hputnam_uri_edu/conda/envs/barrnap
barrnap -h
conda deactivate 
```

Run tRNAscan and barrnap on all genomes. In the scripts folder: `nano tRNAscan_barrnap.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/trnascan

echo "Starting apoc tRNA ID"
echo "Apoc tRNAscan"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 

tRNAscan-SE \
-E \
-o Apoc-tRNA.out \
-f Apoc-tRNA_struct.out \
-m Apoc-tRNA_stats.out \
-b Apoc-tRNA.bed \
-j Apoc-tRNA.gff3 \
-a Apoc-tRNA.fasta \
-d \
/work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta 

echo "Apoc tRNA ID complete, Starting nvec tRNA ID"
echo "Nvec tRNAscan"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec

tRNAscan-SE \
-E \
-o Nvec-tRNA.out \
-f Nvec-tRNA_struct.out \
-m Nvec-tRNA_stats.out \
-b Nvec-tRNA.bed \
-j Nvec-tRNA.gff3 \
-a Nvec-tRNA.fasta \
-d \
/work/pi_hputnam_uri_edu/genomes/Nvec/Nvec200.fasta 

echo "Nvec tRNA ID complete, Starting ahya tRNA ID"
echo "Ahya tRNAscan"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya

tRNAscan-SE \
-E \
-o Ahya-tRNA.out \
-f Ahya-tRNA_struct.out \
-m Ahya-tRNA_stats.out \
-b Ahya-tRNA.bed \
-j Ahya-tRNA.gff3 \
-a Ahya-tRNA.fasta \
-d \
/work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta

echo "Ahya tRNA ID complete, all spp done"

conda deactivate
conda activate /work/pi_hputnam_uri_edu/conda/envs/barrnap

echo "Starting apoc rRNA ID"
echo "Apoc barrnap"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc 

barrnap \
--kingdom euk \
--outseq Apoc-rRNA.fasta \
/work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta 

echo "Apoc rRNA ID complete, Starting nvec rRNA ID"
echo "Nvec barrnap"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec

barrnap \
--kingdom euk \
--outseq Nvec-rRNA.fasta \
/work/pi_hputnam_uri_edu/genomes/Nvec/Nvec200.fasta 

echo "Nvec rRNA ID complete, Starting ahya rRNA ID"
echo "Ahya barrnap"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya

barrnap \
--kingdom euk \
--outseq Ahya-rRNA.fasta \
/work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta

echo "Ahya rRNA ID complete, all spp done"

conda deactivate
```

Submitted batch job 43338150. Check number of tRNAs and rRNAs identified, and cat together. 

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc
grep -c ">" Apoc-rRNA.fasta 
228
grep -c ">" Apoc-tRNA.fasta 
11196
cat Apoc-rRNA.fasta Apoc-tRNA.fasta > Apoc_tRNA_rRNA.fasta 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec
grep -c ">" Nvec-rRNA.fasta 
1097
grep -c ">" Nvec-tRNA.fasta 
9910
cat Nvec-rRNA.fasta Nvec-tRNA.fasta > Nvec_tRNA_rRNA.fasta  

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya
grep -c ">" Ahya-rRNA.fasta 
1738
grep -c ">" Ahya-tRNA.fasta 
8629
cat Ahya-rRNA.fasta Ahya-tRNA.fasta > Ahya_tRNA_rRNA.fasta
```

Install [sortmerna](https://github.com/sortmerna/sortmerna) to clear tRNA and rRNA reads 

```
cd /work/pi_hputnam_uri_edu/conda/envs
module load conda/latest # need to load before making any conda envs
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create --prefix /work/pi_hputnam_uri_edu/conda/envs/sortmerna sortmerna
conda activate /work/pi_hputnam_uri_edu/conda/envs/sortmerna
sortmerna -h
conda deactivate 
```

Run sortmerna. In the scripts folder: `nano sortmerna.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/sortmerna

echo "Starting apoc sortmerna"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/

for f in *no-dust
do
mkdir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/${f/-collapsed.filt.no-dust}
sortmerna \
--ref Apoc_tRNA_rRNA.fasta \
--reads ${f} \
--workdir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/${f/-collapsed.filt.no-dust} \
--fastx \
--other
done

echo "Apoc sortmerna complete, starting nvec sortmerna"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/

for f in *no-dust
do
mkdir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/${f/-collapsed.filt.no-dust}
sortmerna \
--ref Nvec_tRNA_rRNA.fasta \
--reads ${f} \
--workdir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/${f/-collapsed.filt.no-dust} \
--fastx \
--other
done

echo "Nvec sortmerna complete, starting ahya sortmerna"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/

for f in *no-dust
do
mkdir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/${f/-collapsed.filt.no-dust}
sortmerna \
--ref Ahya_tRNA_rRNA.fasta \
--reads ${f} \
--workdir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/${f/-collapsed.filt.no-dust} \
--fastx \
--other
done

echo "Ahya sortmerna complete"
conda deactivate 
```

Submitted batch job 43348919. Success. The sequences to move forward with are stored in the file `other.fq` under each sample folder created. Okay for some reason the apoc3 fastq file was in a different folder and did not get run with all the others. Write a script to do all the previous things for apoc3. `nano apoc3_piRNA_prep.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Apoc 3 piRNA read prep"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/

perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_collapse.pl -i apoc_3_S32_L001_R1_001_trim.fastq -o apoc_3_S32_L001_R1_001_trim.fastq.collapsed

perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_length-filter.pl -i apoc_3_S32_L001_R1_001_trim.fastq.collapsed -o apoc_3_S32_L001_R1_001_trim.fastq.collapsed.filt -min 25 -max 35

perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_duster.pl -i apoc_3_S32_L001_R1_001_trim.fastq.collapsed.filt

module load conda/latest # need to load before making any conda envs
conda activate /work/pi_hputnam_uri_edu/conda/envs/sortmerna

echo "Starting apoc 3 sortmerna"

for f in apoc_3*no-dust
do
mkdir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/${f/-collapsed.filt.no-dust}
sortmerna \
--ref Apoc_tRNA_rRNA.fasta \
--reads ${f} \
--workdir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/${f/-collapsed.filt.no-dust} \
--fastx \
--other
done

echo "Apoc 3 sortmerna complete"
```

Submitted batch job 43350906. Success. Time to run sRNAmapper. In the scripts folder: `nano sRNAmapper.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting apoc putative piRNA mapping"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna

for f in apoc_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/sRNAmapper.pl \
-input ${f}/out/other.fq \
-genome /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta \
-alignments best 
done

echo "Apoc putative piRNA mapping complete, starting nvec putative piRNA mapping "

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna

for f in nvec_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/sRNAmapper.pl \
-input ${f}/out/other.fq \
-genome /work/pi_hputnam_uri_edu/genomes/Nvec/Nvec200.fasta \
-alignments best 
done

echo "Nvec putative piRNA mapping complete, starting ahya putative piRNA mapping "

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna

for f in ahya_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/sRNAmapper.pl \
-input ${f}/out/other.fq \
-genome /work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta \
-alignments best 
done

echo "Ahya putative piRNA mapping complete"
```

Submitted batch job 43350967. While this runs, run repeatmasker on the genomes. In the scripts folder: `nano repeatmasker.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load RepeatMasker/4.1.5-foss-2022a

echo "Apoc repeatmasker"

RepeatMasker -norna -dir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/ /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta

echo "Nvec repeatmasker"

RepeatMasker -norna -dir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/ /work/pi_hputnam_uri_edu/genomes/Nvec/Nvec200.fasta

echo "Ahya repeatmasker"

RepeatMasker -norna -dir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/ /work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta
```

Submitted batch job 43362252. 

Extended scratch workspace on 9/25/25 for 30 days. My `43350967` job timed out, it takes a lot longer to align the reads with sRNAmapper than I anticipated. It ended in the middle of nvec 3, so I need to run Nvec 3 and 4 and Ahya 1, 2, 3, and 4. I am going to run these all as separate jobs. 

`nano sRNAmapper_nvec_3.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/nvec_3_S34_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

echo "Starting nvec 3 putative piRNA mapping "

perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/sRNAmapper.pl \
-input other.fq \
-genome /work/pi_hputnam_uri_edu/genomes/Nvec/Nvec200.fasta \
-alignments best

echo "Nvec 3 putative piRNA mapping complete"
```

Submitted batch job 44107785. This stopped running -- I needed to remove the 'done'. Submitted batch job 47008570. Waiting to see how this one does before starting the others. 

If it errors out again with time, try: 

```
#SBATCH -t 4-00:00:00
#SBATCH -q long
```

Come back to nvec 4 and the ahya samples after nvec 3 sRNA mapper runs. Okay Nvec 3 ran in about 2 days. I'm going to start nvec 4 and ahya 1. `nano sRNAmapper_nvec_4.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 5-00:00:00
#SBATCH -q long
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/nvec_4_S35_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

echo "Starting nvec 4 putative piRNA mapping "

perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/sRNAmapper.pl \
-input other.fq \
-genome /work/pi_hputnam_uri_edu/genomes/Nvec/Nvec200.fasta \
-alignments best

echo "Nvec 4 putative piRNA mapping complete"
```

Submitted batch job 47247790. This ran successfully in about 2 days. 

Run ahya 1 as well to see how long it will take to run those samples. `nano sRNAmapper_ahya_1.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 5-00:00:00
#SBATCH -q long
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_1_S27_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

echo "Starting ahya 1 putative piRNA mapping"

perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/sRNAmapper.pl \
-input other.fq \
-genome /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta \
-alignments best 

echo "Ahya 1 putative piRNA mapping complete"
```

Submitted batch job 47887852

`nano sRNAmapper_ahya_2.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 5-00:00:00
#SBATCH -q long
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_2_S28_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

echo "Starting ahya 2 putative piRNA mapping"

perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/sRNAmapper.pl \
-input other.fq \
-genome /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta \
-alignments best 

echo "Ahya 2 putative piRNA mapping complete"
```

Submitted batch job 47888085

`nano sRNAmapper_ahya_3.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 5-00:00:00
#SBATCH -q long
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_3_S29_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

echo "Starting ahya 3 putative piRNA mapping"

perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/sRNAmapper.pl \
-input other.fq \
-genome /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta \
-alignments best

echo "Ahya 3 putative piRNA mapping complete"
```

Submitted batch job 47888099

`nano sRNAmapper_ahya_4.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 5-00:00:00
#SBATCH -q long
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_4_S30_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

echo "Starting ahya 4 putative piRNA mapping"

perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/sRNAmapper.pl \
-input other.fq \
-genome /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta \
-alignments best

echo "Ahya 4 putative piRNA mapping complete"
```

Submitted batch job 47888107. 

Move forward with the Apoc piRNA data. After sRNA mapping, I will use the reallocate script for origin assignment for reads with multiple mappings (following Ashey et al. 2025 github [code](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/18-PEVE-piRNA-proTRAC.Rmd) from Javi). `nano reallocate_apoc.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting apoc origin assessment for reads with multiple mapping locations"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna

for f in apoc_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/reallocate.pl ${f}/out/other.fq.map 10000 1000 b 0
done

echo "Apoc origin assessment complete"
```

Submitted batch job 47068517. Complete after a few hours! Time to run proTRAC, which predicts piRNA clusters in each sample. `nano protrac_apoc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting protrac for apoc"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna

for f in apoc_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/proTRAC_2.4.2.pl \
-map ${f}/out/other.fq.map.weighted-10000-1000-b-0 \
-genome /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta \
-repeatmasker /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/apoculata.genome.fasta.out \
-geneset /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.gtf
done

echo "Apoc protrac complete!"
```

Submitted batch job 47074553. Success! Looking at the slurm output file:

```
## Apoc 2 -- proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m21d17h2m44s
Predicted Cluster: Ap14 14605235-14607384 directionality: mono:plus
Predicted Cluster: Ap12 7323688-7324994 directionality: mono:plus
Predicted Cluster: Ap9 26361024-26368010 directionality: mono:plus
Predicted Cluster: Ap7 6545108-6552425 directionality: mono:plus
Predicted Cluster: Ap6 22072131-22076641 directionality: mono:minus
Predicted Cluster: Ap13 5529006-5531171 directionality: mono:minus
Predicted Cluster: Ap2 28449661-28453270 directionality: mono:plus
Predicted Cluster: Ap2 28609004-28613642 directionality: mono:minus
Predicted Cluster: Ap1 11987043-11991524 directionality: mono:minus
Predicted Cluster: Ap1 32291280-32293852 directionality: mono:minus
Predicted Cluster: Ap4 8249079-8254709 directionality: mono:minus
Predicted Cluster: Ap4 14674394-14681635 directionality: mono:minus
Predicted Cluster: Ap4 15216264-15217852 directionality: mono:minus
Predicted Cluster: Ap4 16367577-16368633 directionality: mono:plus
Predicted Cluster: Ap4 16429408-16431300 directionality: mono:minus
Predicted Cluster: Ap4 19551586-19556522 directionality: mono:plus
Predicted Cluster: Ap4 20112688-20114149 directionality: mono:minus
Predicted Cluster: Ap4 21144289-21149927 directionality: mono:plus
Predicted Cluster: Ap4 22095192-22098356 directionality: mono:minus
Predicted Cluster: Ap4 23905724-23909328 directionality: mono:minus
Predicted Cluster: Ap4 23967567-23975496 directionality: mono:plus
Predicted Cluster: Ap4 24418700-24423869 directionality: mono:minus
Predicted Cluster: Ap4 24817260-24822015 directionality: mono:minus
Predicted Cluster: Ap4 25548882-25551839 directionality: mono:plus
Predicted Cluster: Ap4 26247700-26249983 directionality: mono:plus
Predicted Cluster: Ap4 28311087-28315846 directionality: mono:plus
Predicted Cluster: Ap4 28359164-28360486 directionality: mono:plus
Predicted Cluster: Ap4 28792033-28806012 directionality: mono:plus
Predicted Cluster: Ap4 29339861-29348771 directionality: mono:plus
Predicted Cluster: Ap4 29847190-29856948 directionality: mono:minus
Predicted Cluster: Ap4 29921571-29923015 directionality: mono:plus
Predicted Cluster: Ap4 30057962-30062685 directionality: mono:minus
Predicted Cluster: Ap4 32416234-32418457 directionality: mono:minus
Predicted Cluster: Ap4 32616138-32620445 directionality: mono:plus
Predicted Cluster: Ap4 35477812-35481957 directionality: mono:minus
Total size of 35 predicted piRNA clusters: 154641 bp (0.034%)

## Apoc 3 - proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m21d17h3m59s
Predicted Cluster: Ap14 22536285-22543883 directionality: mono:plus
Predicted Cluster: Ap4 23914752-23918489 directionality: mono:minus
Predicted Cluster: Ap4 24422101-24423714 directionality: mono:minus
Predicted Cluster: Ap4 25549401-25551567 directionality: mono:plus
Predicted Cluster: Ap4 28288220-28291892 directionality: mono:minus
Predicted Cluster: Ap4 28792034-28796894 directionality: mono:plus
Predicted Cluster: Ap4 28802423-28804635 directionality: mono:plus
Predicted Cluster: Ap4 29340005-29345838 directionality: mono:plus
Predicted Cluster: Ap4 29848336-29853006 directionality: mono:minus
Total size of 9 predicted piRNA clusters: 36370 bp (0.008%)

## Apoc 4 - proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m21d17h5m13s
Predicted Cluster: Ap7 27892450-27894310 directionality: mono:minus
Predicted Cluster: Ap2 28611868-28612921 directionality: mono:minus
Predicted Cluster: Ap4 14677878-14679887 directionality: mono:minus
Predicted Cluster: Ap4 19821254-19827332 directionality: mono:minus
Predicted Cluster: Ap4 24422101-24423714 directionality: mono:minus
Predicted Cluster: Ap4 24809883-24811688 directionality: mono:minus
Predicted Cluster: Ap4 28288593-28290523 directionality: mono:minus
Predicted Cluster: Ap4 28311357-28313299 directionality: mono:plus
Predicted Cluster: Ap4 28796093-28803974 directionality: mono:plus
Predicted Cluster: Ap4 29851252-29853436 directionality: mono:minus
Predicted Cluster: Ap4 30061453-30062633 directionality: mono:minus
Predicted Cluster: Ap4 32416462-32418457 directionality: mono:minus
Total size of 12 predicted piRNA clusters: 31542 bp (0.007%)
```

Interesting that there are way more clusters in Apoc 2...what if I ran protrac just on the mapped reads without reallocating? Changed it in the script and reran. Submitted batch job 47077635. Very similar results, slightly more clusters with the reallocated ones. Run ping pong ID on all reads. `nano pingpong_apoc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting ping pong for apoc"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna

for f in apoc_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_pingpong.pl \
-i ${f}/out/other.fq.map.weighted-10000-1000-b-0 \
-o ${f}.pp
done

echo "Apoc ping pong complete!"
```

Submitted batch job 47083778

Let's look at the cluster overlap between the 3 replicates. First, create bed files from cluster fasta files.

Apoc 2: 

```
salloc 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m21d17h2m44s

input_fasta="clusters.fasta"
output_bed="clusters_apoc_2.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Apoc 3:

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m21d17h3m59s

input_fasta="clusters.fasta"
output_bed="clusters_apoc_3.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Apoc 4

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m21d17h5m13s

input_fasta="clusters.fasta"
output_bed="clusters_apoc_4.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Move all cluster bed files to same folder 

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna
mkdir proTRAC_bed
cd proTRAC_bed

mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m21d17h2m44s/clusters_apoc_2.bed .
mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m21d17h3m59s/clusters_apoc_3.bed .
mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m21d17h5m13s/clusters_apoc_4.bed .
```

Merge clusters across samples using bedtools. Following Javi's work, I selected clusters that are present in at least two reps and overlap is at least 50% of the cluster and is reciprocal (ie A overlaps 50% of B and B overlaps 50% of A)

```
FILES=(*.bed)
NUM_FILES=${#FILES[@]}
# Loop over each file
for (( i=0; i<$NUM_FILES; i++ )); do
  for (( j=i+1; j<$NUM_FILES; j++ )); do
    FILE1=${FILES[$i]}
    FILE2=${FILES[$j]}
    
    intersectBed -a $FILE1 -b $FILE2 -f 0.5 -r > ${i}_${j}_merged.bed
  done
done

cat *_merged.bed | sort -k1,1 -k2,2n | bedtools merge -i - > apoc.merged.clusters.bed
```

Here are the merged clusters:

```
Ap4     24422100        24423714
Ap4     25549400        25551567
Ap4     28288592        28290523
Ap4     28796092        28803974
Ap4     29340004        29345838
Ap4     32416461        32418457
```

Very odd that they are all present in chromosome 4 within a couple million base pairs of each other...

The ping pong analysis completed! Here's an example of the results file: 

```
Results for apoc_3_S32_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out/other.fq.map.weighted-10000-1000-b-0
Total scores for 5' overlap:
overlap [nt]	read pairs
1	1.5803821756036e-05
2	
3	
4	
5	
6	
7	
8	
9	
10	8.00126787427719
11	
12	
13	
14	
15	
16	
17	
18	
19	1
20	
21	
22	2
23	
24	
25	
26	
27	
28	
29	
30	
31	
32	
33	
34	

Z-score calculation:
Average background (overlaps 1-9+11-20): 0.0526324107274609
Variance background: 0.0498614083015296
Standard deviation: 0.223296682244787
Ping-Pong Z-Score: 35.5967468197136

Z-score >= 1.6449 (one-tailed hypothesis) -> p < 0.05
Z-score >= 2.3264 (one-tailed hypothesis) -> p < 0.01
```

The output is essentially overlap scores between complementary reads -- how often a piRNA 5' end overlaps another piRNA 5' end by a given number of nucleotides. The important value is at 10 nt overlap, as the sign of ping pong amplification is 10 nt overlap where the secondary piRNA starts 10 nt down from the primary piRNA. It's all quite confusing to me. But we see that there is overlap at 10 nt, which is more than any of the other nucleotides and suggests active ping-pong processing. The z score is quantifiying how strong the 10 nt overlap is compared to background overlaps. We see the average background is quite low, while ping pong is quite high, which further supports the ping pong amplification processing here. This provides more evidence that these reads/clusters are bona fide piRNAs.  

Extended scratch workspace for 30 more days on 10/25/25. 

```
ws_extend cnidarian_sperm 30
Info: could not read email from users config ~/.ws_user.conf.
Info: reminder email will be sent to local user account
Info: extending workspace.
Info: changed mail address to jillashey_uri_edu
Info: changed reminder setting.
/scratch3/workspace/jillashey_uri_edu-cnidarian_sperm
remaining extensions  : 3
remaining time in days: 30
```

Hooray everything has finally finshed mapping!! After sRNA mapping, I will use the reallocate script for origin assignment for reads with multiple mappings for both Ahya and Nvec. 

`nano reallocate_nvec.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting nvec origin assessment for reads with multiple mapping locations"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna

for f in nvec_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/reallocate.pl ${f}/out/other.fq.map 10000 1000 b 0
done

echo "nvec origin assessment complete"
```

Submitted batch job 48116373

`nano reallocate_ahya.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting ahya origin assessment for reads with multiple mapping locations"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna

for f in ahya_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/reallocate.pl ${f}/out/other.fq.map 10000 1000 b 0
done

echo "ahya origin assessment complete"
```

Submitted batch job 48116386. 

These finished running, now time to run protrac! `nano protrac_nvec.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting protrac for nvec"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna

for f in nvec_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/proTRAC_2.4.2.pl \
-map ${f}/out/other.fq.map.weighted-10000-1000-b-0 \
-genome /work/pi_hputnam_uri_edu/genomes/Nvec/Nvec200.fasta \
-repeatmasker /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/Nvec200.fasta.out \
-geneset /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/genome/NV2g.20240221.gff
done

echo "Nvec protrac complete!"
```

Submitted batch job 48117305. Success! Look at slurm output file. 

```
## Nvec 1 -- proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h5m50s
Predicted Cluster: chr1 17237284-17238894 directionality: mono:plus
Predicted Cluster: chr1 20737155-20739413 directionality: mono:plus
Predicted Cluster: chr2 19149473-19151911 directionality: mono:minus
Predicted Cluster: chr4 11413352-11417730 directionality: mono:plus
Predicted Cluster: chr4 12785446-12787930 directionality: mono:minus
Predicted Cluster: chr7 9542183-9552840 directionality: mono:minus
Predicted Cluster: chr7 11470990-11472872 directionality: mono:minus
Predicted Cluster: chr8 13768282-13774843 directionality: mono:plus
Predicted Cluster: chr9 4104018-4108597 directionality: mono:minus
Predicted Cluster: chr11 8498144-8508380 directionality: mono:plus
Predicted Cluster: chr12 1380645-1385812 directionality: mono:plus
Predicted Cluster: chr12 6761627-6764087 directionality: mono:minus
Predicted Cluster: chr12 6786479-6789048 directionality: mono:minus

## Nvec 2 -- 	proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h7m6s
Predicted Cluster: chr4 12291022-12294057 directionality: mono:minus
Predicted Cluster: chr4 12785446-12787935 directionality: mono:minus
Predicted Cluster: chr11 6084013-6091893 directionality: mono:minus
Predicted Cluster: chr11 8499046-8510667 directionality: mono:plus
Predicted Cluster: chr12 6786480-6789047 directionality: mono:minus
Predicted Cluster: chr13 7916562-7920742 directionality: mono:plus

## Nvec 3 -- proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h8m22s
Predicted Cluster: chr1 21551145-21561022 directionality: mono:minus
Predicted Cluster: chr8 858000-863372 directionality: mono:minus
Predicted Cluster: chr9 2242772-2251021 directionality: mono:minus
Predicted Cluster: chr9 9295023-9299988 directionality: mono:plus
Predicted Cluster: chr11 6080374-6097935 directionality: mono:minus
Predicted Cluster: chr11 6111102-6119625 directionality: mono:minus
Predicted Cluster: chr11 8499924-8511124 directionality: mono:plus
Predicted Cluster: chr12 1364983-1369018 directionality: mono:plus
Predicted Cluster: chr12 1378283-1386019 directionality: mono:plus
Predicted Cluster: chr13 7392022-7396886 directionality: mono:minus
Predicted Cluster: chr13 7423036-7439952 directionality: mono:minus
Predicted Cluster: chr13 9246461-9247755 directionality: mono:plus
Predicted Cluster: chr13 11846002-11854946 directionality: mono:plus
Predicted Cluster: chr13 11879268-11887968 directionality: mono:plus
Predicted Cluster: chr15 1521213-1535022 directionality: mono:minus
Predicted Cluster: chr15 1561014-1568949 directionality: mono:minus

## Nvec 4 -- proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h9m37s
Predicted Cluster: chr2 19148529-19151000 directionality: mono:minus
Predicted Cluster: chr11 8499926-8508480 directionality: mono:plus
Predicted Cluster: chr12 1379417-1382857 directionality: mono:plus
Predicted Cluster: chr12 12443567-12452679 directionality: mono:minus
```

Similar numbers of clusters to Apoc. Run ping pong ID on all reads. `nano pingpong_nvec.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting ping pong for nvec"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna

for f in nvec_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_pingpong.pl \
-i ${f}/out/other.fq.map.weighted-10000-1000-b-0 \
-o ${f}.pp
done

echo "Nvec ping pong complete!"
```

Submitted batch job 48117451. Looking at the ping pong output, XXXXX. 

Let's look at the cluster overlap between the 4 reps. Create bed files from cluster fasta files.

Nvec 1:

```
salloc 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h5m50s

input_fasta="clusters.fasta"
output_bed="clusters_nvec_1.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Nvec 2:

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h7m6s

input_fasta="clusters.fasta"
output_bed="clusters_nvec_2.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Nvec 3:

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h8m22s

input_fasta="clusters.fasta"
output_bed="clusters_nvec_3.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Nvec 4:

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h9m37s

input_fasta="clusters.fasta"
output_bed="clusters_nvec_4.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Move all cluster bed files to same folder 

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna
mkdir proTRAC_bed
cd proTRAC_bed

mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h5m50s/clusters_nvec_1.bed .

mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h7m6s/clusters_nvec_2.bed .

mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h8m22s/clusters_nvec_3.bed .

mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m26d19h9m37s/clusters_nvec_4.bed .
```

Merge clusters across samples using bedtools. Following Javi's work, I selected clusters that are present in at least two reps and overlap is at least 50% of the cluster and is reciprocal (ie A overlaps 50% of B and B overlaps 50% of A)

```
module load uri/main
module load bedtools2/2.31.1

FILES=(*.bed)
NUM_FILES=${#FILES[@]}
# Loop over each file
for (( i=0; i<$NUM_FILES; i++ )); do
  for (( j=i+1; j<$NUM_FILES; j++ )); do
    FILE1=${FILES[$i]}
    FILE2=${FILES[$j]}
    
    intersectBed -a $FILE1 -b $FILE2 -f 0.5 -r > ${i}_${j}_merged.bed
  done
done

cat *_merged.bed | sort -k1,1 -k2,2n | bedtools merge -i - > nvec.merged.clusters.bed
```

Here are the merged clusters: 

```
chr11   8499045 8510667
chr12   1380644 1385812
chr12   6786479 6789047
chr2    19149472        19151000
chr4    12785445        12787930
```

Run protrac for Ahya. `nano protrac_ahya.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting protrac for ahya"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna

for f in ahya_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/proTRAC_2.4.2.pl \
-map ${f}/out/other.fq.map.weighted-10000-1000-b-0 \
-genome /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta \
-repeatmasker /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta.out \
-geneset /work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthuns.genes.gff
done

echo "Ahya protrac complete!"
```

Submitted batch job 48140845. The Ahya V1 genome has a coding GFF (`Ahyacinthus.coding.gff3`) and a gene GFF (`Ahyacinthuns.genes.gff`). I'm going to use the gene GFF. Ran successfully. Looking at slurm output file: 

```
## Ahya 1 --
Predicted Cluster: chr1 23805002-23812976 directionality: mono:plus
Predicted Cluster: chr2 174-5987 directionality: mono:minus
Predicted Cluster: chr2 11387569-11391448 directionality: mono:plus
Predicted Cluster: chr2 19927537-19933446 directionality: mono:plus
Predicted Cluster: chr2 20721045-20728689 directionality: mono:plus
Predicted Cluster: chr3 5894666-5898237 directionality: mono:plus
Predicted Cluster: chr3 12591079-12593814 directionality: mono:plus
Predicted Cluster: chr3 25654405-25659422 directionality: mono:minus
Predicted Cluster: chr4 10052163-10060993 directionality: mono:minus
Predicted Cluster: chr4 10805199-10811026 directionality: bi:plus-minus
Predicted Cluster: chr5 6139145-6146999 directionality: mono:minus
Predicted Cluster: chr5 6228164-6249992 directionality: mono:minus
Predicted Cluster: chr5 9867157-9875682 directionality: bi:plus-minus
Predicted Cluster: chr5 10120447-10128784 directionality: mono:plus
Predicted Cluster: chr5 15693366-15701818 directionality: bi:plus-minus
Predicted Cluster: chr6 10853002-10861833 directionality: mono:minus
Predicted Cluster: chr7 8617006-8622163 directionality: mono:minus
Predicted Cluster: chr7 26303828-26305916 directionality: bi:minus-plus
Predicted Cluster: chr8 8099288-8100497 directionality: mono:plus
Predicted Cluster: chr8 9163997-9165256 directionality: mono:minus
Predicted Cluster: chr8 16521447-16526985 directionality: mono:minus
Predicted Cluster: chr8 16543084-16548383 directionality: mono:minus
Predicted Cluster: chr9 10635273-10641316 directionality: bi:plus-minus
Predicted Cluster: chr9 20767689-20774996 directionality: bi:plus-minus
Predicted Cluster: chr9 21206153-21214762 directionality: mono:minus
Predicted Cluster: chr10 1489825-1496287 directionality: bi:plus-minus
Predicted Cluster: chr10 3476068-3485921 directionality: mono:minus
Predicted Cluster: chr10 22459324-22461596 directionality: mono:plus
Predicted Cluster: chr10 25014073-25022388 directionality: mono:plus
Predicted Cluster: chr11 5700676-5705543 directionality: mono:plus
Predicted Cluster: chr11 6497888-6502159 directionality: mono:minus
Predicted Cluster: chr11 9889042-9894279 directionality: mono:plus
Predicted Cluster: chr11 12532048-12539891 directionality: mono:minus
Predicted Cluster: chr11 23481923-23486590 directionality: mono:plus
Predicted Cluster: chr12 11272332-11273942 directionality: mono:minus
Predicted Cluster: chr12 16229398-16234946 directionality: mono:minus
Predicted Cluster: chr13 3333348-3337919 directionality: mono:minus
Predicted Cluster: chr13 5591097-5597882 directionality: mono:plus
Predicted Cluster: chr13 23920835-23927662 directionality: mono:minus
Predicted Cluster: chr14 9854228-9861010 directionality: mono:minus
Predicted Cluster: chr14 10950061-10955002 directionality: mono:minus
Predicted Cluster: chr14 11016252-11023023 directionality: mono:minus
Predicted Cluster: chr14 12238193-12247000 directionality: mono:minus
Predicted Cluster: chr14 12512046-12520014 directionality: mono:plus
Predicted Cluster: sc009 116014-118978 directionality: mono:minus
Predicted Cluster: sc026 288140-299881 directionality: mono:minus
Predicted Cluster: sc032 3019-12954 directionality: mono:minus
Predicted Cluster: sc214 39028-45919 directionality: mono:minus
Predicted Cluster: sc314 19642-27493 directionality: mono:plus

## Ahya 2 --
Predicted Cluster: chr2 173-7648 directionality: mono:minus
Predicted Cluster: chr2 19928300-19930266 directionality: mono:plus
Predicted Cluster: chr2 21342133-21346549 directionality: mono:minus
Predicted Cluster: chr3 5894847-5898618 directionality: mono:plus
Predicted Cluster: chr3 15968285-15971277 directionality: mono:plus
Predicted Cluster: chr3 25655012-25659400 directionality: mono:minus
Predicted Cluster: chr4 10055285-10060652 directionality: mono:minus
Predicted Cluster: chr4 10807669-10812134 directionality: mono:minus
Predicted Cluster: chr4 12717885-12721889 directionality: mono:plus
Predicted Cluster: chr5 6068277-6074866 directionality: mono:minus
Predicted Cluster: chr5 6082459-6089699 directionality: mono:minus
Predicted Cluster: chr5 6228860-6241193 directionality: mono:minus
Predicted Cluster: chr5 10121158-10125619 directionality: mono:plus
Predicted Cluster: chr5 32495982-32499420 directionality: mono:plus
Predicted Cluster: chr8 9163998-9165217 directionality: mono:minus
Predicted Cluster: chr8 16521490-16527007 directionality: mono:minus
Predicted Cluster: chr8 16543317-16548005 directionality: mono:minus
Predicted Cluster: chr9 21210530-21214614 directionality: mono:minus
Predicted Cluster: chr9 21238079-21244556 directionality: mono:minus
Predicted Cluster: chr10 2235047-2239488 directionality: mono:plus
Predicted Cluster: chr10 3478663-3485872 directionality: mono:minus
Predicted Cluster: chr10 19918786-19921758 directionality: mono:minus
Predicted Cluster: chr11 24240166-24252506 directionality: mono:plus
Predicted Cluster: chr14 9855710-9859900 directionality: mono:minus
Predicted Cluster: chr14 13519167-13524011 directionality: mono:plus
Predicted Cluster: sc026 290051-298252 directionality: mono:minus
Predicted Cluster: sc032 4320-10008 directionality: bi:minus-plus
Predicted Cluster: sc103 970-6936 directionality: mono:minus
Predicted Cluster: sc314 19630-23740 directionality: mono:plus

## Ahya 3 --
Predicted Cluster: chr1 2906547-2909858 directionality: mono:minus
Predicted Cluster: chr1 20388557-20397007 directionality: mono:minus
Predicted Cluster: chr1 21490363-21498207 directionality: bi:plus-minus
Predicted Cluster: chr2 19928234-19933515 directionality: mono:plus
Predicted Cluster: chr3 5892792-5898423 directionality: mono:plus
Predicted Cluster: chr3 10973176-10977971 directionality: bi:plus-minus
Predicted Cluster: chr3 15968285-15972928 directionality: mono:plus
Predicted Cluster: chr4 7489511-7494027 directionality: mono:plus
Predicted Cluster: chr4 10052358-10060993 directionality: mono:minus
Predicted Cluster: chr5 6228164-6240022 directionality: mono:minus
Predicted Cluster: chr5 9861144-9875921 directionality: mono:minus
Predicted Cluster: chr5 10118137-10129017 directionality: mono:plus
Predicted Cluster: chr5 32453010-32460018 directionality: mono:plus
Predicted Cluster: chr6 9234303-9242470 directionality: mono:minus
Predicted Cluster: chr6 17937341-17946011 directionality: bi:minus-plus
Predicted Cluster: chr6 21522002-21530816 directionality: mono:plus
Predicted Cluster: chr6 25092025-25094664 directionality: mono:plus
Predicted Cluster: chr7 1392023-1398495 directionality: mono:plus
Predicted Cluster: chr8 9147323-9151516 directionality: bi:plus-minus
Predicted Cluster: chr8 9162754-9165256 directionality: mono:minus
Predicted Cluster: chr8 16521488-16527708 directionality: mono:minus
Predicted Cluster: chr8 16543153-16548377 directionality: mono:minus
Predicted Cluster: chr8 16775067-16780016 directionality: mono:plus
Predicted Cluster: chr8 27351514-27354952 directionality: mono:plus
Predicted Cluster: chr9 20766388-20770995 directionality: mono:plus
Predicted Cluster: chr9 21206083-21214818 directionality: mono:minus
Predicted Cluster: chr10 1022091-1029968 directionality: mono:plus
Predicted Cluster: chr10 3477038-3485921 directionality: mono:minus
Predicted Cluster: chr10 24614082-24621766 directionality: mono:plus
Predicted Cluster: chr11 11939045-11943897 directionality: mono:plus
Predicted Cluster: chr12 1926120-1934885 directionality: mono:plus
Predicted Cluster: chr12 2284111-2292994 directionality: mono:plus
Predicted Cluster: chr12 11272333-11273942 directionality: mono:minus
Predicted Cluster: chr13 13959274-13962228 directionality: mono:minus
Predicted Cluster: chr14 11014211-11023023 directionality: mono:minus
Predicted Cluster: sc007 202804-206412 directionality: mono:minus
Predicted Cluster: sc022 105244-106409 directionality: mono:minus
Predicted Cluster: sc032 1217-12953 directionality: mono:minus
Predicted Cluster: sc115 37093-41903 directionality: mono:minus
Predicted Cluster: sc214 37064-45919 directionality: bi:plus-minus
Predicted Cluster: sc469 12041-18833 directionality: mono:plus

## Ahya 4 --
Predicted Cluster: chr1 21492132-21498143 directionality: bi:minus-plus
Predicted Cluster: chr2 108-5943 directionality: mono:minus
Predicted Cluster: chr2 21043016-21052016 directionality: mono:minus
Predicted Cluster: chr3 5894624-5897787 directionality: mono:plus
Predicted Cluster: chr3 11539554-11547283 directionality: mono:minus
Predicted Cluster: chr3 15968285-15972928 directionality: mono:plus
Predicted Cluster: chr4 9312314-9316951 directionality: mono:minus
Predicted Cluster: chr5 6068139-6076904 directionality: mono:minus
Predicted Cluster: chr5 6241127-6249985 directionality: mono:minus
Predicted Cluster: chr5 9868243-9876990 directionality: mono:minus
Predicted Cluster: chr6 10805163-10813996 directionality: mono:minus
Predicted Cluster: chr8 9146963-9150640 directionality: bi:plus-minus
Predicted Cluster: chr8 16521490-16527634 directionality: mono:minus
Predicted Cluster: chr8 16626153-16632031 directionality: mono:minus
Predicted Cluster: chr8 16775065-16780646 directionality: mono:plus
Predicted Cluster: chr9 20769502-20773911 directionality: mono:plus
Predicted Cluster: chr9 21206097-21212016 directionality: mono:minus
Predicted Cluster: chr9 21238107-21246779 directionality: mono:minus
Predicted Cluster: chr10 3477110-3485774 directionality: mono:minus
Predicted Cluster: chr11 12582678-12588013 directionality: mono:minus
Predicted Cluster: chr11 24245099-24251987 directionality: mono:plus
Predicted Cluster: chr13 5667030-5673887 directionality: mono:plus
Predicted Cluster: chr14 12567517-12574222 directionality: mono:plus
Predicted Cluster: chr14 14424200-14432903 directionality: bi:minus-plus
Predicted Cluster: sc026 293015-297872 directionality: mono:minus
Predicted Cluster: sc032 4061-10008 directionality: mono:minus
Predicted Cluster: sc178 5793-10354 directionality: mono:plus
```

Run ping pong ID on all reads. `nano pingpong_ahya.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0

echo "Starting ping pong for ahya"

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna

for f in ahya_*_L001_R1_001_trim.fastq.collapsed.filt.no-dust
do
perl /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/ngs_toolbox/TBr2_pingpong.pl \
-i ${f}/out/other.fq.map.weighted-10000-1000-b-0 \
-o ${f}.pp
done

echo "Ahya ping pong complete!"
```

Submitted batch job 48141938

Look at cluster overlap between the 4 reps. Create bed files from cluster fasta files.

Ahya 1:

```
salloc 

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m27d14h18m54s

input_fasta="clusters.fasta"
output_bed="clusters_ahya_1.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Ahya 2:

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m27d14h20m26s

input_fasta="clusters.fasta"
output_bed="clusters_ahya_2.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Ahya 3:

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m27d14h21m36s

input_fasta="clusters.fasta"
output_bed="clusters_ahya_3.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Ahya 4:

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m27d14h23m12s

input_fasta="clusters.fasta"
output_bed="clusters_ahya_4.bed"

# Extract relevant info from header and convert to BED
# BED format: chrom  start(0-based)  end(1-based)  .  0  strand(+/-)

grep "^>" "$input_fasta" | while read -r header; do
  # Remove leading '>'
  header=${header#>}
  
  # Extract chromosome (first word)
  chrom=$(echo "$header" | awk '{print $1}')
  
  # Extract coordinates (second word), split by '-'
  coords=$(echo "$header" | awk '{print $2}')
  start=$(echo "$coords" | cut -d'-' -f1)
  end=$(echo "$coords" | cut -d'-' -f2)
  
  # Convert to 0-based start for BED
  start=$((start - 1))
  
  # Extract strand from "directionality:" field (+ or -)
  # directionality: mono:plus means '+'
  # directionality: mono:minus means '-'
  strand=$(echo "$header" | grep -oP "directionality: mono:\K(plus|minus)")
  if [ "$strand" = "plus" ]; then
    strand="+"
  else
    strand="-"
  fi
  
  # Print BED line
  echo -e "${chrom}\t${start}\t${end}\t.\t0\t${strand}"
done > "$output_bed"
```

Move all cluster bed files to same folder 

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna
mkdir proTRAC_bed
cd proTRAC_bed

mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m27d14h18m54s/clusters_ahya_1.bed .

mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m27d14h20m26s/clusters_ahya_2.bed .

mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m27d14h21m36s/clusters_ahya_3.bed .

mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_other.fq.map.weighted-10000-1000-b-0_2025y10m27d14h23m12s/clusters_ahya_4.bed .
```

Merge clusters across samples using bedtools. Following Javi's work, I selected clusters that are present in at least two reps and overlap is at least 50% of the cluster and is reciprocal (ie A overlaps 50% of B and B overlaps 50% of A)

```
module load uri/main
module load bedtools2/2.31.1

FILES=(*.bed)
NUM_FILES=${#FILES[@]}
# Loop over each file
for (( i=0; i<$NUM_FILES; i++ )); do
  for (( j=i+1; j<$NUM_FILES; j++ )); do
    FILE1=${FILES[$i]}
    FILE2=${FILES[$j]}
    
    intersectBed -a $FILE1 -b $FILE2 -f 0.5 -r > ${i}_${j}_merged.bed
  done
done

cat *_merged.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ahya.merged.clusters.bed
```

Here are the merged clusters: 

```

chr1    21492131        21498143
chr10   3477037 3485921
chr11   24245098        24251987
chr12   11272332        11273942
chr14   9855709 9859900
chr14   11016251        11023023
chr2    172     5987
chr2    19928233        19933446
chr3    5894623 5898423
chr3    15968284        15972928
chr3    25655011        25659400
chr4    10052357        10060993
chr4    10807668        10811026
chr5    6068276 6074866
chr5    6228163 6241193
chr5    9867156 9875921
chr5    10120446        10128784
chr8    9147322 9150640
chr8    9163996 9165256
chr8    16521487        16527634
chr8    16543152        16548377
chr8    16775066        16780016
chr9    20769501        20773911
chr9    21206096        21214762
chr9    21238106        21244556
sc026   290050  298252
sc032   3018    12953
sc214   39027   45919
sc314   19641   23740
```

Trying to determine exactly how many piRNA reads were identified. 

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_1_S27_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

# Input file: other.fq.map.weighted-10000-1000-b-0 (example name)
input="other.fq.map.weighted-10000-1000-b-0"

# Extract sequences (field 3)
cut -f3 $input > seqs.txt

# Compute counts of 1U (1st base U) and 10A (10th base A) relative to total reads
total_reads=$(wc -l < seqs.txt)

# Count 1U reads
count_1U=$(awk 'substr($0,1,1)=="T" || substr($0,1,1)=="U" {count++} END {print count+0}' seqs.txt)

# Count 10A reads (check length >= 10 to avoid short seqs)
count_10A=$(awk 'length($0)>=10 && substr($0,10,1)=="A" {count++} END {print count+0}' seqs.txt)

echo "Total reads (25-35 nt): $total_reads"
echo "Reads with 1U at position 1: $count_1U ($(awk "BEGIN {printf \"%.2f\", ($count_1U/$total_reads)*100}"))%"
echo "Reads with A at position 10: $count_10A ($(awk "BEGIN {printf \"%.2f\", ($count_10A/$total_reads)*100}"))%"

Total reads (25-35 nt): 1,388,856
Reads with 1U at position 1: 1130221 (81.38)%
Reads with A at position 10: 642051 (46.23)%
```

I don't know...I think I just need to make a decision and go from there. 

Look at the intersection of piRNA clusters and TE. In R, I did the following: 

Read in repeatmasker out file 
```{r}
out_file <- read.table("~/Desktop/GFFs/ahya/Ahyacinthus.chrsV1.fasta.out", comment.char = "", header = F, fill = T, skip = 2)
```

Retain only chromosome, start, stop, TE name, Family Name and percent divergence
```{r}
bed_file <- out_file[,c(5:7,11,2)]

#remove empty lines 
bed_file <- bed_file %>%
  drop_na(.)
```

Write out bed file
```{r}
write.table(bed_file, file = "~/Desktop/GFFs/ahya/Ahya_repeatmasker.bed", sep = "\t", col.names = F, row.names = F, quote = F)
```

Upload bed file to `/scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya`. Then intersect the TE file with the piRNA cluster bed file. 

```
module load bedtools2/2.31.1

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_bed

bedtools intersect -wo -a ahya.merged.clusters.bed -b /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahya_repeatmasker.bed > Ahya_piRNA_cluster_TEs_intersect.txt
```

I need to make bed files for each piRNA reads per rep and then merge them and then intersect them with the other piRNA bed files AND repeat masker. Make bed file for piRNAs first. 

```
## Ahya 1
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_1_S27_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

awk '{
  chr = $1;
  start = $2 - 1;                  # convert 1-based to 0-based start for BED
  end = start + length($3);       # end is start + length of sequence in field 3
  print chr"\t"start"\t"end;
}' other.fq.map.weighted-10000-1000-b-0 > other.fq.map.weighted-10000-1000-b-0.ahya_1.bed

## Ahya 2
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_2_S28_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

awk '{
  chr = $1;
  start = $2 - 1;                  # convert 1-based to 0-based start for BED
  end = start + length($3);       # end is start + length of sequence in field 3
  print chr"\t"start"\t"end;
}' other.fq.map.weighted-10000-1000-b-0 > other.fq.map.weighted-10000-1000-b-0.ahya_2.bed

## Ahya 3
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_3_S29_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

awk '{
  chr = $1;
  start = $2 - 1;                  # convert 1-based to 0-based start for BED
  end = start + length($3);       # end is start + length of sequence in field 3
  print chr"\t"start"\t"end;
}' other.fq.map.weighted-10000-1000-b-0 > other.fq.map.weighted-10000-1000-b-0.ahya_3.bed

## Ahya 4
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_4_S30_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

awk '{
  chr = $1;
  start = $2 - 1;                  # convert 1-based to 0-based start for BED
  end = start + length($3);       # end is start + length of sequence in field 3
  print chr"\t"start"\t"end;
}' other.fq.map.weighted-10000-1000-b-0 > other.fq.map.weighted-10000-1000-b-0.ahya_4.bed
```

Move new bed files into their own folder and assess overlap across reps

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna
mkdir piRNA_bed
cd piRNA_bed
mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_1_S27_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out/other.fq.map.weighted-10000-1000-b-0.ahya_1.bed .
mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_2_S28_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out/other.fq.map.weighted-10000-1000-b-0.ahya_2.bed .
mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_3_S29_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out/other.fq.map.weighted-10000-1000-b-0.ahya_3.bed .
mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/ahya_4_S30_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out/other.fq.map.weighted-10000-1000-b-0.ahya_4.bed .

module load bedtools2/2.31.1

FILES=(*.bed)
NUM_FILES=${#FILES[@]}
# Loop over each file
for (( i=0; i<$NUM_FILES; i++ )); do
  for (( j=i+1; j<$NUM_FILES; j++ )); do
    FILE1=${FILES[$i]}
    FILE2=${FILES[$j]}
    
    intersectBed -a $FILE1 -b $FILE2 -f 0.5 -r > ${i}_${j}_merged.bed
  done
done

cat *_merged.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ahya.merged.piRNA.bed
wc -l ahya.merged.piRNA.bed 
364493 ahya.merged.piRNA.bed
```

Make fasta for piRNAs 
```
bedtools getfasta -fi /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.chrsV1.fasta -bed ahya.merged.piRNA.bed -fo ahya.merged.piRNA.fasta
```

See % of seqs that start with U/T and that have A at position 10. 

```
awk '/^>/ {getline seq; if(substr(seq,1,1)=="U" || substr(seq,1,1)=="T") ustart++; if(substr(seq,10,1)=="A") apos10++; total++} END {print "Total seqs:", total; print "Start with U/T:", ustart; print "A at pos 10:", apos10; print "% U/T start:", 100*ustart/total; print "% A pos 10:", 100*apos10/total}' ahya.merged.piRNA.fasta
Total seqs: 364493
Start with U/T: 182922
A at pos 10: 130219
% U/T start: 50.1853
% A pos 10: 35.7261
```

Intersect piRNA bed file with the repeatmasker bed file 

```
bedtools intersect -wo -a ahya.merged.piRNA.bed -b /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahya_repeatmasker.bed > Ahya_piRNAs_TEs_intersect.txt
```

Intersect piRNA bed file and piRNA cluster bed file with the gene and coding gffs. 

```
## piRNAs 
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/piRNA_bed

bedtools intersect -wo -a ahya.merged.piRNA.bed -b /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthuns.genes.gff > Ahya_piRNAs_genes_intersect.txt

bedtools intersect -wo -a ahya.merged.piRNA.bed -b /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.coding.gff3 > Ahya_piRNAs_coding_intersect.txt

## Clusters 
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/sortmerna/proTRAC_bed

bedtools intersect -wo -a ahya.merged.clusters.bed -b /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthuns.genes.gff > Ahya_piRNA_cluster_genes_intersect.txt

bedtools intersect -wo -a ahya.merged.clusters.bed -b /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya/Ahyacinthus.coding.gff3 > Ahya_piRNA_cluster_coding_intersect.txt
```

Do the same for Apoc. Look at the intersection of piRNA clusters and TE. In R, I did the following:

Read in repeatmasker out file. 

```
out_file <- read.table("~/Desktop/GFFs/apoc/apoculata.genome.fasta.out", comment.char = "", header = F, fill = T, skip = 2)
```

Retain only chromosome, start, stop, TE name, Family Name and percent divergence

```
bed_file <- out_file[,c(5:7,11,2)]

#remove empty lines 
bed_file <- bed_file %>%
  drop_na(.)
```

Write out bed file

```
write.table(bed_file, file = "~/Desktop/GFFs/apoc/Apoc_repeatmasker.bed", sep = "\t", col.names = F, row.names = F, quote = F)
```

Upload bed file to `/scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/ahya`. Then intersect the TE file with the piRNA cluster bed file.

```
module load bedtools2/2.31.1

cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/proTRAC_bed

bedtools intersect -wo -a apoc.merged.clusters.bed -b /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/Apoc_repeatmasker.bed > Apoc_piRNA_cluster_TEs_intersect.txt
```

I need to make bed files for each piRNA reads per rep and then merge them and then intersect them with the other piRNA bed files AND repeat masker. Make bed file for piRNAs first.

```
## Apoc 2
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/apoc_2_S31_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

awk '{
  chr = $1;
  start = $2 - 1;                  # convert 1-based to 0-based start for BED
  end = start + length($3);       # end is start + length of sequence in field 3
  print chr"\t"start"\t"end;
}' other.fq.map.weighted-10000-1000-b-0 > other.fq.map.weighted-10000-1000-b-0.apoc_2.bed

## Apoc 3
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/apoc_3_S32_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

awk '{
  chr = $1;
  start = $2 - 1;                  # convert 1-based to 0-based start for BED
  end = start + length($3);       # end is start + length of sequence in field 3
  print chr"\t"start"\t"end;
}' other.fq.map.weighted-10000-1000-b-0 > other.fq.map.weighted-10000-1000-b-0.apoc_3.bed

## Apoc 4
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/apoc_4_S33_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out

awk '{
  chr = $1;
  start = $2 - 1;                  # convert 1-based to 0-based start for BED
  end = start + length($3);       # end is start + length of sequence in field 3
  print chr"\t"start"\t"end;
}' other.fq.map.weighted-10000-1000-b-0 > other.fq.map.weighted-10000-1000-b-0.apoc_4.bed
```

Move new bed files into their own folder and assess overlap across reps

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna
mkdir piRNA_bed
cd piRNA_bed
mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/apoc_2_S31_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out/other.fq.map.weighted-10000-1000-b-0.apoc_2.bed .
mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/apoc_3_S32_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out/other.fq.map.weighted-10000-1000-b-0.apoc_3.bed .
mv /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/apoc_4_S33_L001_R1_001_trim.fastq.collapsed.filt.no-dust/out/other.fq.map.weighted-10000-1000-b-0.apoc_4.bed .

FILES=(*.bed)
NUM_FILES=${#FILES[@]}
# Loop over each file
for (( i=0; i<$NUM_FILES; i++ )); do
  for (( j=i+1; j<$NUM_FILES; j++ )); do
    FILE1=${FILES[$i]}
    FILE2=${FILES[$j]}
    
    intersectBed -a $FILE1 -b $FILE2 -f 0.5 -r > ${i}_${j}_merged.bed
  done
done

cat *_merged.bed | sort -k1,1 -k2,2n | bedtools merge -i - > apoc.merged.piRNA.bed
wc -l apoc.merged.piRNA.bed 
5537 apoc.merged.piRNA.bed
```

Much less piRNAs than the Ahya...

Make fasta for piRNAs 
```
bedtools getfasta -fi /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta -bed apoc.merged.piRNA.bed -fo apoc.merged.piRNA.fasta
```

See % of seqs that start with U/T and that have A at position 10. 

```
awk '/^>/ {getline seq; if(substr(seq,1,1)=="U" || substr(seq,1,1)=="T") ustart++; if(substr(seq,10,1)=="A") apos10++; total++} END {print "Total seqs:", total; print "Start with U/T:", ustart; print "A at pos 10:", apos10; print "% U/T start:", 100*ustart/total; print "% A pos 10:", 100*apos10/total}' apoc.merged.piRNA.fasta
Total seqs: 5537
Start with U/T: 2233
A at pos 10: 1742
% U/T start: 40.3287
% A pos 10: 31.4611
```

Intersect piRNA bed file with the repeatmasker bed file

```
bedtools intersect -wo -a apoc.merged.piRNA.bed -b /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/Apoc_repeatmasker.bed > Apoc_piRNAs_TEs_intersect.txt
```

Intersect piRNA bed file and piRNA cluster bed file with the gene gff.

```
## piRNAs
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/piRNA_bed

bedtools intersect -wo -a apoc.merged.piRNA.bed -b /work/pi_hputnam_uri_edu/genomes/Apoc/apoc_GFFannotation.gene.gff > Apoc_piRNAs_genes_intersect.txt

## piRNA clusters
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/sortmerna/proTRAC_bed

bedtools intersect -wo -a apoc.merged.clusters.bed -b /work/pi_hputnam_uri_edu/genomes/Apoc/apoc_GFFannotation.gene.gff > Apoc_piRNA_cluster_genes_intersect.txt
```



























### tRNA derived RNAs 

Nick and I have been talking about other smRNAs to identify. He studies tRNA derived RNAs (tDRs), which are small RNAs that originate from tRNAs. Based on evidence from [model systems](https://pmc.ncbi.nlm.nih.gov/articles/PMC10766869/), these RNAs can perform selective translation regulation (miRNA-like mechanism and interaction with RNA binding proteins) or global translation regulation (eIF4F sequestration and interaction with ribosomes). However, even in model systems, work is pretty limited and often is more experimental than bioinformatical. Nick has a [paper](https://www.biorxiv.org/content/10.1101/2025.04.14.648817v1#:~:text=Abstract,paternal%20non%2Dgenetic%20inheritance%20mechanistically.) on biorxiv that exakined tDRs as a non-genetic inheritance mechanism via sperm in Celegans, which is super cool. 

The tools that I have found for tDR bioinformatic analysis are [tsRNAsearch](https://academic.oup.com/bioinformatics/article/37/23/4424/6320783) (github [here](https://github.com/GiantSpaceRobot/tsRNAsearch); published 2021) and [tDRmapper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0800-0) (github [here](https://github.com/sararselitsky/tDRmapper); published 2015). In the tsRNAsearch paper, it compared the two tools and found "The comparison between tsRNAsearch and tDRmapper using the heavily trimmed data resulted in an average Pearson’s r = 0.70 and a standard deviation of ±0.15." Pretty comparable. 

I want to try out tDRmapper first. From the [github](https://github.com/sararselitsky/tDRmapper/tree/master), it looks like there is this primary script `Scripts/TdrMappingScripts.pl`, which references other scripts in the scripts folder. That is kinda confusing. I uploaded the perl and R scripts here: `/work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts/tDRmapper_scripts`. Let's try to run! `nano tDRmapper_apoc_2.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Perl/5.40.0-GCCcore-14.2.0
module load R/4.3.2-gfbf-2023a

echo "Starting apoc 2 putative tRNA mapping"

cd /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

perl tDRmapper_scripts/TdrMappingScripts.pl /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/Apoc-tRNA.fasta /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc/apoc_2_S31_L001_R1_001_trim.fastq

echo "Apoc 2 putative tRNA mapping complete!"
```

Submitted batch job 47018721. Failed immediately with the error "no fasta detected"...Changed so that I `cd`ed right into the apoc directory. Submitted batch job 47021871. Still getting the same error...Here are some solutions from perplexity: 

- If the file uses CRLF (\r\n), Perl may interpret > as part of a different string. Fix by: `dos2unix Apoc-tRNA.fasta` - did this. 
- Make sure its uncompressed plain text with `file Apoc-tRNA.fasta` - did this, printed `Apoc-tRNA.fasta: ASCII text`. 

Let's see if it works. Submitted batch job 47023410. Ran for a little longer but still ended with same error. 

Maybe just pivot to tsRNAsearch...They provide installation code on their github. I'm going to run it as a job since it may take a while. `nano install_tsRNAsearch.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 50:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

module load uri/main
module load Miniconda3/4.9.2
module load R/4.3.2-gfbf-2023a
module load Nextflow/22.04.0

echo "Installing tsRNAsearch via conda"

cd /work/pi_hputnam_uri_edu/conda/envs
git clone https://github.com/GiantSpaceRobot/tsRNAsearch.git
conda env create -f tsRNAsearch/environment.yml
conda activate tsrnasearch_env

echo "Installing R packages for tsRNAsearch"

Rscript tsRNAsearch/bin/InstallPackages.R

echo "Test run tsRNAsearch"
nextflow run tsRNAsearch --species mouse --input_dir tsRNAsearch/ExampleData --output_dir Results
```

Submitted batch job 47067938. Okay this ran, but the environment still had conflicts and failed...maybe message Unity help desk about this. 


tsRNAsearch is a nextflow pipeline, and nextflow is already installed on the HPC. 


http://rna.sysu.edu.cn/tsRFun/analysisResults/JC8G8OKDOQ_FUNC/index_FUNC.php

https://rna.sysu.edu.cn/tsRFun/analysisResults/0/index.php


I also found [tRFTarget](http://trftarget.net/online_targets), which can be used to identify the targets of tDRs. They have an online tool which seems very nifty. 

Maybe try [tsRFinder](https://github.com/zhlingl/tsRFun/tree/main/tsRFinder)

or mintMAP

How about I just try to map the reads to the tRNA fasta??








### Other stuff 

Extract chromosome names from fasta files 

```
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc
sed -n '/^>/p' /work/pi_hputnam_uri_edu/genomes/Apoc/apoculata.genome.fasta > apoc_chroms.txt

cd ../nvec/genome
sed -n '/^>/p' Nvec200.fasta > nvec_chroms.txt

cd ../../ahya
sed -n '/^>/p' Ahyacinthus.chrsV1.fasta > ahya_chroms.txt
```

Apoc has 14 chromosomes, Nvec has 30 chromosomes, and Ahya has 907 chromosomes (more like contigs). 

Extract data for length/nucleotide plotting from each spp. 

```
## Apoc 
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/apoc
awk 'NR%4==2 {len=length($0); 
              A=gsub(/A/,"",$0);
              C=gsub(/C/,"",$0);
              G=gsub(/G/,"",$0);
              U=gsub(/T/,"",$0);
              print len, A, C, G, U}' apoc_2_S31_L001_R1_001_trim.fastq > apoc_2_base_counts.txt
              
awk 'NR%4==2 {len=length($0); 
              A=gsub(/A/,"",$0);
              C=gsub(/C/,"",$0);
              G=gsub(/G/,"",$0);
              U=gsub(/T/,"",$0);
              print len, A, C, G, U}' apoc_3_S32_L001_R1_001_trim.fastq > apoc_3_base_counts.txt

awk 'NR%4==2 {len=length($0); 
              A=gsub(/A/,"",$0);
              C=gsub(/C/,"",$0);
              G=gsub(/G/,"",$0);
              U=gsub(/T/,"",$0);
              print len, A, C, G, U}' apoc_4_S33_L001_R1_001_trim.fastq > apoc_4_base_counts.txt

## Nvec 
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/nvec
awk 'NR%4==2 {len=length($0); 
              A=gsub(/A/,"",$0);
              C=gsub(/C/,"",$0);
              G=gsub(/G/,"",$0);
              U=gsub(/T/,"",$0);
              print len, A, C, G, U}' nvec_1_S25_L001_R1_001_trim.fastq > nvec_1_base_counts.txt

awk 'NR%4==2 {len=length($0); 
              A=gsub(/A/,"",$0);
              C=gsub(/C/,"",$0);
              G=gsub(/G/,"",$0);
              U=gsub(/T/,"",$0);
              print len, A, C, G, U}' nvec_2_S26_L001_R1_001_trim.fastq > nvec_2_base_counts.txt

awk 'NR%4==2 {len=length($0); 
              A=gsub(/A/,"",$0);
              C=gsub(/C/,"",$0);
              G=gsub(/G/,"",$0);
              U=gsub(/T/,"",$0);
              print len, A, C, G, U}' nvec_3_S34_L001_R1_001_trim.fastq > nvec_3_base_counts.txt

awk 'NR%4==2 {len=length($0); 
              A=gsub(/A/,"",$0);
              C=gsub(/C/,"",$0);
              G=gsub(/G/,"",$0);
              U=gsub(/T/,"",$0);
              print len, A, C, G, U}' nvec_4_S35_L001_R1_001_trim.fastq > nvec_4_base_counts.txt
```

5' end 

```
## Apoc 
awk 'NR%4==2 {print length($1), substr($1,1,1)}' apoc_2_S31_L001_R1_001_trim.fastq > apoc_2_first_base_by_length.txt

awk 'NR%4==2 {print length($1), substr($1,1,1)}' apoc_3_S32_L001_R1_001_trim.fastq > apoc_3_first_base_by_length.txt

awk 'NR%4==2 {print length($1), substr($1,1,1)}' apoc_4_S33_L001_R1_001_trim.fastq > apoc_4_first_base_by_length.txt

cat *_first_base_by_length.txt > apoc_combined_first_base.txt

sort apoc_combined_first_base.txt | uniq -c | awk '{print $2, $3, $1}' > apoc_first_base_counts.txt

## Nvec 
awk 'NR%4==2 {print length($1), substr($1,1,1)}' nvec_1_S25_L001_R1_001_trim.fastq > nvec_1_first_base_by_length.txt

awk 'NR%4==2 {print length($1), substr($1,1,1)}' nvec_2_S26_L001_R1_001_trim.fastq > nvec_2_first_base_by_length.txt

awk 'NR%4==2 {print length($1), substr($1,1,1)}' nvec_3_S34_L001_R1_001_trim.fastq > nvec_3_first_base_by_length.txt

awk 'NR%4==2 {print length($1), substr($1,1,1)}' nvec_4_S35_L001_R1_001_trim.fastq > nvec_4_first_base_by_length.txt

cat *_first_base_by_length.txt > nvec_combined_first_base.txt

sort nvec_combined_first_base.txt | uniq -c | awk '{print $2, $3, $1}' > nvec_first_base_counts.txt
```

Lengths

```
## Apoc 
awk 'NR%4==2 {print length($1)}' apoc_2_S31_L001_R1_001_trim.fastq > apoc_2_lengths.txt
awk 'NR%4==2 {print length($1)}' apoc_3_S32_L001_R1_001_trim.fastq > apoc_3_lengths.txt
awk 'NR%4==2 {print length($1)}' apoc_4_S33_L001_R1_001_trim.fastq > apoc_4_lengths.txt
cat *_lengths.txt > apoc_combined_lengths.txt

## Nvec
awk 'NR%4==2 {print length($1)}' nvec_1_S25_L001_R1_001_trim.fastq > nvec_1_lengths.txt
awk 'NR%4==2 {print length($1)}' nvec_2_S26_L001_R1_001_trim.fastq > nvec_2_lengths.txt
awk 'NR%4==2 {print length($1)}' nvec_3_S34_L001_R1_001_trim.fastq > nvec_3_lengths.txt
awk 'NR%4==2 {print length($1)}' nvec_4_S35_L001_R1_001_trim.fastq > nvec_4_lengths.txt
cat *_lengths.txt > nvec_combined_lengths.txt

```

### Orthologous genes 

Use [Broccoli](https://github.com/rderelle/Broccoli) to identify orthologous groups to see if there are any similar proteins being targeted by miRNAs and/or being overlapped by piRNAs. Using Zoe's [code](https://github.com/zdellaert/multi-sp-snRNA/blob/main/reference_genes/Marker-Genes.md) as a reference. I downloaded all protein fasta files here: `/scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/protein_fastas`. Rename all the protein fasta files so they are named `spp_proteins.faa`

In the scripts folder: `nano broccoli_max_likelihood.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/cnidarian_sperm_smRNA/scripts

#load conda and activate conda environment
module load conda/latest
conda activate /work/pi_hputnam_uri_edu/conda/envs/env-broccoli

#load additional programs needed to run broccoli
module load uri/main
module load diamond/2.1.7
module load all/FastTree/2.1.11-GCCcore-12.3.0

mkdir -p /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/broccoli_max_likelihood
cd /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/broccoli_max_likelihood

python /work/pi_hputnam_uri_edu/conda/envs/env-broccoli/Broccoli/broccoli.py -dir /scratch3/workspace/jillashey_uri_edu-cnidarian_sperm/protein_fastas -phylogenies 'ml' -ext '.faa' -path_fasttree FastTree

conda deactivate
```

Submitted batch job 48680375




 