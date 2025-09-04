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

MultiQC results are here: XXXXX

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

MultiQC results are here: XXXXX

Data looks good. Ahya 2 dropped quite a bit in number of reads but that sample also had the most adapter content. 

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

Submitted batch job 42222961. Still got the same result...I wonder if I can relax the shortstack parameters? I'm going to look more at the shortstack github [wiki](https://github.com/MikeAxtell/ShortStack/wiki/Vignette-%231-%3A-A-%22complete%22-run). Interestingly, the recommended inputs are "one or more untrimmed FASTQ files of small RNA-seq data". Shortstack has an `--autotrim` flag, which (according to the wiki): "tells ShortStack that the reads from --readfile need to be adapter trimmed, and that ShortStack should infer the adapter sequence. This is the recommended way to trim the adapters off of sRNA-seq data." Let's try running the same code but with the raw reads instead. 

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

Submitted batch job 42223404. Still the same problem...only IDing three miRNAs. It seems to be IDing other known miRNAs but just not marking them as miRNAs? Looking at the `Results.txt` file: 

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

Submitted batch job 42225432. Still not identifying a lot of miRNAs. In the results.txt file, I am seeing though that we are getting hits to the known miRNAs (including the ones from the Apoc adults). 

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

Submitted batch job 42239097. Doesn't look much different bleh. 
















 