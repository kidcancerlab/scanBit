#!/bin/bash
#SBATCH --job-name=count_R0075
#SBATCH --output=count-%j.out
#SBATCH --error=count-%j.out
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=general,himem

ml purge
ml miniforge3
eval "$(conda shell.bash hook)"

conda activate scanBit_xkcd_1337 

# identify SNVs on chr19 and create a bed file of their positions
# This is so we can grab only informative reads and keep the bam size small
# This is based on previous analyses
bcftools view \
    ../../../24_validate_snvs/output/snv/mouse/mergedmouse_snvs_c30.bcf \
    chr19 \
  | grep -v "^#" \
  | cut -f 1,2 \
  | head -n 1000 \
  > chr19_snvs.bed


# Make a tiny bam for B6
samtools view \
    -b \
    /home/gdrobertslab/lab/Counts_2/S0066/possorted_genome_bam.bam \
    chr19 \
  | samtools view \
    -L chr19_snvs.bed \
    -b \
  > B6_chr19_scRNA.bam
samtools index B6_chr19_scRNA.bam

# Make a tiny bam for BALBC
samtools view \
    -b \
    /home/gdrobertslab/lab/Counts_2/S0074/possorted_genome_bam.bam \
    chr19 \
  | samtools view \
    -L chr19_snvs.bed \
    -b \
  > BALBc_chr19_scRNA.bam
samtools index BALBc_chr19_scRNA.bam


# Make tiny reference genome with only chr19
cp /home/gdrobertslab/lab/GenRef/10x-mm10/fasta/genome.fa .
samtools faidx genome.fa
samtools faidx genome.fa chr19 > mm10_chr19.fa
rm genome.fa genome.fa.fai

# This is putting tiny placeholder chromosomes into the fasta file so there is
# not a mismatch between the bam header and the fasta
samtools view \
    -H \
    BALBc_chr19_scRNA.bam \
  | grep "^@SQ" \
  | grep -v chr19 \
  | perl -pe 's/.+SN:/>/' \
  | perl -pe 's/\t.+/\nATGC/' \
  | cat mm10_chr19.fa - \
  > mm10_chr19_plus.fa

# Get rid of a few intermediate files
rm mm10_chr19.fa chr19_snvs.bed

conda deactivate
