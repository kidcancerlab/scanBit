
variant_cutoff=1

region_ranges=$( \
    tail -n +2 output/merged_c5_variable_variants.tsv \
    | awk '$7 >= 1 {print $1 ":" $2-1 "-" $2}' \
    | paste -s -d " " \
    )

region_list=$( \
    tail -n +2 output/merged_c5_variable_variants.tsv \
    | awk '$7 >= 1 {print $1 ":" $2}' \
    | paste -s -d "," \
    )

bam_file=/home/gdrobertslab/lab/Counts_2/S0058/possorted_genome_bam.bam
cell_barcode=ACCCAAAAGGTGCTGA-1
ref_genome=/home/gdrobertslab/lab/GenRef/10x-human/fasta/genome.fa
ploidy=GRCh38

run_mpileup() {
    local cell_barcode=$1

    samtools view \
            -u \
            --threads 5 \
            --tag CB:${cell_barcode} \
            ${bam_file} \
            ${region_ranges} \
        | samtools view \
            -u \
            --threads 5 \
            --tag xf:25 \
        | samtools sort \
        | bcftools mpileup \
            --threads 5 \
            --max-depth 8000 \
            --regions ${region_list} \
            --annotate FORMAT/DP \
            -O u \
            --fasta-ref ${ref_genome} \
            - \
        | bcftools call \
            --threads 5 \
            --annotate GQ \
            -O u \
            --multiallelic-caller \
            --ploidy ${ploidy} \
        | bcftools filter \
            --threads 5 \
            -e "FORMAT/DP>0" \
            --SnpGap 10 \
            -O u \
        | bcftools view \
            --threads 5 \
            --no-header \
            --exclude-types indels \
            -O v \
        > test2.vcf

}


