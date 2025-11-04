#!/bin/sh
#SBATCH --account=placeholder_account
#SBATCH --output=placeholder_job_log
#SBATCH --error=placeholder_job_log
#SBATCH --job-name=cell_var_calls
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
###SBATCH --mem=2G Figure this out
#SBATCH --wait
placeholder_job_header_other

placeholder_batch_other

set -e ### stops bash script if line ends with error

conda activate scanBit_xkcd_1337

start_time=$(date +%s)

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID} Beginning: $(date '+%Y-%m-%d %H:%M:%S')

variant_cutoff=placeholder_variant_cutoff

# genomic region ranges to pass to samtools view
# *** Make sure the 8th column is cross_minus_within_dist ***
# I should determine this programatically
export region_ranges=$( \
    tail -n +2 placeholder_variable_variants \
    | awk -v vc=$variant_cutoff '$8 >= vc {print $1 ":" $2-1 "-" $2}' \
    | paste -s -d " " \
    )

# genomic regions for variants of interest to
export region_list=$( \
    tail -n +2 placeholder_variable_variants \
    | awk -v vc=$variant_cutoff '$8 >= vc {print $1 ":" $2}' \
    | paste -s -d "," \
    )

export temp_dir=placeholder_temp_dir

cell_barcode_array=(placeholder_cell_barcodes)

run_mpileup() {
    local cell_barcode=$1
    # We'll likely have some cells with no reads at any of the loci, so this
    # will keep the program from choking on that
    set +e

    samtools view \
            -u \
            --threads 2 \
            --tag CB:${cell_barcode} \
            placeholder_bam_file \
            ${region_ranges} \
        | samtools view \
            -u \
            --threads 2 \
            --tag xf:25 \
        | samtools sort \
        > ${temp_dir}/sorted_${cell_barcode}.bam

    samtools index ${temp_dir}/sorted_${cell_barcode}.bam

    bcftools mpileup \
            --threads 2 \
            --max-depth 8000 \
            --regions ${region_list} \
            --annotate FORMAT/DP \
            -O u \
            --fasta-ref placeholder_ref_fasta \
            ${temp_dir}/sorted_${cell_barcode}.bam \
        | bcftools call \
            --threads 2 \
            --annotate GQ \
            -O u \
            --multiallelic-caller \
            --ploidy placeholder_ploidy \
        | bcftools filter \
            --threads 2 \
            -e "FORMAT/DP<=placeholder_min_depth" \
            --SnpGap 10 \
            -O u \
        | bcftools view \
            --threads 2 \
            --no-header \
            --exclude-types indels \
            -O v \
        | awk \
            -v cb=$cell_barcode \
            'BEGIN {OFS = "\t"}
            {print cb, $1, $2, $4, $5, $10}' \
        | perl -pe 's/:.+//' \
        > ${temp_dir}/${cell_barcode}.tsv

    rm \
        ${temp_dir}/sorted_${cell_barcode}.bam \
        ${temp_dir}/sorted_${cell_barcode}.bam.bai
}

export -f run_mpileup

parallel -j 20 run_mpileup ::: ${cell_barcode_array[@]}

all_files=(${cell_barcode_array[@]})

for i in ${!all_files[@]}
do
    all_files[$i]=${temp_dir}/${all_files[$i]}.tsv
done

echo -e "cell_barcode\tchr\tpos\tref\talt\tcall" > placeholder_out_file
cat ${all_files[@]} >> placeholder_out_file

rm ${all_files[@]}

end_time=$(date +%s)

elapsed_seconds=$((end_time - start_time))

echo Done: $(date '+%Y-%m-%d %H:%M:%S')
echo Elapsed seconds: $elapsed_seconds

conda deactivate
