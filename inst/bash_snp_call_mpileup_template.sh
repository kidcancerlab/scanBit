#!/bin/sh
set -e ### stops bash script if line ends with error

placeholder_batch_other

conda activate scanBit_xkcd_1337

start_time=$(date +%s)

echo ${HOSTNAME} Beginning: $(date '+%Y-%m-%d %H:%M:%S')

cell_file_array=(placeholder_cell_files)

if [ ! -d placeholder_bcf_dir ]
then
    mkdir -p placeholder_bcf_dir
fi

# the bam files are written out with the filename determined by the label in
# your cellid_bam_table. We use this label to name the output files and to
# change the sample name during the mpileup call so that when we merge the bcfs
# each column will have a unique name

run_mpileup() {
    local cell_file=$1

    local label=${cell_file%%_cell_ids.txt}
    label=${label##*/}

    local orig_samp_name=$(samtools view -H placeholder_bam_file \
        | grep "SM:" \
        | head -n 1 \
        | perl -ne '/\tSM:(.+?)\t/; print $1')

    echo ${label} ${orig_samp_name}

    samtools view \
            -u \
            --threads 5 \
            --tag-file CB:${cell_file} \
            placeholder_bam_file \
        | samtools view \
            -u \
            --threads 5 \
            --tag xf:25 \
        | bcftools mpileup \
            --threads 5 \
            --max-depth 8000 \
            -a FORMAT/DP \
            -O u \
            -f placeholder_ref_fasta \
            -s "${orig_samp_name} ${label}" \
            - \
        | bcftools call \
            --threads 5 \
            -a GQ \
            -O u \
            -m \
            placeholder_ploidy \
        | bcftools filter \
            --threads 5 \
            -g 10 \
            -e "FORMAT/DP<placeholder_min_depth" \
            -O u \
        | bcftools view \
            --threads 5 \
            --exclude-types indels \
            -O b \
            -o placeholder_bcf_dir/${label}.bcf

    bcftools index \
        --threads 5 \
        placeholder_bcf_dir/${label}.bcf
}

export -f run_mpileup

# limiting to 20% of CPU count since I'm telling bcftools to use 5 cpus.
parallel -j 20% run_mpileup ::: ${cell_file_array[@]}

end_time=$(date +%s)

elapsed_seconds=$((end_time - start_time))

echo Done: $(date '+%Y-%m-%d %H:%M:%S')
echo Elapsed seconds: $elapsed_seconds

conda deactivate