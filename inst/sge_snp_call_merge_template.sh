#!/bin/bash
#$ -cwd
#$ -j y
#$ -o placeholder_job_log
#$ -l mem_free=10G
#$ -q placeholder_sge_q                 # this needs to be user specified
#$ -pe thread 4             # this needs to be user specified
#$ -N placeholder_job_name
#$ -sync y                  # Wait for finish to continue
placeholder_job_header_other

placeholder_batch_other

set -e ### stops bash script if line ends with error

start_time=$(date +%s)

echo ${HOSTNAME} Beginning: $(date '+%Y-%m-%d %H:%M:%S')

conda activate scanBit_xkcd_1337

bcftools merge \
        --threads 4 \
        -O u \
        placeholder_bcf_dir/*.bcf \
    | bcftools view \
        -i 'GT[*]="alt"' \
        -O b \
        --output placeholder_bcf_out

bcftools index \
    --threads 4 \
    placeholder_bcf_out

end_time=$(date +%s)

elapsed_seconds=$((end_time - start_time))

echo Done: $(date '+%Y-%m-%d %H:%M:%S')
echo Elapsed: $(printf "%(%H:%M:%S)T" "$elapsed_seconds")

conda deactivate