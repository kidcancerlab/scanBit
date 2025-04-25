#!/bin/sh
#SBATCH --output=placeholder_job_log
#SBATCH --error=placeholder_job_log
#SBATCH --job-name=merge_bcfs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --wait
placeholder_job_header_other

placeholder_batch_other

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

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
