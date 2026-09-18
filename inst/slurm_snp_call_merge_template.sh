#!/bin/bash
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

start_time=$(date +%s)

echo ${HOSTNAME} Beginning: $(date '+%Y-%m-%d %H:%M:%S')

placeholder_conda_prefixconda activate scanBit_xkcd_1337

merge_bcfs() {
  bcftools merge \
        --threads 3 \
        -O u \
        placeholder_bcf_dir/*.bcf \
    | bcftools view \
        -i 'GT[*]="alt"' \
        -O b \
    --output placeholder_bcf_out

bcftools index \
    --threads 3 \
    placeholder_bcf_out
}

export -f merge_bcfs

touch placeholder_bcf_out

placeholder_apptainermerge_bcfsplaceholder_end_apptainer

end_time=$(date +%s)

elapsed_seconds=$((end_time - start_time))

echo Done: $(date '+%Y-%m-%d %H:%M:%S')
echo Elapsed seconds: $elapsed_seconds

placeholder_conda_prefixconda deactivate