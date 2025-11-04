#!/bin/sh
#SBATCH --account=placeholder_account
#SBATCH --output=placeholder_job_log
#SBATCH --error=placeholder_job_log
#SBATCH --job-name scanBit_var
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --wait
placeholder_job_header_other

placeholder_batch_other

set -e ### stops bash script if line ends with error

start_time=$(date +%s)

echo ${HOSTNAME} Beginning: $(date '+%Y-%m-%d %H:%M:%S')

conda activate scanBit_xkcd_1337

python \
    placeholder_py_script \
        --processes 2 \
        --bcf placeholder_bcf_input \
        --group_1 placeholder_group_1 \
        --group_2 placeholder_group_2 \
        --min_snvs_for_cluster placeholder_min_snvs \
        --max_prop_missing placeholder_max_missing \
        --out_file placeholder_out_file

end_time=$(date +%s)

elapsed_seconds=$((end_time - start_time))

echo Done: $(date '+%Y-%m-%d %H:%M:%S')
echo Elapsed seconds: $elapsed_seconds

conda deactivate
