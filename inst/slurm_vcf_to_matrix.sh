#!/bin/sh
#SBATCH --output=placeholder_job_log
#SBATCH --error=placeholder_job_log
#SBATCH --job-name scanBit_dist
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --mem=10G # Figure this out later
#SBATCH --wait
placeholder_job_header_other

placeholder_batch_other

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

conda activate scanBit_xkcd_1337

python \
    placeholder_py_script \
        --processes 10 \
        --bcf placeholder_bcf_input \
        --min_snvs_for_cluster placeholder_min_snvs \
        --max_prop_missing placeholder_max_missing \
        --n_bootstrap placeholder_n_bootstrap \
        --bootstrap_threshold placeholder_bootstrap_threshold \
        --figure_file placeholder_fig_file placeholder_verbose \
    > placeholder_groups_output
