#!/bin/sh
#SBATCH --account=placeholder_account
#SBATCH --output=placeholder_slurm_out
#SBATCH --error=placeholder_slurm_out
#SBATCH --job-name rrrSNVs_dist
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --mem=10G # Figure this out later
#SBATCH --partition=himem,general
#SBATCH --wait
placeholder_sbatch_other

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

ml purge
ml load miniforge3/24.3.0

eval "$(conda shell.bash hook)"
conda activate rrrSNVs_xkcd_1337

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
