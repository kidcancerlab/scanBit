#!/bin/bash
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

start_time=$(date +%s)

echo ${HOSTNAME} Beginning: $(date '+%Y-%m-%d %H:%M:%S')

placeholder_conda_prefixconda activate scanBit_xkcd_1337

touch placeholder_n_comps_file
touch placeholder_fig_file
touch placeholder_groups_output

placeholder_apptainerpython \
    placeholder_py_script \
        --processes 10 \
        --bcf placeholder_bcf_input \
        --min_snvs_for_cluster placeholder_min_snvs \
        --max_prop_missing placeholder_max_missing \
        --n_bootstrap placeholder_n_bootstrap \
        --bootstrap_threshold placeholder_bootstrap_threshold \
        --linkage_method placeholder_linkage_method \
        --n_comps_file placeholder_n_comps_file \
        --dist_method placeholder_dist_method \
        --figure_file placeholder_fig_file placeholder_verbose \
    > placeholder_groups_outputplaceholder_end_apptainer

end_time=$(date +%s)

elapsed_seconds=$((end_time - start_time))

echo Done: $(date '+%Y-%m-%d %H:%M:%S')
echo Elapsed seconds: $elapsed_seconds

conda deactivate