#!/bin/bash
#$ -cwd
#$ -j y                                 # Error logging to single file
#$ -o placeholder_job_log
#$ -l mem_free=1G
#$ -q placeholder_sge_q
#$ -pe placeholder_sge_thread 10
#$ -N placeholder_job_name
#$ -sync y
placeholder_job_header_other

placeholder_batch_other

set -e ### stops bash script if line ends with error

start_time=$(date +%s)

echo ${HOSTNAME} Beginning: $(date '+%Y-%m-%d %H:%M:%S')

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

end_time=$(date +%s)

elapsed_seconds=$((end_time - start_time))

echo Done: $(date '+%Y-%m-%d %H:%M:%S')
echo Elapsed seconds: $elapsed_seconds

conda deactivate