#!/bin/bash

placeholder_batch_other

set -e ### stops bash script if line ends with error

start_time=$(date +%s)

echo ${HOSTNAME} Beginning: $(date '+%Y-%m-%d %H:%M:%S')

placeholder_conda_prefixconda activate scanBit_xkcd_1337

placeholder_apptainerpython \
  placeholder_py_script \
    --processes 2 \
    --bcf placeholder_bcf_input \
    --group_1 placeholder_group_1 \
    --group_2 placeholder_group_2 \
    --min_snvs_for_cluster placeholder_min_snvs \
    --max_prop_missing placeholder_max_missing \
    --out_file placeholder_out_fileplaceholder_end_apptainer

end_time=$(date +%s)

elapsed_seconds=$((end_time - start_time))

echo Done: $(date '+%Y-%m-%d %H:%M:%S')
echo Elapsed seconds: $elapsed_seconds

placeholder_conda_prefixconda deactivate
