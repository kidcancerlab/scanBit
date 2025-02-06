#!/bin/bash
set -e

# conda remove -n rrrSNVs_xkcd_1337 -y --all

# remove r conda channel
# conda config --remove channels r

# set conda channel priority to strict
conda config --set channel_priority strict

mamba create -y -n rrrSNVs_xkcd_1337 python=3.12
conda activate rrrSNVs_xkcd_1337

conda env config vars set rrrSnvs_version=0.3.0

mamba install -y samtools bcftools

pip install pysam numpy pandas matplotlib scipy

conda env export --no-builds > inst/conda.yml
