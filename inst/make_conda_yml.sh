#!/bin/bash
set -e

# conda remove -n scanBit_xkcd_1337 -y --all

# remove r conda channel
# conda config --remove channels r

# set conda channel priority to strict
conda config --set channel_priority strict

mamba create -y -n scanBit_xkcd_1337 python=3.12
conda activate scanBit_xkcd_1337

conda env config vars set scanBit_version=0.5.0

mamba install -y samtools bcftools

pip install pysam numpy pandas matplotlib scipy

conda env export --name scanBit_xkcd_1337 --no-builds > inst/conda.yml
