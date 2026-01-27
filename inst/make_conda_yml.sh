#!/bin/bash
set -e

# conda remove -n scanBit_xkcd_1337 -y --all

# remove r conda channel
# conda config --remove channels r

# set conda channel priority to strict
conda config --set channel_priority strict

conda create -y -n scanBit_xkcd_1337 python=3.12
conda activate scanBit_xkcd_1337

conda env config vars set scanBit_version=1.0.0

conda install -y 'samtools>=1.13' bcftools pysam numpy pandas matplotlib scipy

conda env export --name scanBit_xkcd_1337 --no-builds | grep -v "^prefix"> inst/conda.yml

conda env export --name scanBit_xkcd_1337 --from-history | grep -v "^prefix"> inst/conda_general.yml
