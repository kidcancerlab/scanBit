#!/bin/bash
set -e

# conda remove -n rrrSNVs_xkcd_1337 --all

conda create -y -n rrrSNVs_xkcd_1337 python=3.12
conda activate rrrSNVs_xkcd_1337

conda env config vars set rrrSnvs_version=0.2.0

pip install pysam numpy pandas matplotlib scipy

conda env export > inst/conda.yml
