#!/bin/bash
conda create -y -n rrrSNVs_xkcd_1337 python=3.12
conda activate rrrSNVs_xkcd_1337

pip install pysam numpy pandas

conda env export > inst/conda.yml
