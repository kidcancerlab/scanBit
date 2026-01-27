# SCAN-BIT
Single Cell Altered Nucleotide Based Inference of Tumor

---

## Purpose
This package provides functions to analyze single cell RNA-seq data to identify single nucleotide variants (SNVs) and use these to infer tumor phylogenies and identify tumor cells.

## Installation
You can install the development version of scanBit from GitHub with:
```R
# install.packages("devtools") # optional if you don't have devtools installed
devtools::install_github("kidcancerlab/scanBit")
```

After installation, I would recommend setting up the required conda environment.
```R
scanBit::confirm_conda_env()
```
If it works, it will return TRUE.

If the conda setup does not work, you will be prompted to try it with something like:
>It looks like you tried to make the conda environment with specific version requirements and it failed. Would you like to try making it with less stringent version requirements? please enter yes or no\n"

This will try to install the conda environment with looser version requirements.

If that fails for some reason, you can try to do it manually using some variant on the following code:
```bash
conda config --set channel_priority strict

conda create -y -n scanBit_xkcd_1337 python=3.12
conda activate scanBit_xkcd_1337

## scanBit_version _must_ match the version of scanBit you have installed in R
conda env config vars set scanBit_version=0.7.1

conda install -y 'samtools>=1.13' bcftools pysam numpy pandas matplotlib scipy
```
## Usage
See the vignette for detailed usage instructions:
https://github.com/kidcancerlab/scanBit/blob/master/vignettes/basics.Rmd

## Errors
If you run into any errors, please run the following commands and send the
output to me along with the error messages produced.

```
traceback()
sessionInfo()
```
## To do:
- Add in the option to print out a matrix showing how many variants are being used for hierarchical clustering
    - Update vignette with usage for this
