# scanBit 0.8.0
- Updated README with info on environment setup and usage.
- Updated vignette with detailed instructions on environment setup.
- Changed file extension for printing trees to .pdf with argument to specify what is wanted.
- Added some code to identify which variants are most variant across specified populations of cells. I also put in code to call variants at just these locations across individual cells to speed up the process. Just a word of warning, this is still fairly slow and most of your cells will lack calls at many of these positions due to low coverage.

# scanBit 0.7.1
- Fixed a bug with vcfToMatrix.py that I found immediately after the last release. I had moved the argparse code and `args` was no longer global.

# scanBit 0.7.0
- Changed the underlying code for how batch jobs are submitted to allow for slurm, sge or bash jobs.
- Updated documentation and vignette.

# scanBit 0.6.0
- Completely refactored how I pull the reads for each cell from the bam file. I'm now using samtools directly in the pipe to do variant calling which speeds things up a lot.

# scanBit 0.5.0
- Finished work to rename the package
- Completely refactored the code for bootstrapping in vcfToMatrix.py. The old version was not working properly.

# rrrSnvs 0.4.0
- Added option to specify other sbatch options
- Updated add_snv_group_to_sobj and label_tumor_cells to handle the case when the code does not create a *dist.txt file due to insufficient cells, reads, etc..
- Updated vignette

# rrrSnvs 0.3.0
- Modified the distance calculation and output to report both all groups as well as the top level split for each tree.
- Added a new plotting function assemble_tree_plots() to read in the pngs written to the output folder and assemble them into a single plot. This will eventually be used to also add in dimplots and other visualizations for a final output figure.
- Updated the conda environment to include bcftools and samtools.

# rrrSnvs 0.2.0
- Updated the conda environment to include all the necessary packages for the
package.

# rrrSnvs v0.1.0
- Made a new repo for this code and stripped out the rrrSingleCellUtils code.
