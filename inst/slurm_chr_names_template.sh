#!/bin/sh
#SBATCH --output=placeholder_job_log
#SBATCH --error=placeholder_job_log
#SBATCH --job-name=chr_names
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --wait
placeholder_job_header_other

placeholder_batch_other

set -e ### stops bash script if line ends with error

conda activate scanBit_xkcd_1337

chr_file=placeholder_chr_file

bam_file_array=(placeholder_bam_files)

# we use index file to speed this up
ref_file=placeholder_ref_fasta
ref_file=${ref_file}.fai

start_time=$(date +%s)

echo ${HOSTNAME} Beginning: $(date '+%Y-%m-%d %H:%M:%S')

echo input: ${bam_file_array[@]} ${ref_file}, output: ${chr_file}

# Get fasta headers
echo -ne "ref\t" > "${chr_file}"
cut -f 1 \
  "${ref_file}" \
  | perl -pe 's/\n/,/' \
  | perl -pe 's/,$/\n/' \
  >> "${chr_file}"

for this_bam_file in "${bam_file_array[@]}"
do
  echo -ne "${this_bam_file}\t" >> "${chr_file}"

  samtools idxstats \
      "${this_bam_file}" \
    | cut -f 1 \
    | grep -v "^*" \
    | perl -pe 's/\n/,/' \
    | perl -pe 's/,$/\n/' \
  >> "${chr_file}"
done

end_time=$(date +%s)

elapsed_seconds=$((end_time - start_time))

echo Done: $(date '+%Y-%m-%d %H:%M:%S')
echo Elapsed seconds: $elapsed_seconds

conda deactivate