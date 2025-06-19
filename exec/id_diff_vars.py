import argparse
import sys
import warnings
from pysam import VariantFile
import numpy as np

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--bcf',
                    type = str,
                    default=  '/home/gdrobertslab/mvc002/analyses/roberts/24_Osteo_atlas/output/id_tumor/snvs/human_X00003/mergedX00003_c30.bcf',
                    help = 'BCF file with multiple samples as columns')
parser.add_argument('--min_snvs_for_cluster',
                    type = int,
                    default = 250,
                    help = 'minimum number of SNVs for a cluster to be included')
parser.add_argument('--max_prop_missing',
                    type = float,
                    default = 0.9,
                    help = 'max proportion of missing data allowed at a single locus')
parser.add_argument('--group_1',
                    type = str,
                    default = 'cluster_0,cluster_2,cluster_3',
                    help = 'Comma delimited list of members of first group')
parser.add_argument('--group_2',
                    type = str,
                    default = 'cluster_1,cluster_4,cluster_7,cluster_5,cluster_6',
                    help = 'Comma delimited list of members of second group')
parser.add_argument('--out_file',
                    '-o',
                    type = str,
                    default = "variant_diff_table.tsv",
                    help = 'File to write out with variant differences by group')
parser.add_argument('--processes',
                    '-p',
                    type = int,
                    default = 1,
                    help = 'number of processes to use for parallel processing')
parser.add_argument('--verbose',
                    action = 'store_true',
                    help = 'print out extra information')

args = parser.parse_args()

################################################################################
### Code

def main():
    if args.verbose:
        print("Calculating distance matrix", file = sys.stderr)

    write_variant_diff_file(
        bcf_file = args.bcf,
        min_snvs_for_cluster = args.min_snvs_for_cluster,
        max_prop_missing = args.max_prop_missing
    )

    if args.verbose:
        print("Done!", file=sys.stderr)
    return()


########
### functions

def write_variant_diff_file(bcf_file,
                            min_snvs_for_cluster,
                            max_prop_missing):
    """
    Processes a BCF file to compute and output per-variant genotype
    difference statistics between two sample groups.

    Args:
        bcf_file (str): Path to the input BCF file containing variant data.
        min_snvs_for_cluster (int): Minimum number of SNVs required for a
            cluster to be considered.
        max_prop_missing (float): Maximum allowed proportion of missing
            genotype data per variant (0.0 - 1.0).

    Returns:
        differences (np.ndarray): 3D array of absolute genotype differences
            between all sample pairs for each variant.
        samples (tuple): Tuple of sample names included in the analysis.

    Side Effects:
        Writes a tab-delimited file 'test.txt' containing, for each variant:
            - Chromosome and position
            - Mean genotype difference between two predefined sample groups
              (cross_group_mean_diffs)
            - Mean genotype difference within the groups
              (within_group_mean_diffs)
            - The difference between cross-group and within-group means
              (delta_mean_diffs)

    Notes:
        - Assumes the existence of two sample groups.
        - Filters out variants with missing data above the specified threshold.
        - Requires external functions: pad_len_1_genotype and
          filter_diff_matrix.
        - Uses numpy and pysam.VariantFile for processing.
    """

    dist_key_dict = {'00':            0,
                     '01':            1,
                     '10':            1,
                     '11':            2,
                     '(None, None)':  np.nan}
    ### Check if bcf index exists
    bcf_in = VariantFile(bcf_file, threads = args.processes)
    samples = np.array(bcf_in.header.samples)
    records = tuple(x for x in list(bcf_in.fetch()) if len(x.alts) == 1)
    bcf_in.close()

    # Precompute the genotype tuples for all samples
    genotype_tuples = np.array([
        [tuple(pad_len_1_genotype(rec.samples[sample]['GT'])) for sample in samples]
        for rec in records
    ])

    genotype_locations = np.array([(rec.chrom, str(rec.pos)) for rec in records])
    alt_genotypes = np.array([rec.alts[0] for rec in records])

    # Convert genotype tuples to strings and look up in dist_key_dict
    genotype_matrix = np.array([
        [dist_key_dict.get(''.join(map(str, gt)), np.nan) for gt in sample_genotypes]
        for sample_genotypes in genotype_tuples
    ])

    # Filter out variant positions seen in less than x% of samples
    percent_missing = np.sum(np.isnan(genotype_matrix), axis=1) / len(samples)
    genotype_matrix = genotype_matrix[percent_missing <= max_prop_missing]
    genotype_locations = genotype_locations[percent_missing <= max_prop_missing]
    alt_genotypes = alt_genotypes[percent_missing <= max_prop_missing]

    # Filter out clusters/samples with too few variants called
    keep_sample_mask = (
        (
            genotype_matrix.shape[0]
            - np.apply_along_axis(np.isnan, 0, genotype_matrix).sum(axis = 0)
        )
        > min_snvs_for_cluster
    )
    genotype_matrix = genotype_matrix[:, keep_sample_mask]
    samples = samples[np.where(keep_sample_mask)]

    differences = np.abs(genotype_matrix[:, :, np.newaxis]
                         - genotype_matrix[:, np.newaxis, :])

    # convert the diagonal of the 2d array of each location to nan so when I
    # calculate the within-group mean it's not skewing results with extra zeros
    # I think this actually creates problems when one group only has a single
    # sample/cluster or when only one is present at a given genome position
    # differences[:, np.eye(len(samples), dtype = bool)] = np.nan

    print_combined_table(
        differences,
        samples,
        genotype_locations,
        alt_genotypes,
        genotype_matrix
    )

def parse_groups(samples):
    group_1_names = args.group_1.split(',')
    group_2_names = args.group_2.split(',')
    for this_name in group_1_names + group_2_names:
        if this_name not in samples:
            print(f'"{this_name}" not found in the bcl file samples after',
                  'filtering out samples with too much missing data!',
                  'You can increase max_prop_missing if you think this is wrong.')

    group_1 = np.array([list(samples).index(s) for s in group_1_names if s in samples])
    group_2 = np.array([list(samples).index(s) for s in group_2_names if s in samples])

    return group_1, group_2

def print_combined_table(differences,
                         samples,
                         genotype_locations,
                         alt_genotypes,
                         genotype_matrix):
    group_1, group_2 = parse_groups(samples)

    cross_group_mean_diffs = mean_along_axis(
        differences[:, group_1[:, None], group_2[None, :]],
        (1, 2))

    within_group_1_diffs = differences[:, group_1[:, None], group_1[None, :]]
    within_group_1_diffs = within_group_1_diffs.reshape(
        within_group_1_diffs.shape[0],
        -1)
    within_group_2_diffs = differences[:, group_2[:, None], group_2[None, :]]
    within_group_2_diffs = within_group_2_diffs.reshape(
        within_group_2_diffs.shape[0],
        -1)

    within_group_1_mean = mean_along_axis(within_group_1_diffs, 1)
    within_group_2_mean = mean_along_axis(within_group_2_diffs, 1)

    within_group_mean_diffs = mean_along_axis(
        np.concatenate(
            (within_group_1_diffs, within_group_2_diffs),
            axis = 1),
        1)

    delta_mean_diffs = cross_group_mean_diffs - within_group_mean_diffs

    mean_alts_group_1 = mean_along_axis(genotype_matrix[:, group_1], 1)
    mean_alts_group_2 = mean_along_axis(genotype_matrix[:, group_2], 1)

    # If you change the order of these columns make sure you update the column
    # info in *_cell_var_calls.sh and in cross_minus_within_dist_indx below
    combined_table = np.concatenate((genotype_locations,
                                     alt_genotypes[:, None],
                                     within_group_1_mean[:, None],
                                     within_group_2_mean[:, None],
                                     cross_group_mean_diffs[:, None],
                                     within_group_mean_diffs[:, None],
                                     delta_mean_diffs[:, None],
                                     mean_alts_group_1[:, None],
                                     mean_alts_group_2[:, None]),
                                     axis = 1)

    # Get rid of positions where the difference is nan
    combined_table = combined_table[~np.isnan(delta_mean_diffs)]

    # Rearrange the table by the delta_mean_diffs in descending order
    cross_minus_within_dist_indx = 7
    combined_table = combined_table[np.argsort(combined_table[:, cross_minus_within_dist_indx])[::-1]]

    np.savetxt(args.out_file,
               combined_table,
               delimiter = '\t',
               fmt = '%s',
               comments = '',
               header = ('chr\tpos\t'
                         + 'alt_gt\t'
                         + 'group_1_within_dist\t'
                         + 'group_2_within_dist\t'
                         + 'cross_group_dist\t'
                         + 'within_groups_mean_dist\t'
                         + 'cross_minus_within_dist\t'
                         + 'mean_group_1_alts\t'
                         + 'mean_group_2_alts'))

def mean_along_axis(differences_array, axis):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Mean of empty slice",
            category=RuntimeWarning
        )
        try:
            output = np.nanmean(differences_array, axis)
        except Exception as e:
            print(f"Calculating the mean of differences failed with error {e}")

    return output

def pad_len_1_genotype(gt):
    if len(gt) == 1:
        return (gt[0], 0)
    return gt

def filter_diff_matrix(differences, samples, min_snvs_for_cluster):
    n_snps_per_sample = np.sum(~np.isnan(differences), axis=0).diagonal()
    samples_to_keep = np.where(n_snps_per_sample >= min_snvs_for_cluster)[0].tolist()
    samples = [samples[i] for i in samples_to_keep]
    differences = differences[:, samples_to_keep, :]
    differences = differences[:, :, samples_to_keep]
    return differences, samples

################################################################################
### main

if __name__ == '__main__':
    main()
