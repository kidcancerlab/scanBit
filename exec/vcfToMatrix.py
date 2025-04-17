import argparse
import sys
from pysam import VariantFile
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import multiprocessing
from itertools import repeat, chain

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--bcf',
                    type = str,
                    default=  '/home/gdrobertslab/mvc002/analyses/roberts/24_Osteo_atlas/output/id_tumor/snvs/mouse_S0169/mergedS0169_c5.bcf',
                    help = 'BCF file with multiple samples as columns')
parser.add_argument('--figure_file',
                    '-o',
                    type = str,
                    default = 'dendrogram.pdf',
                    help = 'output file name of plot. Id suggest either png or pdf')
parser.add_argument('--min_snvs_for_cluster',
                    type = int,
                    default = 1000,
                    help = 'minimum number of SNVs for a cluster to be included')
parser.add_argument('--max_prop_missing',
                    type = float,
                    default = 0.9,
                    help = 'max proportion of missing data allowed at a single locus')
parser.add_argument('--n_bootstrap',
                    type = int,
                    default = 1000,
                    help = 'number of bootstrap samples to use')
parser.add_argument('--bootstrap_threshold',
                    type = float,
                    default = 0.99,
                    help = 'threshold for collapsing clusters')
parser.add_argument('--verbose',
                    action = 'store_true',
                    help = 'print out extra information')
parser.add_argument('--processes',
                    '-p',
                    type = int,
                    default = 1,
                    help = 'number of processes to use for parallel processing')
parser.add_argument('--fig_width',
                    type = float,
                    default = 6,
                    help = 'width of the figure in inches')
parser.add_argument('--fig_height',
                    type = float,
                    default = 6,
                    help = 'height of the figure in inches')
parser.add_argument('--fig_dpi',
                    type = int,
                    default = 300,
                    help = 'dpi of the figure')

args = parser.parse_args()

################################################################################
### Code

def main():
    if args.verbose:
        print("Calculating distance matrix", file = sys.stderr)
    differences, samples = get_diff_matrix_from_bcf(
        bcf_file = args.bcf,
        min_snvs_for_cluster = args.min_snvs_for_cluster,
        max_prop_missing = args.max_prop_missing)

    prop_diff_matrix = calc_proportion_dist_matrix(differences)

    hclust_out = hierarchical_clustering(prop_diff_matrix)

    original_cluster_dict = get_cluster_dict(hclust_out, samples)

    if args.verbose:
        print("Performing bootstrapping", file=sys.stderr)
    # What I want to get back from this function:
    # A proportion for each node in original_cluster_dict.keys()
    # Internally:
    # An array of 1/0 where the value indicates if that node existed in the bootstrap and had the appropriate sub-clusters
    # Use vstack to stack these into a 2D numpy array
    # Each sample will be a row
    # The columns of the array should be ordered according to original_cluster_dict.keys() to allow for downstream analysis
    bootstrap_values, node_members = \
        calculate_bootstrap_values(
            original_cluster_dict,
            differences,
            samples,
            args.n_bootstrap,
            threads = args.processes
        )

    plt.figure(figsize=(args.fig_width, args.fig_height))
    plot = plot_dendro_with_bootstrap_values(hclust_out, bootstrap_values, samples)
    plot.tight_layout()
    plot.savefig(args.figure_file, dpi = args.fig_dpi)

    collapsed_clusters = collapse_clusters(node_members,
                                           bootstrap_values,
                                           threshold=args.bootstrap_threshold)

    top_lvl_clusters = collapse_top_lvl_clusters(node_members,
                                                 bootstrap_values,
                                                 threshold=args.bootstrap_threshold)

    print_cluster_names(collapsed_clusters, top_lvl_clusters)

    if args.verbose:
        print("Done!", file=sys.stderr)
    return()


########
### functions

####### Can I clear out genotypes with all missing data? #################
def get_diff_matrix_from_bcf(bcf_file,
                             min_snvs_for_cluster,
                             max_prop_missing):
    dist_key_dict = {'00':            0,
                     '01':            1,
                     '10':            1,
                     '11':            2,
                     '(None, None)':  np.nan}
    ### Check if bcf index exists
    bcf_in = VariantFile(bcf_file, threads = args.processes)
    samples = tuple(bcf_in.header.samples)
    records = tuple(x for x in list(bcf_in.fetch()) if (len(x.alts) == 1))
    bcf_in.close()

    # Precompute the genotype tuples for all samples
    genotype_tuples = np.array([
        [tuple(pad_len_1_genotype(rec.samples[sample]['GT'])) for sample in samples]
        for rec in records
    ])

    # Convert genotype tuples to strings and look up in dist_key_dict
    genotype_matrix = np.array([
        [dist_key_dict.get(''.join(map(str, gt)), np.nan) for gt in sample_genotypes]
        for sample_genotypes in genotype_tuples
    ])

    # Filter out variant positions seen in less than x% of samples
    percent_missing = np.sum(np.isnan(genotype_matrix), axis=1) / len(samples)
    genotype_matrix = genotype_matrix[percent_missing <= max_prop_missing]
    differences = np.abs(genotype_matrix[:, :, np.newaxis]
                         - genotype_matrix[:, np.newaxis, :])
    differences, samples = filter_diff_matrix(differences,
                                              samples,
                                              min_snvs_for_cluster)

    if differences.shape[1] == 0:
        sys.exit(
            f'No clusters with at least {min_snvs_for_cluster} SNVs for {args.bcf}'
        )

    return differences, samples

def pad_len_1_genotype(gt):
    if len(gt) == 1:
        return (gt[0], 0)
    else:
        return gt

def filter_diff_matrix(differences, samples, min_snvs_for_cluster):
    n_snps_per_sample = np.sum(~np.isnan(differences), axis=0).diagonal()
    samples_to_keep = np.where(n_snps_per_sample >= min_snvs_for_cluster)[0].tolist()
    samples = [samples[i] for i in samples_to_keep]
    differences = differences[:, samples_to_keep, :]
    differences = differences[:, :, samples_to_keep]
    return differences, samples

def calc_proportion_dist_matrix(differences, bootstrap=False):
    if bootstrap:
        differences = differences[np.random.choice(differences.shape[0],
                                                   size=differences.shape[0],
                                                   replace=True)]
    n_comps_matrix = np.sum(~np.isnan(differences), axis=0)
    # Sum up differences while ignoring np.nan values
    sum_differences = np.nansum(differences, axis=0)
    # Calculate the proportion of differences
    prop_diff_matrix = sum_differences / (n_comps_matrix * 2)
    return prop_diff_matrix

def hierarchical_clustering(distance_matrix,
                            linkage_method='ward'):
    distance_matrix = squareform(distance_matrix)
    Z = linkage(distance_matrix, method=linkage_method)
    return Z

class cluster_hierarchy:
    def __init__(self, cluster_no, cluster_member_nums,
                 cluster_member_names, sub_cluster_nums,
                 sub_cluster_names):
        self.cluster_no = cluster_no
        self.cluster_member_nums = cluster_member_nums
        self.cluster_member_names = cluster_member_names
        self.sub_cluster_numbers = sub_cluster_nums
        self.sub_cluster_names = sub_cluster_names

def build_initial_cluster_dict(samples):
    n_samples = len(samples)
    cluster_dict = {}

    # Initialize the cluster_dict with the first n_samples clusters
    for i in range(n_samples):
        cluster_dict[samples[i]] = cluster_hierarchy(
            cluster_no=i,
            cluster_member_nums=[i],
            cluster_member_names=[samples[i]],
            sub_cluster_nums=[i],
            sub_cluster_names=[samples[i]]
        )
    return cluster_dict

def get_cluster_dict(hclust, samples):
    n_samples = len(samples)

    # going to create a dict, where the key is a pasted list of sorted cluster
    # names and the value is a cluster_hierarchy object
    # I can use this to create a mapping of cluster names to their members
    cluster_dict = build_initial_cluster_dict(samples)
    # Now do the more complex clusters
    # In this context, a node is a point in the hierarchy where two clusters split
    n_nodes = hclust.shape[0]
    node_names = list(samples)
    this_cluster_num = n_samples
    for i in range(n_nodes):
        # numbers and names of the two clusters under this node
        node_member_nums = hclust[i, :2].astype(int).tolist()
        node_member_names = [node_names[i] for i in node_member_nums]
        # get sample names of all samples below this node
        subcluster_names = \
            (
                cluster_dict.get(node_member_names[0]).cluster_member_names,
                cluster_dict.get(node_member_names[1]).cluster_member_names
            )
        # get number ids of all samples beneath this node
        subcluster_numbers = \
            (
                cluster_dict.get(node_member_names[0]).cluster_member_nums,
                cluster_dict.get(node_member_names[1]).cluster_member_nums
            )

        all_names = sorted(chain(*subcluster_names))
        all_numbers = sorted(chain(*subcluster_numbers))

        # here we concatenate the sorted member names with "-" to make single string to use as key
        cluster_dict['-'.join(all_names)] = cluster_hierarchy(
            cluster_no=this_cluster_num,
            cluster_member_nums=all_numbers,
            cluster_member_names=all_names,
            sub_cluster_nums=subcluster_numbers,
            sub_cluster_names=subcluster_names
        )
        this_cluster_num += 1
        node_names.append('-'.join(all_names))

    return cluster_dict

def bootstrap_worker(rand_seed, differences, true_cluster_dict, samples):
    np.random.seed(rand_seed)
    prop_diff_matrix_boot = calc_proportion_dist_matrix(differences,
                                                        bootstrap = True)

    hclust_out_boot = hierarchical_clustering(prop_diff_matrix_boot)

    boot_cluster_dict = get_cluster_dict(hclust_out_boot, samples)

    # This gets a numpy array ordered by keys in true_cluster_dict
    # For each conserved node in true_cluster_dict, the array has 1, otherwise 0
    call_vals = compare_orig_boot_dicts(true_cluster_dict, boot_cluster_dict)

    return call_vals

def get_subcluster_set(cluster_dict, node_name):
    subclusters = cluster_dict.get(node_name).sub_cluster_names
    subclusters_cat = ["-".join(sorted(x)) for x in subclusters]
    return(set(subclusters_cat))

def compare_orig_boot_dicts(true_cluster_dict, boot_cluster_dict):
    node_names = list(true_cluster_dict.keys())

    conserved_nodes = np.zeros(shape = len(node_names))
    for i in range(len(node_names)):
        this_node_name = node_names[i]
        if this_node_name in boot_cluster_dict.keys():
            true_subclusters_set = get_subcluster_set(true_cluster_dict,
                                                      this_node_name)

            boot_subclusters_set = get_subcluster_set(boot_cluster_dict,
                                                      this_node_name)

            # use "set()" so I can get for equivalence even if they're not ordered the same
            if true_subclusters_set == boot_subclusters_set:
                conserved_nodes[i] = 1

    return conserved_nodes


# Calculate bootstrap values for each node
def calculate_bootstrap_values(true_cluster_dict,
                               differences,
                               samples,
                               n_bootstraps,
                               threads):
    with multiprocessing.Pool(processes=threads) as pool:
        bootstrap_conserved_node_vals = pool.starmap(
            bootstrap_worker,
            zip(
                range(n_bootstraps),
                repeat(differences),
                repeat(true_cluster_dict),
                repeat(samples)
            )
        )

    bootstrap_2d_array = np.vstack(bootstrap_conserved_node_vals)

    # Sum each "column" which is each node
    prop_diff = bootstrap_2d_array.sum(axis = 0) / n_bootstraps
    # Force the tip nodes (nodes consisting of a single sample) to have bootstrap 0
    prop_diff[:len(samples)] = 0

    # node_number_matching_keys = \
    #     [true_cluster_dict.get(x).cluster_no for x in true_cluster_dict.keys()]#[len(samples):]


    node_members_matching_keys = \
        [true_cluster_dict.get(x).sub_cluster_names for x in true_cluster_dict.keys()]#[len(samples):]

    return prop_diff, node_members_matching_keys

def plot_dendro_with_bootstrap_values(hclust,
                                      bootstrap_values,
                                      samples):
    # Don't want bootstrap values for tips here
    bootstrap_values = bootstrap_values[len(samples):]

    # The nodes in bootstrap_values are ordered by min to max distance
    # They also exclude the values from the tips (one sample)
    dend_plot = dendrogram(hclust,
                           labels=samples,
                           leaf_rotation=90.,
                           leaf_font_size=6.,
                           color_threshold=0)
    icoords = dend_plot['icoord']
    dcoords = dend_plot['dcoord']

    # Need to sort the coordinates in the dendrogram so that they match the order of in hclust and bootstrap_values
    # This sorts the coordinates by the y value low -> high
    right_order = np.argsort([x[1] for x in dcoords])
    icoords = [icoords[x] for x in right_order]
    dcoords = [dcoords[x] for x in right_order]

    for i, (icoord, dcoord) in enumerate(zip(icoords, dcoords)):
        x = ((icoord[1] + icoord[2]) * 0.5)
        y = dcoord[1]
        support = bootstrap_values[i] * 100 # original data is proportion
        plt.text(x, y, f'{support:.2f}%', va='bottom', ha='center', fontsize=6)

    return(plt)

def make_parent_key(true_clusters):
    parent_key = {}
    for i in range(len(true_clusters)):
        if (len(true_clusters[i]) == 1):
            # Putting this in list so it's not passed by reference
            this_group = list(true_clusters[i])
        else:
            this_group = sorted(true_clusters[i][0] + true_clusters[i][1])
        parent_group = [x for x in range(len(true_clusters)) if this_group in true_clusters[x]]
        # We want this to be a single number, not a list
        parent_group = parent_group[0] if len(parent_group) > 0 else ''
        parent_key[i] = parent_group
    return parent_key

def collapse_clusters(true_clusters,
                      bootstrap_values,
                      threshold = args.bootstrap_threshold):
    # This is a dict where the key is cluster number and value is parent number
    parent_key = make_parent_key(true_clusters)

    # Store which parent nodes are forks here
    parent_dict = {}
    # Store groups here
    group_list = []

    root_node_num = len(bootstrap_values) - 1
    # loop over each cluster, starting with the largest
    # We're going to build a dict where the lowest nodes where all higher nodes
    # have bootstrap > threshhold is going to define groups
    # We move down the tree from the top
    # Where we have a node where the parent is a fork (bootstrap > threshhold)
    # and the bootstrap value is not, all tips below that are a group
    for cluster_num in range(len(bootstrap_values) - 1, -1, -1):
        # want to only consider where it's the root node or the parent node has valid bootstrap value
        if cluster_num == root_node_num or parent_key.get(cluster_num) in parent_dict:
            if bootstrap_values[cluster_num] >= threshold:
                # this node is a fork, so the downstream tips are not a group
                parent_dict[cluster_num] = "fork"
                # # if there are only tips below this node (sub-lists are all len of 1)
                ####if max([len(x) for x in true_clusters[cluster_num]]) == 1:
                ####    group_list.extend(true_clusters[cluster_num])
            else:
                # this node is not a fork, but it's parent is, so downstream tips are all a single group
                # We merge the two sub_groups below this node
                if len(true_clusters[cluster_num]) == 1:
                    group_list.append(true_clusters[cluster_num])
                else:
                    group_list.append(sorted(true_clusters[cluster_num][0] + \
                                             true_clusters[cluster_num][1]))

    # Keep only keys from "cluster" nodes in parent_dict
    clusters = list(group_list)

    return clusters

def collapse_top_lvl_clusters(true_clusters,
                              bootstrap_values,
                              threshold = args.bootstrap_threshold):
    top_lvl = len(bootstrap_values) - 1
    # top_lvl bootstrap value is the first division
    # We want either one or two lists inside another list
    top_clusters = []
    if bootstrap_values[top_lvl] >= threshold:
        top_clusters = list(true_clusters[top_lvl])
    else:
        top_clusters.append(true_clusters[top_lvl][0] + true_clusters[top_lvl][1])

    return top_clusters

def print_cluster_names(collapsed_clusters, top_lvl_clusters):
    to_print_dict = {}
    # We're going to put both datasets into the same dictionary using the sample
    # name as the key with sub-dictionaries for all_groups and top_lvl_groups
    # Then we will loop over the dictionary to print out the group names for
    # each cluster name
    for i in range(len(collapsed_clusters)):
        for sample in collapsed_clusters[i]:
            to_print_dict[sample] = {'all_groups': 'group' + str(i)}

    for i in range(len(top_lvl_clusters)):
        for sample in top_lvl_clusters[i]:
            to_print_dict[sample]['top_lvl_groups'] = 'top_lvl_group' + str(i)

    for this_key in to_print_dict.keys():
        print(f'{this_key}\t{to_print_dict[this_key]["all_groups"]}\t{to_print_dict[this_key]["top_lvl_groups"]}')

################################################################################
### main

if __name__ == '__main__':
    main()
