#' Generate a tree using SNPs derived from single cell clusters
#'
#' @param cellid_bam_table A tibble with three columns: cell_barcode,
#'   cell_group and bam_file. The cell_barcode column should contain the cell
#'   barcode, the cell_group column should contain the cluster label and the
#'   bam_file column should contain the path to the bam file for that cell.
#' @param temp_dir The directory to write temporary files to.
#' @param output_dir The directory to write output distance files to.
#' @param output_base_name The prefix to use with the output files.
#' @param slurm_base The directory to write slurm output files to.
#' @param sbatch_base The prefix to use with the sbatch job file.
#' @param account The hpc account to use.
#' @param ploidy Either a file path to a ploidy file, or a string indicating
#'   the ploidy. #### PUT EXAMPLES HERE #### For ploidy file format see bcftools
#'   call --ploidy-file option description here:
#'   "https://samtools.github.io/bcftools/bcftools.html"
#' @param ref_fasta The path to the reference fasta file.
#' @param min_depth The minimum depth to use when calling SNPs. Can be a vector
#'   of numbers.
#' @param min_snvs_per_cluster The minimum number of sites that must be covered by
#'   a cell_group to be included in the tree.
#' @param max_prop_missing_at_site The maximum proportion of missing data
#'   allowed at a site.
#' @param n_bootstraps The number of bootstraps to use when evaluating grouping
#'   of samples in the hierarchical clustering tree.
#' @param bootstrap_cutoff The bootstrap cutoff to use when collapsing the tree.
#'   This is the proportion of bootstraps that must support a grouping for it to
#'   be considered a valid grouping.
#' @param verbose Whether to print out verbose output or not.
#' @param submit Whether to submit the sbatch jobs to the cluster or not.
#' @param cleanup Whether to clean up the temporary files after execution.
#'
#' @details for ploidy, GRCh37 is hg19, GRCh38 is hg38, X, Y, 1, mm10_hg19 is
#'   our mixed species reference with species prefixes on chromosomes, mm10 is
#'   mm10
#'
#' @return TRUE if the function succeeded. Output is written to output_dir.
#' @export
#'
#' @examples
#' \dontrun{
#' placeholder for now
#' }
get_snp_tree <- function(cellid_bam_table,
                         temp_dir = tempdir(),
                         output_dir,
                         output_base_name = "hier_tree",
                         slurm_base = "/slurmOut",
                         sbatch_base = "sbatch_",
                         account = "gdrobertslab",
                         ploidy,
                         ref_fasta,
                         min_depth = 5,
                         min_snvs_per_cluster = 500,
                         max_prop_missing_at_site = 0.9,
                         n_bootstraps = 1000,
                         bootstrap_cutoff = 0.99,
                         verbose = TRUE,
                         submit = TRUE,
                         cleanup = TRUE) {
    check_cellid_bam_table(cellid_bam_table)

    ## Check that the conda command is available
    check_cmd("conda")
    # Check that conda environment rrrSNVs_xkcd_1337 exists, and if not,
    # create it
    confirm_conda_env()

    # Check that temp_dir exists, and if not, create it
    if (!dir.exists(temp_dir)) {
        dir.create(temp_dir)
    }
    message("Using temporary directory: ", temp_dir)


    # Warn that this is going to take a while
    message("Hold onto your hat and get a coffee, this will take a while.")

    bam_files <- unique(cellid_bam_table$bam_file)

    results <-
        parallel::mclapply(seq_len(length(bam_files)),
                           mc.cores = length(bam_files),
                           mc.preschedule = FALSE,
                           function(i) {
            bam_name <- bam_files[i]

            sub_cellid_bam_table <-
                cellid_bam_table %>%
                dplyr::filter(bam_file == bam_name) %>%
                dplyr::select(cell_barcode, cell_group)

            call_snps(
                cellid_bam_table = sub_cellid_bam_table,
                bam_to_use = bam_name,
                bam_out_dir = paste0(
                    temp_dir,
                    "/split_bams_",
                    i,
                    "/"
                ),
                bcf_dir = paste0(
                    temp_dir,
                    "/split_bcfs_",
                    i
                ),
                slurm_base = slurm_base,
                sbatch_base = sbatch_base,
                account = account,
                ploidy = ploidy,
                min_depth = min_depth,
                ref_fasta = ref_fasta,
                temp_dir = temp_dir,
                submit = submit,
                cleanup = cleanup
            )
        })

    # The output from the previous step is a folder for each bam file located
    # in temp_dir/split_bcfs_{i}_c{min_depth}/. Next merge all the bcf files and
    # calculate a distance matrix using my slow python script
    # We do this separately for each min_depth provided

    parallel::mclapply(min_depth,
                       mc.cores = 100,
                       mc.preschedule = FALSE,
                       function(this_min_depth) {
        merge_bcfs(
            bcf_in_dir = paste0(temp_dir,
                                "/split_bcfs_[0-9]*_c",
                                this_min_depth,
                                "/"),
            out_bcf = paste0(output_dir,
                            "/merged",
                            output_base_name,
                            "_c",
                            this_min_depth,
                            ".bcf"),
            submit = submit,
            slurm_base = paste0(slurm_base, "_merge-%j.out"),
            sbatch_base = sbatch_base,
            account = account,
            temp_dir = temp_dir,
            cleanup = cleanup
        )
    })

    # Read in the merged bcf and make a tree for each min_depth
    # Number of cores doesn't matter here, we're just submitting slurm jobs
    parallel::mclapply(min_depth,
                       mc.cores = 101,
                       mc.preschedule = FALSE,
                       function(this_min_depth) {
        group_clusters_by_dist(
            merged_bcf_file =
                paste0(
                    output_dir,
                    "/merged",
                    output_base_name,
                    "_c",
                    this_min_depth,
                    ".bcf"
                ),
            output_file =
                paste0(
                    output_dir, "/",
                    output_base_name,
                    "_",
                    this_min_depth,
                    "_dist.txt"
                ),
            min_snvs_per_cluster = min_snvs_per_cluster,
            max_prop_missing_at_site = max_prop_missing_at_site,
            n_bootstraps = n_bootstraps,
            bootstrap_cutoff = bootstrap_cutoff,
            tree_figure_file =
                paste0(
                    output_dir, "/",
                    output_base_name,
                    "_",
                    this_min_depth,
                    "_tree.png"
                ),
            verbose = verbose,
            sbatch_base = sbatch_base,
            account = account,
            slurm_base = slurm_base,
            temp_dir = temp_dir,
            submit = submit
        )
    })

    return(TRUE)
}

#' Call SNPs for a single bam file
#'
#' @param cellid_bam_table A table with columns cell_id, cell_group and bam_file
#' @param bam_to_use The bam file to use
#' @param bam_out_dir The directory to write the bam files to
#' @param bcf_dir The directory to write the bcf files to
#' @param slurm_base The base name for the slurm output files. Don't include path
#' @param sbatch_base The prefix to use with the sbatch job file
#' @param account The hpc account to use in slurm scripts
#' @param ploidy The ploidy of the organism. See details for more information
#' @param ref_fasta The reference fasta file to use
#' @param min_depth The minimum depth to use when calling snps
#' @param temp_dir The directory to write temporary files to
#' @param submit Whether or not to submit the slurm scripts
#' @param cleanup Whether or not to clean up the bam files afterwards
#'
#' @return 0 if the snps were called successfully
#'
#' @details GRCh37 is hg19, GRCh38 is hg38, X, Y, 1, mm10_hg19 is our mixed
#'
#' @noRd
call_snps <- function(cellid_bam_table,
                      bam_to_use,
                      bam_out_dir,
                      bcf_dir,
                      slurm_base = paste0(getwd(), "/slurmOut_call-%j.txt"),
                      sbatch_base = "sbatch_",
                      account = "gdrobertslab",
                      ploidy,
                      ref_fasta,
                      min_depth = 5,
                      temp_dir,
                      submit = TRUE,
                      cleanup = TRUE) {
    # Check that the bam and bcf directories exist, and if not, create them
    if (!dir.exists(bam_out_dir)) {
        dir.create(bam_out_dir)
    }

    # write out the cell ids to a file with two columns: cell_id, cell_group
    cell_file <- paste0(bam_out_dir, "/cell_ids.txt")
    readr::write_tsv(dplyr::select(cellid_bam_table, cell_barcode, cell_group),
                     file = cell_file,
                     col_names = FALSE,
                     progress = FALSE)

    # call getBarcodesFromBam.py on the bam file and the cell id file by reading
    # in a template and substituting out the placeholder fields
    # write the bam files to a subfolder for each source bam
    py_file <-
        paste0(find.package("rrrSnvs"),
               "/exec/getBarcodesFromBam.py")

    replace_tibble_split <-
        dplyr::tribble(
            ~find,                      ~replace,
            "placeholder_account",      account,
            "placeholder_slurm_out",    paste0(
                                            temp_dir, "/",
                                            slurm_base,
                                            "_splitbams-%j.out"
                                        ),
            "placeholder_cell_file",    cell_file,
            "placeholder_bam_file",     bam_to_use,
            "placeholder_bam_dir",      paste0(bam_out_dir, "/"),
            "placeholder_py_file",      py_file
        )

    # We're going to split the bam file and store the output in
    # bam_out_dir/split_bams/. We'll then call mpileup on each of these
    # bam_out_dir is passed by get_snp_tree() and is the temp_dir/split_bams_{i}/
    result <-
        use_sbatch_template(replace_tibble_split,
                            "snp_call_splitbams_template.job",
                            warning_label = "Bam splitting",
                            submit = submit,
                            file_dir = temp_dir,
                            temp_prefix = paste0(sbatch_base, "split_"))

    ploidy <- pick_ploidy(ploidy)

    array_max <-
        list.files(path = bam_out_dir, pattern = ".bam") %>%
        length() - 1

    # Since min_depth can be a vector of unknown length, we are going to loop
    # through each element and call mpileup on each split_bams folder
    # Due to this, we need to append the min_depth used to the output bcf folder
    # Since this is just submitting slurm jobs, we don't need to worry about
    # how many cores we use
    parallel::mclapply(min_depth,
                       mc.cores = 100,
                       mc.preschedule = FALSE,
                       function(this_min_depth) {

        replace_tibble_snp <-
            dplyr::tribble(
                ~find,                    ~replace,
                "placeholder_account",    account,
                "placeholder_slurm_out",  paste0(
                                                temp_dir, "/",
                                                slurm_base,
                                                "_mpileup-%j.out"
                                          ),
                "placeholder_array_max",  as.character(array_max),
                "placeholder_bam_dir",    bam_out_dir,
                "placeholder_ref_fasta",  ref_fasta,
                "placeholder_ploidy",     ploidy,
                "placeholder_bcf_dir",    paste0(bcf_dir, "_c", this_min_depth),
                "placeholder_min_depth",  as.character(this_min_depth)
            )

        # Call mpileup on each split_bams folder using a template and substituting
        # out the placeholder fields and index the individual bcf files
        result <-
            use_sbatch_template(replace_tibble_snp,
                                "snp_call_mpileup_template.job",
                                warning_label = "Calling SNPs",
                                submit = submit,
                                file_dir = temp_dir,
                                temp_prefix = paste0(sbatch_base, "mpileup_"))
        })
    # Delete contents of the split_bams folder
    # I should change how I do this, this may be a bit dangerous
    if (cleanup) {
        system(paste0(
            "rm ",
            bam_out_dir,
            "*.bam"
        ))
    }
    return(0)
}

#' Transform the ploidy argument into a valid argument for bcftools
#'
#' @param ploidy Either a file path to a ploidy file, or a string indicating
#'  the ploidy.
#'
#' @return A string that can be passed to use_sbatch_template() to fill in
#'  placeholder_ploidy
#'
#' @details GRCh37 is hg19, GRCh38 is hg38, X, Y, 1, mm10_hg19 is our mixed
#' species reference with species prefixes on chromosomes, mm10_hg38 is our
#' mixed reference from 10x, mm10 is mm10
#'
#' @noRd
pick_ploidy <- function(ploidy) {
    if (file.exists(ploidy)) {
        return(paste("--ploidy-file", ploidy))
    } else if (file.exists(paste0(find.package("rrrSnvs"),
                                  "/extdata/",
                                  ploidy,
                                  "_ploidy.txt"))) {
        return(paste0("--ploidy-file ",
                      find.package("rrrSnvs"),
                      "/extdata/",
                      ploidy,
                      "_ploidy.txt"))
    } else if (ploidy %in% c("GRCh37", "GRCh38", "X", "Y", "1")) {
        return(paste("--ploidy", ploidy))
    } else {
        warning("Ploidy argument not valid. Did you mean to pass a file path or ",
                "spell something wrong?",
                immediate. = TRUE)
        stop()
    }
}

#' Merge bcfs generated by call_snps() and write out a distance matrix
#'
#' @param bcf_in_dir The directory containing the bcfs to merge
#' @param out_bcf The path to the output bcf file
#' @param submit Whether to submit the job to slurm
#' @param account The cluster account to use in the slurm script
#' @param slurm_base The name of the slurm output file
#' @param sbatch_base The prefix to use with the sbatch job file
#' @param temp_dir The directory to write temporary files to
#' @param cleanup Whether to delete the bcfs after merging
#'
#' @return 0 if successful
#'
#' @noRd
merge_bcfs <- function(bcf_in_dir,
                       out_bcf,
                       submit = TRUE,
                       account = "gdrobertslab",
                       slurm_base = "slurmOut",
                       sbatch_base = "sbatch_",
                       temp_dir,
                       cleanup = TRUE) {

    # use template to merge bcfs and write out a distance matrix, substituting
    # out the placeholder fields
    replace_tibble_merge <-
        dplyr::tribble(
            ~find,                      ~replace,
            "placeholder_account",      account,
            "placeholder_slurm_out",    paste0(
                                            temp_dir, "/",
                                            slurm_base,
                                            "_merge-%j.out"
                                        ),
            "placeholder_bcf_out",      out_bcf,
            "placeholder_bcf_dir",      bcf_in_dir
        )

    # Call mpileup on each split_bams folder using a template and substituting
    # out the placeholder fields and index the individual bcf files
    result <-
        use_sbatch_template(replace_tibble_merge,
                            "snp_call_merge_template.job",
                            warning_label = "Calling SNPs",
                            submit = submit,
                            file_dir = temp_dir,
                            temp_prefix = paste0(sbatch_base, "merge_"))

    # remove individual bcf files
    if (cleanup) {
        unlink(bcf_in_dir, recursive = TRUE)
    }
    return(0)
}

#' Use Phylogenetic Tree from Merged BCF File to Find Similar Clusters
#'
#' This function calculates how samples (from pseudobulked clusters) from a
#' merged BCF file cluster hierarchically by collapsing a bootstrapped
#' hierarchical clustered tree. This provides mathematical testing of
#' genotypic divergence between the samples (clusters).
#'
#' @param merged_bcf_file Character. Path to the merged BCF file.
#' @param output_file Character. Path to the output file for assigned groupings
#'   per sample (cluster). Default is merged_bcf_file + "_dist.txt".
#' @param min_snvs_per_cluster Numeric. Minimum number of SNVs per cluster.
#'   Default is 500.
#' @param max_prop_missing_at_site Numeric. Maximum proportion of missing data
#'   allowed at a variant site. Default is 0.9.
#' @param n_bootstraps Numeric. Number of bootstrap iterations. Default is 1000.
#' @param bootstrap_cutoff Numeric. Bootstrap cutoff threshold to use for
#'   collapsing groups of samples (clusters). Default is 0.99.
#' @param tree_figure_file Character. Path to the output tree figure.
#' @param verbose Logical. If TRUE, enables verbose output. Default is TRUE.
#' @param sbatch_base Character. Base name for SLURM batch scripts. Can put
#'   these in /tmp/ if you want. Default is "./sbatch_dist".
#' @param account Character. SLURM account name. Default is "gdrobertslab".
#' @param slurm_base Character. Path to the SLURM output file.
#' @param temp_dir The directory to write temporary files to.
#' @param submit Logical. If TRUE, submits the SLURM job. Default is TRUE. Used
#'   mostly for testing.
#'
#' @return The result of the SLURM job submission.
#'
#' @examples
#' \dontrun{
#'   group_clusters_by_dist(
#'       merged_bcf_file = "path/to/merged.bcf",
#'       tree_figure_file = "path/to/tree_figure.png"
#'   )
#' }
#'
#' @export
group_clusters_by_dist <- function(
    merged_bcf_file,
    output_file = paste0(merged_bcf_file, "_dist.txt"),
    min_snvs_per_cluster = 500,
    max_prop_missing_at_site = 0.9,
    n_bootstraps = 1000,
    bootstrap_cutoff = 0.99,
    tree_figure_file,
    verbose = TRUE,
    sbatch_base = "sbatch_dist",
    account = "gdrobertslab",
    slurm_base = "slurmOut",
    temp_dir,
    submit = TRUE) {
    py_file <-
        paste0(find.package("rrrSnvs"),
               "/exec/vcfToMatrix.py")

    verbose_setting <-
        dplyr::if_else(verbose, "--verbose", "")

    replace_tibble_dist <-
        dplyr::tribble(
            ~find,                              ~replace,
            "placeholder_account",              account,
            "placeholder_slurm_out",            paste0(
                                                    temp_dir, "/",
                                                    slurm_base,
                                                    "dist_-%j.out"
                                                ),
            "placeholder_py_script",            py_file,
            "placeholder_bcf_input",            merged_bcf_file,
            "placeholder_min_snvs",             as.character(min_snvs_per_cluster),
            "placeholder_max_missing",          as.character(max_prop_missing_at_site),
            "placeholder_n_bootstrap",          as.character(n_bootstraps),
            "placeholder_bootstrap_threshold",  as.character(bootstrap_cutoff),
            "placeholder_fig_file",             tree_figure_file,
            "placeholder_verbose",              verbose_setting,
            "placeholder_groups_output",        output_file
        )

    # Call mpileup on each split_bams folder using a template and substituting
    # out the placeholder fields and index the individual bcf files
    result <-
        use_sbatch_template(replace_tibble_dist,
                            "vcf_to_matrix.sh",
                            warning_label = "Calculating distance from BCF",
                            submit = submit,
                            file_dir = temp_dir,
                            temp_prefix = paste0(sbatch_base, "_dist"))

    return(result)
}

#' Use the metadata in a Seurat object to determine which clusters are
#' predominantly control cell types
#'
#' @param sobject A Seurat object
#' @param normal_celltypes A vector of cell types that are considered non-tumor
#' @param cluster_col The name of the column that was used to define clusters in
#'  get_snp_tree()
#' @param celltype_col The name of the column that has cell type information
#' @param min_prop_control The minimum proportion of cells in a cluster that
#'  are control cell types to consider that cluster a control cluster
#'
#' @return A vector of cluster names that are predominantly control cell types
#'
#' @export
match_celltype_clusters <- function(sobject,
                                    normal_celltypes,
                                    cluster_col,
                                    celltype_col,
                                    min_prop_control = 0.5) {

    control_clusters <-
        sobject@meta.data %>%
        dplyr::select(dplyr::all_of(c(cluster_col, celltype_col))) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(cluster_col))) %>%
        dplyr::filter(sum(get(celltype_col) %in% normal_celltypes) /
                      dplyr::n() > min_prop_control) %>%
        dplyr::pull(cluster_col) %>%
        unique()

    return(control_clusters)
}
