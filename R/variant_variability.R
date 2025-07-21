#' Use Merged BCF File to Quantify Variant Variability Across Groups
#'
#' This function quantifies the variability of variants across two groups of
#' clusters in a merged BCF file. It uses a Python script to calculate the
#' variant variability and returns the results in a specified output file.
#'
#' @inheritParams get_snp_tree
#' @param merged_bcf_file Character. Path to the merged BCF file.
#' @param group_1 Comma delimited string of multiple cluster names for
#'   group 1. No spaces.
#' @param group_2 Comma delimited string of multiple cluster names for
#'   group 2. No spaces.
#' @param output_file Character. Path to the output file cross group . Default
#'   is merged_bcf_file + "_variable.tsv".
#' @param min_snvs_per_cluster Numeric. Minimum number of SNVs per cluster.
#'   Default is 250.
#'
#' @details
#' The output file has the following columns:
#' - chr: Chromosome
#' - pos: Position of the variant
#' - alt_gt: Alternate genotype
#' - group_1_within_dist: Genetic distance within group 1
#' - group_2_within_dist: Genetic distance within group 2
#' - cross_group_dist: Genetic distance between groups
#' - cross_minus_within_dist: Difference between cross group distance and
#'  within group distances
#'
#' @return The result of the batch job submission.
#'
#' @examples
#' \dontrun{
#'   calc_variant_variability(
#'       merged_bcf_file = "path/to/merged.bcf",
#'       output_file = "variable_variants.tsv"
#'   )
#' }
#'
#' @export
calc_variant_variability <- function(
        merged_bcf_file,
        group_1,
        group_2,
        output_file = paste0(merged_bcf_file, "_variable.tsv"),
        min_snvs_per_cluster = 250,
        max_prop_missing_at_site = 0.5,
        verbose = TRUE,
        job_base = "batch",
        account = "gdrobertslab",
        log_base = "slurm_log_var",
        temp_dir = tempdir(),
        submit = TRUE,
        other_job_header_options = "",
        other_batch_options = "",
        sge_q = "all.q",
        sge_parallel_environment = "smp",
        job_scheduler = "slurm") {
    py_file <- file.path(find.package("scanBit"), "exec/id_diff_vars.py")

    verbose_setting <-
        dplyr::if_else(verbose, "--verbose", "")

    if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive = TRUE)
    }

    job_header_other <-
        make_header_other_string(other_job_header_options, job_scheduler)

    replace_tibble_dist <-
        dplyr::tribble(
            ~find,                          ~replace,
            "placeholder_account",          account,
            "placeholder_job_log",          paste0(
                                                temp_dir, "/",
                                                log_base,
                                                "var-",
                                                get_job_id_str(job_scheduler),
                                                ".out"
                                            ),
            "placeholder_sge_q",            sge_q,
            "placeholder_sge_thread",       sge_parallel_environment,
            "placeholder_job_header_other", job_header_other,
            "placeholder_batch_other",      paste0(
                                                other_batch_options,
                                                collapse = "\n"
                                            ),
            "placeholder_py_script",        py_file,
            "placeholder_bcf_input",        merged_bcf_file,
            "placeholder_group_1",          group_1,
            "placeholder_group_2",          group_2,
            "placeholder_min_snvs",         as.character(
                                                min_snvs_per_cluster
                                            ),
            "placeholder_max_missing",      as.character(
                                                max_prop_missing_at_site
                                            ),
            "placeholder_verbose",          verbose_setting,
            "placeholder_out_file",    output_file
        )

    # Call mpileup on each cell id file using a template and substituting
    # out the placeholder fields and index the individual bcf files
    result <-
        use_job_template(
            replace_tibble_dist,
            "bcf_variant_variability.sh",
            warning_label = "Calculating variable variants from BCF",
            submit = submit,
            file_dir = temp_dir,
            temp_prefix = paste0(job_base, "_var_")
        )

    return(result)
}

#' Call Variants on Each Cell From Select Variants
#'
#' This function calls variants on each cell from the selected variants.
#'
#' @inheritParams get_snp_tree
#' @param bam_to_use Character. Path to the BAM file to use for variant calling.
#' @param variable_variants_file Path to the location of the file produced by
#'   calc_variant_variability().
#' @param cell_barcodes A single string with all cell barcodes to look for in
#'   the bam file. Each barcode should be separated by a space.
#' @param output_file Character. Path to the output file of per-cell variant
#'   calls.
#' @param variant_diff_cutoff Numeric threshold to use to selct variants from
#'   variable_variants_file. Any variant with cross_minus_within_dist greater
#'   or equal to the cutoff will be used.
#'
#' @return The result of the batch job submission. The data is printed to the
#' file specified by output_file.
#'
#' @examples
#' \dontrun{
#'   call_variants_single_cell(
#'       bam_to_use = "path/to/merged.bcf",
#'       variable_variants_file,
#'       cell_barcodes,
#'       ref_fasta,
#'       ploidy,
#'       output_file = "variable_variants.tsv"
#'   )
#' }
#'
#' @export
call_variants_single_cell <- function(
        bam_to_use,
        variable_variants_file,
        cell_barcodes,
        output_file = "cell_variant_calls.tsv",
        ref_fasta,
        min_depth = 1,
        variant_diff_cutoff = 1,
        ploidy,
        job_base = "batch",
        account = "gdrobertslab",
        log_base = "slurm_log_var",
        temp_dir = tempfile(pattern = "tempdir"),
        submit = TRUE,
        other_job_header_options = "",
        other_batch_options = "",
        sge_q = "all.q",
        sge_parallel_environment = "smp",
        job_scheduler = "slurm") {

    # Check that temp_dir exists, and if not, create it
    if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive = TRUE)
    }
    message("Using temporary directory: ", temp_dir)

    job_header_other <-
        make_header_other_string(other_job_header_options, job_scheduler)

    replace_tibble_dist <-
        dplyr::tribble(
            ~find,                           ~replace,
            "placeholder_account",           account,
            "placeholder_job_log",           paste0(
                                                temp_dir, "/",
                                                log_base,
                                                "var-",
                                                get_job_id_str(job_scheduler),
                                                ".out"
                                             ),
            "placeholder_sge_q",             sge_q,
            "placeholder_sge_thread",        sge_parallel_environment,
            "placeholder_job_header_other",  job_header_other,
            "placeholder_batch_other",       paste0(
                                                other_batch_options,
                                                collapse = "\n"
                                             ),
            "placeholder_variant_cutoff",    as.character(variant_diff_cutoff),
            "placeholder_variable_variants", variable_variants_file,
            "placeholder_ploidy",            ploidy,
            "placeholder_temp_dir",          temp_dir,
            "placeholder_min_depth",         as.character(min_depth),
            "placeholder_bam_file",          bam_to_use,
            "placeholder_ref_fasta",         ref_fasta,
            "placeholder_cell_barcodes",     cell_barcodes,
            "placeholder_out_file",          output_file
        )

    # Call mpileup on each cell id file using a template and substituting
    # out the placeholder fields and index the individual bcf files
    result <-
        use_job_template(
            replace_tibble_dist,
            "cell_var_calls.sh",
            warning_label = "Calculating single cell variants",
            submit = submit,
            file_dir = temp_dir,
            temp_prefix = paste0(job_base, "_var_")
        )

    return(result)
}
