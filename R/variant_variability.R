#' Use Merged BCF File to Quantify Variant Variability Across Groups
#'
#' This function ...............
#'
#' @inheritParams get_snp_tree
#' @param merged_bcf_file Character. Path to the merged BCF file.
#' @param output_file Character. Path to the output file cross group . Default
#'   is merged_bcf_file + "_variable.tsv".
#' @param min_snvs_per_cluster Numeric. Minimum number of SNVs per cluster.
#'   Default is 250.
#'
#' @return The result of the SLURM job submission.
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
        output_file = paste0(merged_bcf_file, "_variable.tsv"),
        min_snvs_per_cluster = 250,
        max_prop_missing_at_site = 0.5,
        verbose = TRUE,
        job_base = "batch",
        account = "gdrobertslab",
        log_base,
        slurm_base = "slurmOut",
        temp_dir = tempdir(),
        submit = TRUE,
        other_job_header_options = "",
        other_batch_options = "",
        sge_q,
        sge_parallel_environment,
        job_scheduler = "slurm") {
    py_file <- file.path(find.package("scanBit"), "exec/id_diff_vars.py")

    verbose_setting <-
        dplyr::if_else(verbose, "--verbose", "")

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
            "placeholder_min_snvs",         as.character(
                                                min_snvs_per_cluster
                                            ),
            "placeholder_max_missing",      as.character(
                                                max_prop_missing_at_site
                                            ),
            "placeholder_verbose",          verbose_setting,
            "placeholder_groups_output",    output_file
        )

    # Call mpileup on each cell id file using a template and substituting
    # out the placeholder fields and index the individual bcf files
    result <-
        use_sbatch_template(
            replace_tibble_dist,
            "vcf_to_matrix.sh",
            warning_label = "Calculating variable variants from BCF",
            submit = submit,
            file_dir = temp_dir,
            temp_prefix = paste0(job_base, "_var")
        )

    return(result)
}
