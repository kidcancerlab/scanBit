utils::globalVariables(c(".", "avg.exp.scaled", "features.plot", "id",
                         "pct.exp", "score", "x", "y", "p_val_adj",
                         "avg_log2FC", "gene", "Freq", "Var1", "cid",
                         "col.fill", "freq", "label", "label14", "label30",
                         "lt_14", "lt_30", "database", "from", "to", "pearson",
                         "test_ligand", "ligand_target_matrix", "ligands",
                         "ligands_bona_fide", "lr_network", "lr_network_strict",
                         "receptors", "receptors_bona_fide", "value",
                         "weighted_networks", "weighted_networks_lr", "Phase",
                         "Cluster", "Proportion", "weight", "cluster",
                         "sil_width", "sil_vals", "res_vals", "num_clusters",
                         "min_val", "max_val", "median_val", "sd_val",
                         "feature", "Sample_ID", "exp_type", "suffix", "fastqs",
                         "library_type", "link_folder", "tx_id", "tar_folder",
                         "CB", "bam_file", "cell_barcode", "cell_group",
                         "tree_group", "sample_1", "sample_2", "group_count",
                         "n_bams", "snp_dist"))

#' Use a job template to submit a job to the cluster
#'
#' @inheritParams get_snp_tree
#' @param replace_tibble A tibble with two columns, find and replace.
#'   The find column should contain the placeholder text to be replaced and the
#'   replace column should contain the text to replace it with.
#' @param template The name of the template file to use. This should be a file
#'   in the rrrSingleCellUtils/inst folder. The job_scheduler argument will be
#'   prefixed to the template name to find the file.
#' @param file_dir The directory to write the temporary job file to.
#' @param temp_ext The extension to use for the temporary job file.
#' @param temp_prefix The prefix to use for the temporary job file.
#' @param warning_label A string to use in the warning message if the job
#'   submission fails.
#'
#' @return 0 if the job submission was successful, otherwise an error is
#'  thrown
use_job_template <- function(replace_tibble,
                             template,
                             file_dir = tempdir(),
                             temp_ext = ".sh",
                             temp_prefix = "job_",
                             warning_label = "",
                             submit = TRUE,
                             job_scheduler = "slurm") {
    job_template <-
        readr::read_file(
            file.path(
                find.package("scanBit"),
                paste0(
                    job_scheduler,
                    "_",
                    template
                )
            )
        )

    # Replace placeholders with real data
    for (i in seq_len(nrow(replace_tibble))) {
        job_template <-
            stringr::str_replace_all(
                job_template,
                pattern = replace_tibble$find[i],
                replacement = replace_tibble$replace[i]
            )
    }

    temp_file <-
        tempfile(fileext = temp_ext,
                 tmpdir = file_dir,
                 pattern = temp_prefix)

    readr::write_file(job_template, file = temp_file)

    if (submit == TRUE) {
        return_val <- system(paste(get_submit_cmd(job_scheduler), temp_file))
    } else {
        return_val <- 0
    }

    if (return_val != 0) {
        stop(paste0(warning_label,
                    " job submission failed. Error code ",
                    return_val))
    }
    return(0)
}

get_submit_cmd <- function(job_scheduler) {
    job_submit_cmd <-
        switch(
            job_scheduler,
            "slurm" = "sbatch",
            "sge" = "qsub",
            "bash" = "bash",
            stop(paste("Unknown job scheduler:", job_scheduler))
        )
    return(job_submit_cmd)
}

#' Check that a command is available on the system
#'
#' @param cmd The command to check for
#' @return 0 if the command is available, otherwise an error is thrown
check_cmd <- function(cmd) {
    if (Sys.which(cmd) == "") {
        stop(paste(cmd,
                   "command not found, do you need to load a module or add",
                   "this to your PATH?"))
    }
    return(0)
}

#' Generate SBATCH Options String
#'
#' This function takes a string of additional SBATCH options and formats it
#' into a string suitable for inclusion in an SBATCH script. If the input
#' string is empty, it prepends "#SBATCH" to each option and joins them with
#' newline characters. If the input string is not empty, it returns an empty
#' string.
#'
#' @param other_sbatch_options A character string containing additional SBATCH
#'   options. If not empty, the function will format it with "#SBATCH" and
#'   newlines. If empty, the function will return an empty string.
#'
#' @return A formatted character string of SBATCH options or an empty string.
#'
#' @noRd
make_sbatch_other_string <- function(other_sbatch_options) {
    if (length(other_sbatch_options) == 1 && other_sbatch_options == "") {
        sbatch_string <- ""
    } else {
        sbatch_string <-
            paste("#SBATCH", other_sbatch_options, collapse = "\n")
    }

    return(sbatch_string)
}