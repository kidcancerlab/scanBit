#' Add SNV Group to Seurat Object
#'
#' This function adds a new column to the metadata of a Seurat object based on
#' a provided SNV (single nucleotide variant) group file. The new column will
#' contain the SNV group information for each cluster in the Seurat object.
#'
#' @param sobject A Seurat object to which the SNV group information will be added.
#' @param snv_group_file A file path to a tab-separated values (TSV) file
#'   containing the SNV group information. The file should have two columns:
#'   "cluster" and "dist_group". This is output from `get_snp_tree()`.
#' @param new_column A string specifying the name of the new column to be added
#'   to the Seurat object's metadata.
#' @param cell_group A string specifying the name of the column in the Seurat
#'   object's metadata that contains the cluster information.This is the same
#'   column used in `get_snp_tree()` in the table provided to the
#'   cellid_bam_table argument.
#'
#' @return The Seurat object with the new column added to its metadata.
#'
#' @examples
#' \dontrun{
#' sobj <- add_snv_group_to_sobj(
#'     sobject = seurat_object,
#'     snv_group_file = "path/to/snv_group_file_5_dist.txt",
#'     new_column = "snv_group_5"
#' )
#' }
#'
#' @export
add_snv_group_to_sobj <- function(sobject,
                                  snv_group_file,
                                  new_column,
                                  cell_group) {
    cluster_group_key <-
        readr::read_tsv(
            file = snv_group_file,
            col_names = c("cluster", "dist_group"),
            col_types = c("cf")
        ) %>%
        dplyr::pull("dist_group", name = cluster)

    if (!cell_group %in% colnames(sobject@meta.data)) {
        stop("cell_group column: '",
             cell_group,
             "' not found in sobject metadata")
    }

    sobject[[new_column]] <-
        cluster_group_key[sobject@meta.data[[cell_group]]] %>%
        as.vector()

    return(sobject)
}

#' Label Tumor Groups in a Seurat Object Using SNV Groups
#'
#' This function labels tumor cells in a Seurat object based on the specified
#' single nucleotide variant (SNV) group column and normal clusters.
#'
#' @param sobject A Seurat object containing single-cell RNA-seq data.
#' @param cell_group A character string specifying the column name in the
#'   Seurat object metadata that contains the cluster information. Default is
#'   "used_clusters".
#' @param snv_group_col A character string specifying the column name in the
#'   Seurat object metadata that contains the SNV group information derived from
#'   `add_snv_group_to_sobj()`.
#' @param normal_clusters A vector of cluster identifiers that are considered
#'   normal (non-tumor) clusters. These clusters should be present in the
#'   `cell_group` column of the Seurat object metadata.
#' @param tumor_call_column A character string specifying the column name in the
#'   Seurat object metadata where the tumor call labels will be stored. Default
#'   is "snv_tumor_call".
#'
#' @return A Seurat object with an additional column in the metadata indicating
#'   whether each cell belongs to a normal, tumor, or unknown group. The
#'   "unknown" cells indicate there was no snv tumor call for those cells.
#'
#' @details If tumor_call_column already exists in the Seurat object metadata,
#'   it will be overwritten. The function identifies normal SNV groups
#'   based on the provided normal clusters and labels cells accordingly. If no
#'   normal clusters are found, the function will return an error message.
#'
#' @examples
#' \dontrun{
#'   sobject <- label_tumor_cells(
#'       sobject = seurat_obj,
#'       cell_group = "used_clusters",
#'       snv_group_col = "snv_group_20",
#'       normal_clusters = c("cluster1", "cluster2"),
#'       tumor_call_column = "snv_tumor_call"
#'   )
#' }
#'
#' @export
label_tumor_cells <- function(sobject,
                              cell_group = "used_clusters",
                              snv_group_col,
                              normal_clusters,
                              tumor_call_column = "snv_tumor_call") {
    if (tumor_call_column %in% colnames(sobject@meta.data)) {
        message("tumor_call_column already exists in sobject, overwriting")
    }

    if (!is.null(normal_clusters)) {
        norm_snv_groups <-
            table(
                sobject@meta.data[[cell_group]],
                sobject@meta.data[[snv_group_col]]
            ) %>%
            as.data.frame() %>%
            dplyr::rename("groups" = "Var1", "snv_group" = "Var2") %>%
            dplyr::filter(Freq > 0) %>%
            dplyr::group_by(dplyr::pick("snv_group")) %>%
            dplyr::group_split() %>%
            lapply(function(x) {
                dplyr::if_else(
                    any(x$groups %in% normal_clusters),
                    x$snv_group[1],
                    NA)}) %>%
            unlist() %>%
            stats::na.omit() %>%
            as.character()

        if (length(norm_snv_groups) > 0) {
            # if tree_group contains the control, then rename it to "normal"
            # and label the other group(s) as "tumor"
            sobject@meta.data[[tumor_call_column]] <-
                dplyr::if_else(
                    sobject@meta.data[[snv_group_col]] %in% norm_snv_groups,
                    "normal",
                    dplyr::if_else(
                        is.na(sobject@meta.data[[snv_group_col]]),
                        "unknown",
                        "tumor"
                    )
                )
        } else {
            stop("No normal clusters found, check the normal_clusters argument",
                 " and make sure you are specifying the right groups.")
        }
    }

    return(sobject)
}
