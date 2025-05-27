#' Check that the cellid_bam_table input is of the proper format
#'
#' @param cellid_bam_table the cellid_bam_table to check
#' @keywords internal
#'
#' @return 0 if the table is of the proper format, otherwise an error is thrown
check_cellid_bam_table <- function(cellid_bam_table) {
    # Columns should be cell_barcode, cell_group and bam_file
    if (any(!c("cell_barcode",
               "cell_group",
               "bam_file") %in% colnames(cellid_bam_table))) {
        stop("Columns should include cell_barcode, cell_group and bam_file")
    }

    # Check that there are more than one group in cell_group
    if (length(unique(cellid_bam_table$cell_group)) < 2) {
        stop("There should be more than one group in cell_group")
    }

    # Check that there is data in the table
    if (nrow(cellid_bam_table) == 0) {
        stop("No data in cellid_bam_table")
    }

    # check that cell barcode is of format `[ATGC]+-1`
    if (!all(stringr::str_detect(cellid_bam_table$cell_barcode,
                                 "^[ATGC]+-1$"))) {
        warning("!!!\nCell barcodes should be of format `^[ATGC]+-1$`.\n",
                "We're going to try anyways, but this might not work.\n!!!")
    }

    # check that each cell_group is unique to a single bam_file
    # if the same group label is present in multiple bam_files, when we call
    # snps and merge the bcfs, we'll get duplicate column names and errors
    bams_per_group <-
        cellid_bam_table %>%
        dplyr::select(-cell_barcode) %>%
        dplyr::distinct() %>%
        dplyr::group_by(cell_group) %>%
        dplyr::summarize(n_bams = dplyr::n()) %>%
        dplyr::pull(n_bams) %>%
        max()
    if (bams_per_group != 1) {
        stop("Each cell_group should be unique to a single bam file")
    }

    # check that there are no spaces in either the cell_group or bam_file
    if (any(stringr::str_detect(cellid_bam_table$cell_group, " ")) ||
        any(stringr::str_detect(cellid_bam_table$bam_file, " "))) {
        stop("Neither cell_group nor bam_file should not contain spaces")
    }

    # check that all cell_groups start with a letter
    if (!all(stringr::str_detect(cellid_bam_table$cell_group, "^[A-Za-z]"))) {
        stop("All cell_groups should start with a letter")
    }
    return(0)
}