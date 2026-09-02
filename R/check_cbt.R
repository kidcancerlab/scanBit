#' Check that the cellid_bam_table input is of the proper format
#'
#' @param cellid_bam_table the cellid_bam_table to check
#'
#' @noRd
check_cellid_bam_table <- function(cellid_bam_table) {
  check_column_names(cellid_bam_table)

  # Check that there are more than one group in cell_group
  check_multiple_categories(cellid_bam_table$cell_group, "cell_group")

  # Check that there is data in the table
  check_data_in_table(cellid_bam_table)

  # check that cell barcode is of format `[ATGC]+-1`
  check_barcode_format(cellid_bam_table)

  # check that each cell_group is unique to a single bam_file
  # if the same group label is present in multiple bam_files, when we call
  # snps and merge the bcfs, we'll get duplicate column names and errors
  check_bam_cell_group_unique(cellid_bam_table)

  # check that there are no spaces in the bam_file
  check_no_spaces(cellid_bam_table$bam_file, "bam_file")

  # check that all cell_groups start with a letter
  check_starts_with_letter(cellid_bam_table$cell_group, "cell_group")

  # Check cell_groups are safe for filenames
  check_str_safe_filename(cellid_bam_table$cell_group, "cell_group")
}

#' Check column names in cellid_bam_table
#'
#' @param cellid_bam_table the cellid_bam_table to check
#'
#' @noRd
check_column_names <- function(cellid_bam_table) {
  if (
    any(
      !c("cell_barcode", "cell_group", "bam_file") %in%
        colnames(cellid_bam_table)
    )
  ) {
    stop("Columns should include cell_barcode, cell_group and bam_file")
  }
}

#' Check multiple categories in vector
#'
#' @param cell_group vector of group labels
#' @param label label for error message
#'
#' @noRd
check_multiple_categories <- function(cell_group, label) {
  if (length(unique(cell_group)) < 2) {
    stop("There should be more than one group in ", label)
  }
}

#' Check that table contains data
#'
#' @param cellid_bam_table the cellid_bam_table to check
#'
#' @noRd
check_data_in_table <- function(cellid_bam_table) {
  if (nrow(cellid_bam_table) == 0) {
    stop("No data in cellid_bam_table")
  }
}


#' Check cell barcode format
#'
#' @param cellid_bam_table the cellid_bam_table to check
#'
#' @noRd
check_barcode_format <- function(cellid_bam_table) {
  if (!all(stringr::str_detect(cellid_bam_table$cell_barcode, "^[ATGC]+-1$"))) {
    warning(
      "!!!\nCell barcodes should be of format `^[ATGC]+-1$`.\n",
      "We're going to try anyways, but this might not work.\n!!!"
    )
  }
}

#' Check that each cell_group maps to single BAM file
#'
#' @param cellid_bam_table the cellid_bam_table to check
#'
#' @noRd
check_bam_cell_group_unique <- function(cellid_bam_table) {
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
}

#' Check that text does not contain spaces
#'
#' @param text vector of text to check
#' @param label label for error message
#'
#' @noRd
check_no_spaces <- function(text, label) {
  if (any(stringr::str_detect(text, " "))) {
    stop(label, " should not contain spaces")
  }
}

#' Check that text starts with a letter
#'
#' @param text vector of text to check
#' @param label label for error message
#'
#' @noRd
check_starts_with_letter <- function(text, label) {
  if (any(!stringr::str_detect(text, "^[A-Za-z]"))) {
    stop(label, " should start with a letter")
  }
}

#' Check that text is safe for use in filenames
#'
#' @param text vector of text to check
#' @param label label for error message
#'
#' @noRd
check_str_safe_filename <- function(text, label) {
  if (any(!stringr::str_detect(text, "^[A-Za-z0-9_\\-\\.]+$"))) {
    stop(
      label,
      " column should only contain letters, numbers, underscores, hyphens and ",
      "periods"
    )
  }
}
