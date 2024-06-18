test_that("check_cellid_bam_table stops if missing columns", {
    test_cbt <-
        tibble::tibble(cell_barcode = c("ATGG-1", "ATGC-1"),
                       bam_file = c("file1.bam", "file2.bam"))
    expect_error(check_cellid_bam_table(test_cbt),
                 "Columns should include cell_id, cell_group and bam_file")
})

test_that("Table has correct column names", {
    correct_table <-
        data.frame(
            cell_barcode = c("ATGG-1", "ATGC-1", "ATGA-1"),
            cell_group = c("Group1", "Group2", "Group3"),
            bam_file = c("file1.bam", "file2.bam", "file3.bam")
        )
    expect_equal(check_cellid_bam_table(correct_table), 0)
})

# Test for missing columns
test_that("Error on missing columns", {
    # Missing cell_group
    missing_columns_table <-
        data.frame(
            cell_barcode = c("ATGG-1", "ATGC-1", "ATGA-1"),
            bam_file = c("file1.bam", "file2.bam", "file3.bam")
        )
  expect_error(check_cellid_bam_table(missing_columns_table))
})

# Test for minimum number of groups
test_that("Error when less than three unique cell_groups", {
    few_groups_table <-
        data.frame(
            cell_barcode = c("ATGG-1", "ATGC-1", "ATGA-1"),
            cell_group = c("Group1", "Group2", "Group1"),
            bam_file = c("file1.bam", "file2.bam", "file3.bam")
        )

  expect_error(check_cellid_bam_table(few_groups_table))
})

# Test for empty table
test_that("Error on empty table", {
    empty_table <-
        data.frame(
            cell_barcode = character(),
            cell_group = character(),
            bam_file = character()
        )
    expect_error(check_cellid_bam_table(empty_table))
})

# Test for cell barcode format
test_that("Warning for incorrect cell barcode format", {
    incorrect_format_table <-
        data.frame(
            cell_barcode = c("bob_ATGG-1", "bob_ATGC-1", "bob_ATGA-1"),
            cell_group = c("Group1", "Group2", "Group3"),
            bam_file = c("file1.bam", "file2.bam", "file3.bam")
        )

    expect_warning(check_cellid_bam_table(incorrect_format_table))
})

# Test for unique cell_group to bam_file mapping
test_that("Error when cell_group is not unique to a single bam_file", {
    non_unique_group_table <-
        data.frame(
            cell_barcode = c("ATGG-1", "ATGC-1", "ATGA-1", "ATGC-1"),
            cell_group = c("Group1", "Group2", "Group3", "Group2"),
            bam_file = c("file1.bam", "file1.bam", "file1.bam", "file2.bam")
        )

    expect_error(check_cellid_bam_table(non_unique_group_table))
})

# Test for spaces in cell_group or bam_file
test_that("Error when cell_group or bam_file contains spaces", {
    spaces_table <-
        data.frame(
            cell_barcode = c("ATGG-1", "ATGC-1", "ATGA-1"),
            cell_group = c("Group 1", "Group2", "Group3"),
            bam_file = c("file1.bam", "file 2.bam", "file3.bam")
        )
    expect_error(check_cellid_bam_table(spaces_table))
})

# Test for cell_group starting with a letter
test_that("Error when cell_group does not start with a letter", {
    non_letter_start_table <-
        data.frame(
            cell_barcode = c("ATGG-1", "ATGC-1", "ATGA-1"),
            cell_group = c("1Group1", "2Group2", "3Group3"),
            bam_file = c("file1.bam", "file2.bam", "file3.bam")
        )

  expect_error(check_cellid_bam_table(non_letter_start_table))
})
