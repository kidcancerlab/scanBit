#' Download Demo Data
#'
#' Downloads and extracts the scanBit testing-data into the installed package
#' directory. Data can be found here:
#' https://doi.org/10.6084/m9.figshare.33399454
#'
#' @param data_loc Directory in which to download and extract the demo data.
#'   Defaults to the installed scanBit package directory.
#'
#' @noRd
download_demo_data <- function(data_loc = find.package("scanBit")) {
  zip_file <- file.path(data_loc, "scanBitTestingData.zip")

  download.file(
    url = "https://ndownloader.figshare.com/files/68075737",
    destfile = zip_file
  )

  unzip(zip_file, exdir = data_loc)

  file.remove(zip_file)
}

#' Prepare Demo Data
#'
#' Downloads the demo data and updates the cell barcode table so its BAM file
#' paths are relative to the installed scanBit package directory. Data are
#' downloaded from https://doi.org/10.6084/m9.figshare.33399454.
#'
#' @export
prep_demo_data <- function() {
  download_demo_data()
  package_location <- find.package("scanBit")

  c_b_t_file <- file.path(package_location, "cell_barcode_table.qs2")

  qs2::qs_read(c_b_t_file) |>
    dplyr::mutate(bam_file = paste0(package_location, "/", bam_file)) |>
    qs2::qs_save(c_b_t_file)
}
