#' Assemble Tree Plots
#'
#' This function assembles tree plots from PNG images for different depths.
#'
#' @param fig_folder A character string specifying the folder where the PNG
#'   images are stored.
#' @param sample_id A character string specifying the sample identifier.
#' @param depths A numeric vector specifying the depths for which the tree
#'   plots are to be generated.
#'
#' @return A patchwork object containing the assembled tree plots arranged in
#'   a single column.
#'
#' @examples
#' \dontrun{
#' fig_folder <- "path/to/figures"
#' sample_id <- "sample123"
#' depths <- c(10, 20, 30)
#' assembled_plots <- assemble_tree_plots(fig_folder, sample_id, depths)
#' }
#'
#' @export
assemble_tree_plots <- function(fig_folder,
                                sample_id,
                                depths) {
    plot_list <-
        lapply(
            depths,
            function(this_depth) {
                png_path <-
                    paste0(
                        fig_folder, "/",
                        sample_id, "_", this_depth, "_tree.png"
                    )
                if (file.exists(png_path)) {
                    plot_out <-
                        ggplot2::ggplot() +
                        ggpubr::background_image(png::readPNG(png_path)) +
                        ggplot2::coord_fixed() +
                        ggplot2::ggtitle(paste0("Depth: ", this_depth))
                } else {
                    plot_out <-
                        ggplot2::ggplot() +
                        ggplot2::ggtitle(paste0("Depth: ", this_depth)) +
                        ggplot2::theme_void() +
                        ggplot2::geom_text(
                            x = 0.5,
                            y = 0.5,
                            ggplot2::aes(label = "No plot available")
                        )
                }
                return(plot_out)
            }
        )
    assembled_plots <- patchwork::wrap_plots(plot_list, ncol = 1)

    return(assembled_plots)
}