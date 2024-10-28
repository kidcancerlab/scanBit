#' confirm_conda_env
#'
#' @details Checks if the conda environment "rrrSNVs_xkcd_1337" exists and
#' creates it if not.
#'
#' @return NULL
#'
#' @keyword internal
confirm_conda_env <- function() {
    pkg_location <- find.package("rrrSnvs")

    env_found <- conda_env_exists()

    right_env_version <- get_conda_version()



    if (env_found && right_env_version) {

        conda_yml_file <-
            paste0(find.package("rrrSnvs"),
                   "/conda.yml")

        message("Creating required conda environment rrrSNVs_xkcd_1337")
        system(paste0("conda env create -n rrrSNVs_xkcd_1337 --file ",
                      conda_yml_file))
    } else if (env_found && !right_env_version) {


    }
    return()
}

#' Check if the Conda Environment Exists
#'
#' @details This function checks whether the Conda environment named
#' "rrrSNVs_xkcd_1337" exists.
#' It first verifies that the `conda` command is available in the system's PATH.
#' If the `conda` command is not found, it stops with an error message.
#' If the `conda` command is found, it checks for the existence of the specified
#' Conda environment.
#'
#' @return A logical value indicating whether the Conda environment "rrrSNVs_xkcd_1337" exists.
#'
#' @keyword internal
conda_env_exists <- function() {
    if (check_cmd("conda") != 0) {
        stop(
            "conda command not found, do you need to load a module or add ",
            "this to your PATH?"
        )
    }

    return(
        system(
            "conda env list | grep rrrSNVs_xkcd_1337",
            ignore.stdout = TRUE
        ) == 0
    )
}

#' Get the Conda Version for rrrSnvs Package
#'
#' @details This function retrieves the Conda version for the rrrSnvs package
#' by querying the environment configuration variables of the
#' "rrrSNVs_xkcd_1337" Conda environment. It then compares the retrieved Conda
#' version with the package version of rrrSnvs and returns a logical value
#' indicating if they're equal.
#'
#' @return A logical value indicating whether the Conda version matches the package version of rrrSnvs.
#'
#' @keyword internal
get_conda_version <- function() {
    if (conda_env_exists()) {
        package_version <- packageVersion("rrrSnvs")
        conda_version <-
            system(
                "conda env config vars list -n rrrSNVs_xkcd_1337",
                intern = TRUE
            ) |>
            stringr::str_remove("rrrSnvs_version = ")
        return(conda_version == package_version)
    } else {
        stop("Conda environment rrrSNVs_xkcd_1337 not found")
    }
    stop("Shouldn't be possible to get this error message in get_conda_version")
}
