#' confirm_conda_env
#'
#' @details Checks if the conda environment "scanBit_xkcd_1337" exists and
#' creates it if not. Also confirms that the version of the environment matches
#' the package version.
#'
#' @param general Boolean indicating if the conda environment should be
#'   created using specific versions (FALSE) or just using general modules
#'   (TRUE).
#'
#' @return TRUE if the environment is correct or if it was created correctly,
#'   otherwise FALSE.
#'
#' @export
confirm_conda_env <- function(general = FALSE) {
    # Eventually, maybe drop this in the package package

    env_found <- conda_env_exists()

    conda_version_right <- right_conda_version()

    if (env_found && conda_version_right) {
        return(TRUE)
    } else if (!env_found) {
        message("Creating required conda environment scanBit_xkcd_1337")
        env_made <- make_conda_env(general)

        return(env_made)
    } else if (!conda_version_right) {
        system(paste0(
            "conda remove -n scanBit_xkcd_1337 --all -y"
        ))
        env_made <- make_conda_env(general)

        return(env_made)
    }

    # Shouldn't be possible to get here
    return(FALSE)
}

#' Check if the Conda Environment Exists
#'
#' @details This function checks whether the Conda environment named
#' "scanBit_xkcd_1337" exists.
#' It first verifies that the `conda` command is available in the system's PATH.
#' If the `conda` command is not found, it stops with an error message.
#' If the `conda` command is found, it checks for the existence of the specified
#' Conda environment.
#'
#' @return A logical value indicating whether the Conda environment
#'   "scanBit_xkcd_1337" exists.
#'
#' @noRd
conda_env_exists <- function() {
    if (check_cmd("conda") != 0) {
        stop(
            "conda command not found, do you need to load a module or add ",
            "this to your PATH?"
        )
    }

    return_value <- system(
        "conda env list | grep scanBit_xkcd_1337",
        ignore.stdout = TRUE
    )

    return(return_value == 0)
}

#' Get the Conda Version for scanBit Package
#'
#' @details This function retrieves the Conda version for the scanBit package
#' by querying the environment configuration variables of the
#' "scanBit_xkcd_1337" Conda environment. It then compares the retrieved Conda
#' version with the package version of scanBit and returns a logical value
#' indicating if they're equal.
#'
#' @return A logical value indicating whether the Conda version matches the
#'   package version of scanBit.
#'
#' @noRd
right_conda_version <- function() {
    if (conda_env_exists()) {
        package_version <- utils::packageVersion("scanBit")
        conda_version <-
            system(
                "conda env config vars list -n scanBit_xkcd_1337",
                intern = TRUE
            ) |>
            stringr::str_remove("scanBit_version = ")
        return(conda_version == package_version)
    }
    return(FALSE)
}

#' Create a Conda Environment
#'
#' This function creates a Conda environment based on the provided general
#' configuration or parameters.
#'
#' @param general A list or object containing the general configuration settings
#'   required to create the Conda environment.
#'
#' @return Returns TRUE if env was created
#'
#' @noRd
make_conda_env <- function(general) {
    if (general) {
        conda_yml_file <- paste0(find.package("scanBit"), "/conda_general.yml")
    } else {
        conda_yml_file <- paste0(find.package("scanBit"), "/conda.yml")
    }

    return_value <-
        system(paste0(
            "conda env create -n scanBit_xkcd_1337 --file ",
            conda_yml_file
        ))

    if (return_value != 0 && !general) {
        user_choice <- readline(
            paste(
                "It looks like you tried to make the conda environment",
                "with specific version requirements and it failed.",
                "Would you like to try making it with less stringent",
                "version requirements? please enter yes or no\n"
            )
        )
        if (user_choice == "yes") {
            return_value <-
                system(paste0(
                    "conda env create -n scanBit_xkcd_1337 --file ",
                    paste0(find.package("scanBit"), "/conda_general.yml")
                ))
        }
    }

    if (return_value == 0) {
        return(TRUE)
    }

    return(FALSE)
}
