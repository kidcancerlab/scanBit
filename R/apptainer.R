#' Offer to download the Apptainer container
#'
#' Prompts the user to download the scanBit Apptainer container when the
#' requested SIF file is unavailable.
#'
#' @param sif_file Character. Path to the required Apptainer SIF file.
#'
#' @return Invisibly returns the result of [download_container()] when the
#'   user enters `"yes"`; otherwise, an error is thrown.
#' @export
ask_to_download_container <- function(sif_file) {
  user_choice <- readline(
    paste(
      "\n\nYou don't have the Apptainer container",
      sif_file,
      "available. Would you like to download it? please enter yes or no\n"
    )
  ) |>
    tolower()

  if (user_choice == "yes") {
    download_container()
  } else {
    stop("Container not available and user chose not to download it.")
  }
}

#' Download the scanBit Apptainer container
#'
#' Verifies that Apptainer is available and pulls the latest scanBit container
#' into the installed package directory.
#'
#' @return The exit status returned by [base::system()].
#' @noRd
download_container <- function() {
  check_cmd("apptainer")

  system(
    paste0(
      "cd ",
      find.package("scanBit"),
      "; apptainer pull oras://ghcr.io/kidcancerlab/scanbit_container:latest"
    )
  )
}

#' Generate job-template values for Apptainer
#'
#' Creates the Apptainer command and conda-comment prefix used when populating
#' scheduler job templates.
#'
#' @param bound_items Character vector of files or directories to bind into
#'   the container.
#' @param env_vars Character vector of environment variable names to pass to
#'   the container.
#' @param use_apptainer Logical. Whether to generate an Apptainer command.
#'
#' @return A named list containing `apptainer_cmd` and `conda_prefix` strings.
#' @noRd
get_apptainer_placeholders <- function(bound_items, env_vars, use_apptainer) {
  return(list(
    "apptainer_cmd" = make_apptainer_command(
      bound_items,
      env_vars,
      use_apptainer
    ),
    "conda_prefix" = if (use_apptainer) "# " else "",
    "apptainer_cmd_end" = if (use_apptainer) "\"" else ""
  ))
}

#' Build an Apptainer execution command
#'
#' Constructs a multiline command that binds supplied paths, forwards selected
#' environment variables, and runs the scanBit Apptainer container.
#'
#' @param bound_items Character vector of files or directories to bind into
#'   the container.
#' @param env_vars Character vector of environment variable names to pass to
#'   the container.
#' @param use_apptainer Logical. Whether to generate an Apptainer command.
#'
#' @return A character string containing the command, or an empty string when
#'   `use_apptainer` is `FALSE`.
#' @noRd
make_apptainer_command <- function(bound_items, env_vars, use_apptainer) {
  if (use_apptainer) {
    cmd_prefix <- "apptainer exec --home $(pwd) \\\n"

    # Using realpath here so the user can pass relative paths and they will be
    # correctly resolved
    # Using -s so that symbolic links are returned with symlinks rather than
    # being resolved to their literal paths, which would cause a mismatch
    # between the host and container paths
    if (!missing(bound_items)) {
      bindings <- sprintf("  -B $(realpath -s %s) \\", bound_items) |>
        paste0("\n", collapse = "")
    } else {
      bindings <- ""
    }

    # needed syntax: --env MYVAR=A,MYVAR2=B
    if (!missing(env_vars)) {
      env_vars <- paste(
        "  --env ",
        paste0(env_vars, "=", paste0("${", env_vars, "}"), collapse = ","),
        " \\\n"
      )
    } else {
      env_vars <- ""
    }

    sif_file <- paste0(
      "  ",
      find.package("scanBit"),
      "/scanbit_container_latest.sif \\\n",
      "  bash -c \""
    )

    apptainer_cmd <- paste(cmd_prefix, bindings, env_vars, sif_file)
  } else {
    apptainer_cmd <- ""
  }

  return(apptainer_cmd)
}
