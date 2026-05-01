# Convenience wrapper: orchestrates Stage 1 → 2 → 3 in a single call.

#' Convert a Giotto object to a deployable Vitessce GitHub Pages site
#'
#' Orchestrates the full GiottoConnect pipeline:
#' 1. Exports the Giotto object to a SpatialData-Zarr store
#'    ([export_giotto_to_zarr()]).
#' 2. Generates a Vitessce view config JSON
#'    ([create_vitessce_config()]).
#' 3. Scaffolds a Vite/React GitHub Pages repository
#'    ([scaffold_vitessce_site()]).
#'
#' After the function completes, follow the printed instructions to push the
#' output directory to GitHub and enable GitHub Pages.
#'
#' @param gobject A Giotto object (`giotto` class from GiottoClass).
#' @param output_dir Character. Parent directory where the repository folder
#'   and intermediate zarr will be written.
#' @param repo_name Character. Name of the GitHub repository (also the name
#'   of the output folder inside `output_dir`).
#' @param github_username Character. Your GitHub username.
#' @param dataset_name Character. Name for the dataset (used as the zarr
#'   store filename and in the Vitessce config). Default: `"dataset"`.
#' @param spat_unit Character. Spatial unit to export. Default: `"cell"`.
#' @param feat_type Character. Feature type to export. Default: `"rna"`.
#' @param region Character. Name of the SpatialData shapes/points element
#'   that holds cell locations (must match what `giottoToSpatialData()`
#'   produces). Default: `"cell"`.
#' @param obs_set_paths Character vector or `NULL`. Paths to categorical obs
#'   columns for the observation-set panel (e.g. cluster labels). Auto-
#'   detected when `NULL`. Default: `NULL`.
#' @param obs_set_names Character vector or `NULL`. Display names for the
#'   observation sets. Auto-inferred when `NULL`. Default: `NULL`.
#' @param data_mode Character. `"github"` to bundle zarr data inside the
#'   repo (served via GitHub Pages); `"remote"` to point to an external
#'   data host. Default: `"github"`.
#' @param data_base_url Character or `NULL`. Required when
#'   `data_mode = "remote"`. Full base URL to the hosted zarr data, ending
#'   with `/`. Default: `NULL`.
#' @param schema_version Character. Vitessce config schema version.
#'   Default: `"1.0.16"`.
#' @param overwrite Logical. Overwrite existing zarr or repo if present.
#'   Default: `FALSE`.
#' @param python_env Character. Conda environment name.
#'   Default: `"giotto-connect"`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return Invisibly returns the path to the scaffolded repository directory.
#' @export
#'
#' @examples
#' \dontrun{
#' library(GiottoData)
#' gobject <- GiottoData::loadGiottoMini("merfish")
#'
#' setup_giotto_connect_env()
#'
#' result_dir <- giotto_to_vitessce_site(
#'     gobject         = gobject,
#'     output_dir      = "~/vitessce-output",
#'     repo_name       = "merfish-viewer",
#'     github_username = "myuser",
#'     dataset_name    = "mini_merfish"
#' )
#' }
giotto_to_vitessce_site <- function(
    gobject,
    output_dir,
    repo_name,
    github_username,
    dataset_name    = "dataset",
    spat_unit       = "cell",
    feat_type       = "rna",
    region          = "cell",
    obs_set_paths   = NULL,
    obs_set_names   = NULL,
    data_mode       = c("github", "remote"),
    data_base_url   = NULL,
    schema_version  = "1.0.16",
    overwrite       = FALSE,
    python_env      = "giotto-connect",
    verbose         = TRUE
) {
    data_mode <- match.arg(data_mode)

    # -- Determine base_url for the Vitessce config ---------------------------
    if (data_mode == "github") {
        base_url <- paste0(
            "https://", github_username, ".github.io/", repo_name, "/data/"
        )
    } else {
        if (is.null(data_base_url)) {
            stop("[GiottoConnect] data_base_url must be provided when data_mode = 'remote'.")
        }
        base_url <- data_base_url
    }

    # -- Intermediate paths ---------------------------------------------------
    zarr_staging <- file.path(output_dir, "_zarr_temp")
    if (!dir.exists(zarr_staging)) dir.create(zarr_staging, recursive = TRUE)

    zarr_path   <- file.path(zarr_staging, paste0(dataset_name, ".zarr"))
    config_path <- file.path(zarr_staging, "vitessce_config.json")

    # -- Stage 1: Export Giotto → SpatialData-Zarr ----------------------------
    if (verbose) message("[GiottoConnect] === Stage 1/3: Exporting to SpatialData-Zarr ===")

    export_giotto_to_zarr(
        gobject      = gobject,
        output_dir   = zarr_staging,
        dataset_name = dataset_name,
        spat_unit    = spat_unit,
        feat_type    = feat_type,
        overwrite    = overwrite,
        python_env   = python_env,
        verbose      = verbose
    )

    # -- Stage 2: Create Vitessce config --------------------------------------
    if (verbose) message("[GiottoConnect] === Stage 2/3: Creating Vitessce view config ===")

    create_vitessce_config(
        zarr_path      = zarr_path,
        dataset_name   = dataset_name,
        base_url       = base_url,
        schema_version = schema_version,
        region         = region,
        obs_set_paths  = obs_set_paths,
        obs_set_names  = obs_set_names,
        output_file    = config_path,
        python_env     = python_env,
        verbose        = verbose
    )

    # -- Stage 3: Scaffold GitHub Pages site ----------------------------------
    if (verbose) message("[GiottoConnect] === Stage 3/3: Scaffolding site ===")

    scaffold_vitessce_site(
        config_json_path = config_path,
        zarr_dir         = zarr_path,
        output_dir       = output_dir,
        repo_name        = repo_name,
        github_username  = github_username,
        dataset_name     = dataset_name,
        data_mode        = data_mode,
        data_base_url    = data_base_url,
        overwrite        = overwrite,
        verbose          = verbose
    )

    repo_dir <- file.path(output_dir, repo_name)
    invisible(repo_dir)
}
