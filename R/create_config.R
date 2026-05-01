# Stage 2: Generate a Vitessce view config from a SpatialData-Zarr store
#
# Uses the vitessce Python package (via reticulate) to build the config
# programmatically, then serializes it to JSON via jsonlite.
#
# SpatialDataWrapper is used (not AnnDataWrapper) because Stage 1 writes a
# full SpatialData zarr store.

#' Create a Vitessce view config for a SpatialData-Zarr store
#'
#' Builds a Vitessce JSON configuration that points to the provided zarr
#' store and includes three coordinated views: a UMAP scatterplot, a spatial
#' view (centroids / polygon shapes), and an observation-set panel for
#' cluster labels.
#'
#' @param zarr_path Character. Path to the local `.zarr` store produced by
#'   [export_giotto_to_zarr()]. Used to infer the zarr filename and, when
#'   `obs_set_paths = NULL`, to auto-detect categorical obs columns via the
#'   `.zmetadata` file.
#' @param dataset_name Character. Human-readable name for the dataset shown
#'   inside Vitessce.
#' @param base_url Character. Base URL where the zarr data will be served,
#'   e.g. `"https://<user>.github.io/<repo>/data/"`. The zarr filename is
#'   appended automatically.
#' @param schema_version Character. Vitessce config schema version.
#'   Default: `"1.0.16"`.
#' @param region Character. Name of the SpatialData element (shapes or
#'   points) that holds cell locations. Must match the element name used by
#'   `giottoToSpatialData()`. Default: `"cell"`.
#' @param obs_type Character. The observation type label used in Vitessce.
#'   Default: `"cell"`.
#' @param obs_embedding_paths Character vector. Paths to embedding arrays in
#'   the zarr store. Default: `c("obsm/X_umap")`.
#' @param obs_embedding_names Character vector. Display names for the
#'   embeddings. Must be the same length as `obs_embedding_paths`.
#'   Default: `c("UMAP")`.
#' @param obs_set_paths Character vector or `NULL`. Paths to categorical obs
#'   columns to expose as observation sets (e.g. cluster labels). If `NULL`,
#'   paths are auto-detected from the zarr `.zmetadata`. Default: `NULL`.
#' @param obs_set_names Character vector or `NULL`. Display names for the
#'   observation sets. If `NULL`, the last path component is used.
#'   Default: `NULL`.
#' @param output_file Character. Full path where the JSON config will be
#'   written. Default: `<zarr_dir>/vitessce_config.json`.
#' @param python_env Character. Name of the conda environment to use.
#'   Default: `"giotto-connect"`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return Invisibly returns the path to the written JSON config file.
#' @export
#'
#' @examples
#' \dontrun{
#' config_path <- create_vitessce_config(
#'     zarr_path    = "/tmp/mini_merfish.zarr",
#'     dataset_name = "mini_merfish",
#'     base_url     = "https://myuser.github.io/my-repo/data/"
#' )
#' }
create_vitessce_config <- function(
    zarr_path,
    dataset_name,
    base_url,
    schema_version       = "1.0.16",
    region               = "cell",
    obs_type             = "cell",
    obs_embedding_paths  = c("obsm/X_umap"),
    obs_embedding_names  = c("UMAP"),
    obs_set_paths        = NULL,
    obs_set_names        = NULL,
    output_file          = file.path(dirname(zarr_path), "vitessce_config.json"),
    python_env           = "giotto-connect",
    verbose              = TRUE
) {
    # -- Input validation -----------------------------------------------------
    checkmate::assert_string(zarr_path)
    checkmate::assert_string(dataset_name, min.chars = 1L)
    checkmate::assert_string(base_url)
    checkmate::assert_string(schema_version)
    checkmate::assert_string(region)
    checkmate::assert_string(obs_type)
    checkmate::assert_character(obs_embedding_paths, min.len = 1L)
    checkmate::assert_character(obs_embedding_names, len = length(obs_embedding_paths))
    checkmate::assert_flag(verbose)

    if (!dir.exists(zarr_path)) {
        stop("[GiottoConnect] zarr_path does not exist: ", zarr_path)
    }

    # -- Auto-detect obs_set_paths from .zmetadata ----------------------------
    if (is.null(obs_set_paths)) {
        obs_set_paths <- .detect_obs_set_paths(zarr_path, verbose = verbose)
    }
    if (is.null(obs_set_names)) {
        obs_set_names <- basename(obs_set_paths)
    }

    # -- Ensure Python env is ready -------------------------------------------
    .ensure_python()

    # -- Build zarr URL -------------------------------------------------------
    zarr_filename <- basename(zarr_path)
    # Ensure base_url ends with /
    if (!grepl("/$", base_url)) base_url <- paste0(base_url, "/")
    zarr_url <- paste0(base_url, zarr_filename)

    if (verbose) {
        message("[GiottoConnect] Stage 2: Building Vitessce config ...")
        message("[GiottoConnect]   zarr URL: ", zarr_url)
        message("[GiottoConnect]   obs_set_paths: ", paste(obs_set_paths, collapse = ", "))
    }

    # -- Build config via vitessce Python -------------------------------------
    vt <- .gc_python_env$vitessce

    wrapper <- vt$SpatialDataWrapper(
        sdata_url            = zarr_url,
        obs_embedding_paths  = reticulate::r_to_py(as.list(obs_embedding_paths)),
        obs_embedding_names  = reticulate::r_to_py(as.list(obs_embedding_names)),
        obs_set_paths        = reticulate::r_to_py(as.list(obs_set_paths)),
        obs_set_names        = reticulate::r_to_py(as.list(obs_set_names)),
        region               = region,
        obs_type             = obs_type
    )

    vc <- vt$VitessceConfig(
        schema_version = schema_version,
        name           = dataset_name
    )

    ds <- vc$add_dataset(name = dataset_name)
    ds$add_object(wrapper)

    v_scatter <- vc$add_view(vt$ViewType$SCATTERPLOT, dataset = ds)
    v_scatter$use_coordination(
        vc$add_coordination_by_dict(
            reticulate::r_to_py(list(embeddingType = obs_embedding_names[[1]]))
        )[[1]]
    )
    v_spatial <- vc$add_view(vt$ViewType$SPATIAL,   dataset = ds)
    v_sets    <- vc$add_view(vt$ViewType$OBS_SETS,  dataset = ds)

    vc$layout(
        vt$hconcat(
            v_spatial,
            vt$vconcat(v_scatter, v_sets)
        )
    )

    # -- Serialize to JSON ----------------------------------------------------
    if (verbose) message("[GiottoConnect]   Serializing config to JSON ...")

    config_json <- .vitessce_config_to_json(vc, vt)

    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    writeLines(config_json, output_file)

    if (verbose) message("[GiottoConnect] Stage 2 complete: ", output_file)

    invisible(output_file)
}

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' Auto-detect categorical obs set paths from zarr .zmetadata
#' @noRd
.detect_obs_set_paths <- function(zarr_path, verbose = TRUE) {
    zmeta_path <- file.path(zarr_path, ".zmetadata")
    if (!file.exists(zmeta_path)) {
        warning(
            "[GiottoConnect] .zmetadata not found at ", zmeta_path,
            ". Cannot auto-detect obs_set_paths. Please provide them manually."
        )
        return(character(0))
    }

    zmeta <- jsonlite::fromJSON(zmeta_path, simplifyVector = FALSE)
    all_keys <- names(zmeta$metadata)

    # Look for obs/<column_name>/.zarray entries (each obs column has one)
    obs_arrays <- all_keys[grepl("^obs/[^/]+/\\.zarray$", all_keys)]
    col_names  <- sub("^obs/([^/]+)/\\.zarray$", "\\1", obs_arrays)

    # Exclude internal zarr columns (typically start with . or are index-like)
    col_names <- col_names[!grepl("^\\.", col_names)]
    col_names <- col_names[col_names != "_index"]

    paths <- paste0("obs/", col_names)

    if (verbose && length(paths) > 0) {
        message("[GiottoConnect]   Auto-detected obs set paths: ",
                paste(paths, collapse = ", "))
    }

    paths
}

#' Serialize a vitessce VitessceConfig object to a JSON string
#'
#' Tries vc$to_dict() + py_to_r() + jsonlite first. Falls back to vc$export()
#' if the py_to_r() conversion produces an unexpected structure.
#' @noRd
.vitessce_config_to_json <- function(vc, vt) {
    # Primary path: vc$to_dict() -> py_to_r() -> jsonlite
    config_dict <- tryCatch({
        raw <- vc$to_dict()
        reticulate::py_to_r(raw)
    }, error = function(e) NULL)

    if (!is.null(config_dict) && is.list(config_dict)) {
        return(jsonlite::toJSON(config_dict, auto_unbox = TRUE, pretty = TRUE,
                                null = "null", na = "null"))
    }

    # Fallback: export to temp dir, read JSON, return as string
    tmp_dir <- tempfile("vt_export_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    vc$export(to = "files", base_url = ".", out_dir = tmp_dir)

    json_files <- list.files(tmp_dir, pattern = "\\.json$", full.names = TRUE)
    if (length(json_files) == 0L) {
        stop("[GiottoConnect] Vitessce config export produced no JSON files.")
    }

    # The primary config file is typically named after the schema or "config.json"
    cfg_file <- json_files[grepl("config", basename(json_files))]
    if (length(cfg_file) == 0L) cfg_file <- json_files[[1L]]

    paste(readLines(cfg_file, warn = FALSE), collapse = "\n")
}
