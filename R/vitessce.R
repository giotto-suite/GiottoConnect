# Vitessce deployment pipeline
#
# Three-stage pipeline: Giotto object → SpatialData-Zarr → Vitessce view
# config → Vite/React GitHub Pages site.
#
# Stage 1 (export_giotto_to_zarr):   Giotto → SpatialData-Zarr
# Stage 2 (create_vitessce_config):  zarr store → Vitessce JSON config
# Stage 3 (scaffold_vitessce_site):  config + zarr → GitHub Pages repo
# Wrapper (giotto_to_vitessce_site): orchestrates stages 1–3 in one call
#
# Template files for the Vite/React app live in inst/templates/ and use
# {{TOKEN}} placeholders replaced by gsub() — NOT glue/sprintf, because
# JSX uses {expression} syntax throughout which would conflict.
#
# TODO(giotto-upstream): giottoToSpatialData() is the canonical conversion
# path for Stage 1. Any missing slot exports (UMAP in obsm, polygon shapes,
# cell metadata) should be added in GiottoClass rather than worked around here.


# =============================================================================
# Stage 1: Export Giotto → SpatialData-Zarr
# =============================================================================

#' Export a Giotto object to SpatialData-Zarr
#'
#' Converts a Giotto object to a SpatialData store and writes it as a Zarr
#' archive that Vitessce can consume. The output directory will contain a
#' `<dataset_name>.zarr/` folder with a consolidated `.zmetadata` file.
#'
#' @param gobject A Giotto object (`giotto` class from GiottoClass).
#' @param output_dir Character. Path to the directory where the zarr store
#'   will be written.
#' @param dataset_name Character. Base name for the output zarr store
#'   (without `.zarr` extension). Default: `"dataset"`.
#' @param spat_unit Character. Spatial unit to export (e.g., `"cell"`).
#'   Passed to `giottoToSpatialData()`. Default: `"cell"`.
#' @param feat_type Character. Feature type to export (e.g., `"rna"`).
#'   Passed to `giottoToSpatialData()`. Default: `"rna"`.
#' @param include_expression Logical. If `TRUE`, include the expression
#'   matrix in the export. Not yet implemented in v1; reserved for future
#'   use. Default: `FALSE`.
#' @param include_images Logical. If `TRUE`, include spatial images (as
#'   OME-Zarr). Not yet implemented in v1; reserved for future use.
#'   Default: `FALSE`.
#' @param overwrite Logical. If `TRUE`, overwrite an existing zarr store at
#'   the target path. Default: `FALSE`.
#' @param python_env Character. Name of the conda environment to use.
#'   Default: `"giotto-connect"`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return Invisibly returns the path to the written zarr store.
#' @export
#'
#' @examples
#' \dontrun{
#' library(GiottoData)
#' gobject <- GiottoData::loadGiottoMini("merfish")
#' zarr_path <- export_giotto_to_zarr(
#'     gobject    = gobject,
#'     output_dir = tempdir(),
#'     dataset_name = "mini_merfish"
#' )
#' }
export_giotto_to_zarr <- function(
    gobject,
    output_dir,
    dataset_name       = "dataset",
    spat_unit          = "cell",
    feat_type          = "rna",
    include_expression = FALSE,
    include_images     = FALSE,
    overwrite          = FALSE,
    python_env         = "giotto-connect",
    verbose            = TRUE
) {
    # -- Input validation -----------------------------------------------------
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_string(output_dir)
    checkmate::assert_string(dataset_name, min.chars = 1L)
    checkmate::assert_string(spat_unit)
    checkmate::assert_string(feat_type)
    checkmate::assert_flag(include_expression)
    checkmate::assert_flag(include_images)
    checkmate::assert_flag(overwrite)
    checkmate::assert_flag(verbose)

    # -- Prepare output path --------------------------------------------------
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    zarr_path <- file.path(output_dir, paste0(dataset_name, ".zarr"))

    if (dir.exists(zarr_path)) {
        if (!overwrite) {
            stop(
                "[GiottoConnect] Zarr store already exists at: ", zarr_path,
                "\n  Use overwrite = TRUE to replace it."
            )
        }
        unlink(zarr_path, recursive = TRUE)
    }

    # -- Ensure Python env is ready -------------------------------------------
    .ensure_python()

    # -- Convert Giotto → SpatialData -----------------------------------------
    if (verbose) message("[GiottoConnect] Stage 1: Converting Giotto object to SpatialData ...")

    sdata <- GiottoClass::giottoToSpatialData(
        gobject   = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # -- Write SpatialData to zarr --------------------------------------------
    if (verbose) message("[GiottoConnect]   Writing zarr to: ", zarr_path)

    sdata$write(zarr_path)

    # -- Consolidate metadata (required by Vitessce) --------------------------
    if (verbose) message("[GiottoConnect]   Consolidating zarr metadata ...")

    zarr_mod <- .gc_python_env$zarr
    zarr_mod$consolidate_metadata(zarr_path)

    if (verbose) message("[GiottoConnect] Stage 1 complete: ", zarr_path)

    invisible(zarr_path)
}


# =============================================================================
# Stage 2: Generate Vitessce view config
# =============================================================================

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

    zmeta    <- jsonlite::fromJSON(zmeta_path, simplifyVector = FALSE)
    all_keys <- names(zmeta$metadata)

    obs_arrays <- all_keys[grepl("^obs/[^/]+/\\.zarray$", all_keys)]
    col_names  <- sub("^obs/([^/]+)/\\.zarray$", "\\1", obs_arrays)

    col_names <- col_names[!grepl("^\\.", col_names)]
    col_names <- col_names[col_names != "_index"]

    paths <- paste0("obs/", col_names)

    if (verbose && length(paths) > 0L) {
        message("[GiottoConnect]   Auto-detected obs set paths: ",
                paste(paths, collapse = ", "))
    }

    paths
}

#' Serialize a VitessceConfig object to a JSON string
#'
#' Tries vc$to_dict() + py_to_r() + jsonlite first. Falls back to vc$export()
#' if the py_to_r() conversion produces an unexpected structure.
#' @noRd
.vitessce_config_to_json <- function(vc, vt) {
    config_dict <- tryCatch({
        raw <- vc$to_dict()
        reticulate::py_to_r(raw)
    }, error = function(e) NULL)

    if (!is.null(config_dict) && is.list(config_dict)) {
        return(jsonlite::toJSON(config_dict, auto_unbox = TRUE, pretty = TRUE,
                                null = "null", na = "null"))
    }

    # Fallback: export to temp dir, read the JSON file
    tmp_dir <- tempfile("vt_export_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    vc$export(to = "files", base_url = ".", out_dir = tmp_dir)

    json_files <- list.files(tmp_dir, pattern = "\\.json$", full.names = TRUE)
    if (length(json_files) == 0L) {
        stop("[GiottoConnect] Vitessce config export produced no JSON files.")
    }

    cfg_file <- json_files[grepl("config", basename(json_files))]
    if (length(cfg_file) == 0L) cfg_file <- json_files[[1L]]

    paste(readLines(cfg_file, warn = FALSE), collapse = "\n")
}


# =============================================================================
# Stage 3: Scaffold the GitHub Pages site
# =============================================================================

#' Scaffold a Vite/React Vitessce viewer site
#'
#' Creates a directory structure ready to push to GitHub Pages. The output
#' contains a Vite + React 18 application that loads the provided Vitessce
#' config and, optionally, the zarr data.
#'
#' @param config_json_path Character. Path to the `vitessce_config.json`
#'   produced by [create_vitessce_config()].
#' @param zarr_dir Character. Path to the `.zarr` store produced by
#'   [export_giotto_to_zarr()]. Only used when `data_mode = "github"`.
#' @param output_dir Character. Parent directory where the repository folder
#'   will be created.
#' @param repo_name Character. Name of the GitHub repository (also the
#'   output folder name).
#' @param github_username Character. GitHub username, used to construct the
#'   GitHub Pages URL in `package.json` and instructions.
#' @param dataset_name Character or `NULL`. Dataset name used in template
#'   headers. If `NULL`, inferred from `basename(zarr_dir)` without
#'   `.zarr`. Default: `NULL`.
#' @param data_mode Character. One of `"github"` (zarr copied into the repo
#'   under `public/data/`) or `"remote"` (zarr hosted externally; nothing
#'   is copied). Default: `"github"`.
#' @param data_base_url Character or `NULL`. Required when
#'   `data_mode = "remote"`. The base URL where zarr data is hosted.
#'   Default: `NULL`.
#' @param overwrite Logical. If `TRUE`, overwrite an existing repo directory
#'   at the target path. Default: `FALSE`.
#' @param verbose Logical. Print progress messages and post-scaffold
#'   instructions. Default: `TRUE`.
#'
#' @return Invisibly returns the path to the scaffolded repository directory.
#' @export
#'
#' @examples
#' \dontrun{
#' scaffold_vitessce_site(
#'     config_json_path = "/tmp/vitessce_config.json",
#'     zarr_dir         = "/tmp/mini_merfish.zarr",
#'     output_dir       = "~/vitessce-sites",
#'     repo_name        = "merfish-viewer",
#'     github_username  = "myuser"
#' )
#' }
scaffold_vitessce_site <- function(
    config_json_path,
    zarr_dir,
    output_dir,
    repo_name,
    github_username,
    dataset_name  = NULL,
    data_mode     = c("github", "remote"),
    data_base_url = NULL,
    overwrite     = FALSE,
    verbose       = TRUE
) {
    # -- Input validation -----------------------------------------------------
    checkmate::assert_string(config_json_path)
    checkmate::assert_string(output_dir)
    checkmate::assert_string(repo_name, min.chars = 1L)
    checkmate::assert_string(github_username, min.chars = 1L)
    checkmate::assert_flag(overwrite)
    checkmate::assert_flag(verbose)

    data_mode <- match.arg(data_mode)

    if (!file.exists(config_json_path)) {
        stop("[GiottoConnect] config_json_path not found: ", config_json_path)
    }

    if (data_mode == "github") {
        checkmate::assert_string(zarr_dir)
        if (!dir.exists(zarr_dir)) {
            stop("[GiottoConnect] zarr_dir not found: ", zarr_dir)
        }
    }

    if (data_mode == "remote") {
        checkmate::assert_string(data_base_url, null.ok = FALSE)
    }

    if (is.null(dataset_name)) {
        dataset_name <- sub("\\.zarr$", "", basename(zarr_dir))
    }

    # -- Set up directory structure -------------------------------------------
    repo_dir    <- file.path(output_dir, repo_name)

    if (dir.exists(repo_dir)) {
        if (!overwrite) {
            stop(
                "[GiottoConnect] Repo directory already exists: ", repo_dir,
                "\n  Use overwrite = TRUE to replace it."
            )
        }
        unlink(repo_dir, recursive = TRUE)
    }

    src_dir     <- file.path(repo_dir, "src")
    pub_data    <- file.path(repo_dir, "public", "data")
    gh_workflow <- file.path(repo_dir, ".github", "workflows")

    for (d in c(src_dir, pub_data, gh_workflow)) {
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }

    if (verbose) message("[GiottoConnect] Stage 3: Scaffolding site at: ", repo_dir)

    # -- Token replacements ---------------------------------------------------
    gh_pages_url <- paste0("https://", github_username, ".github.io/", repo_name, "/")

    replacements <- list(
        REPO_NAME       = repo_name,
        GITHUB_USERNAME = github_username,
        DATASET_NAME    = dataset_name,
        GH_PAGES_URL    = gh_pages_url
    )

    # -- Locate templates -----------------------------------------------------
    tmpl_dir <- system.file("templates", package = "GiottoConnect")
    if (!nzchar(tmpl_dir) || !dir.exists(tmpl_dir)) {
        tmpl_dir <- file.path(
            dirname(dirname(sys.frame(1L)$ofile %||% "")),
            "inst", "templates"
        )
        if (!dir.exists(tmpl_dir)) {
            stop(
                "[GiottoConnect] Cannot find inst/templates/. ",
                "Install the package or run devtools::load_all() from the package root."
            )
        }
    }

    # -- Write parameterized template files -----------------------------------
    template_map <- list(
        list(file = "App.jsx",        dest = file.path(src_dir, "App.jsx")),
        list(file = "main.jsx",       dest = file.path(src_dir, "main.jsx")),
        list(file = "index.css",      dest = file.path(src_dir, "index.css")),
        list(file = "index.html",     dest = file.path(repo_dir, "index.html")),
        list(file = "vite.config.js", dest = file.path(repo_dir, "vite.config.js")),
        list(file = "package.json",   dest = file.path(repo_dir, "package.json")),
        list(file = "deploy.yml",     dest = file.path(gh_workflow, "deploy.yml"))
    )

    for (tm in template_map) {
        src_path <- file.path(tmpl_dir, tm$file)
        if (!file.exists(src_path)) {
            warning("[GiottoConnect] Template not found: ", src_path)
            next
        }
        filled <- .fill_template(src_path, replacements)
        writeLines(filled, tm$dest)
    }

    # -- Copy vitessce_config.json into src/ ----------------------------------
    file.copy(config_json_path, file.path(src_dir, "vitessce_config.json"))

    # -- Copy zarr data (github mode only) ------------------------------------
    if (data_mode == "github") {
        if (verbose) message("[GiottoConnect]   Copying zarr data to public/data/ ...")
        zarr_dest <- file.path(pub_data, basename(zarr_dir))
        .copy_dir(zarr_dir, zarr_dest)
    } else {
        if (verbose) {
            message("[GiottoConnect]   data_mode = 'remote': zarr NOT copied into repo.")
            message("[GiottoConnect]   Data will be served from: ", data_base_url)
        }
    }

    # -- Write .gitignore -----------------------------------------------------
    writeLines(c("node_modules/", "dist/", ".DS_Store", "*.log"),
               file.path(repo_dir, ".gitignore"))

    if (verbose) {
        .print_deploy_instructions(repo_dir, repo_name, github_username, gh_pages_url)
    }

    invisible(repo_dir)
}

#' Fill a template file by replacing {{TOKEN}} placeholders
#' @noRd
.fill_template <- function(template_path, replacements) {
    text <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
    for (nm in names(replacements)) {
        text <- gsub(
            pattern     = paste0("\\{\\{", nm, "\\}\\}"),
            replacement = replacements[[nm]],
            x           = text,
            fixed       = FALSE
        )
    }
    text
}

#' Recursively copy a directory
#' @noRd
.copy_dir <- function(src, dest) {
    if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)
    all_files <- list.files(src, recursive = TRUE, all.files = TRUE,
                            full.names = TRUE, include.dirs = FALSE)
    rel_paths  <- sub(paste0("^", src, "/?"), "", all_files)
    dest_paths <- file.path(dest, rel_paths)
    dest_dirs  <- unique(dirname(dest_paths))
    for (d in dest_dirs) {
        if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
    file.copy(all_files, dest_paths, overwrite = TRUE)
    invisible(NULL)
}

#' Print next-step deployment instructions to the console
#' @noRd
.print_deploy_instructions <- function(repo_dir, repo_name, github_username,
                                       gh_pages_url) {
    cat("\n")
    cat("[GiottoConnect] Site scaffolded at:", repo_dir, "\n\n")
    cat("Next steps to deploy to GitHub Pages:\n\n")
    cat("  1. cd", repo_dir, "\n")
    cat("  2. git init && git add . && git commit -m 'Initial Vitessce site'\n")
    cat("  3. Create a GitHub repo named '", repo_name, "':\n", sep = "")
    cat("       gh repo create ", github_username, "/", repo_name,
        " --public\n", sep = "")
    cat("  4. git remote add origin https://github.com/", github_username,
        "/", repo_name, ".git\n", sep = "")
    cat("  5. git push -u origin main\n")
    cat("  6. Enable GitHub Pages: repo Settings > Pages > Source: GitHub Actions\n")
    cat("  7. The GitHub Action will run npm install + build + deploy automatically.\n")
    cat("  8. Your site will be live at:", gh_pages_url, "\n\n")
}

#' @noRd
`%||%` <- function(x, y) if (!is.null(x) && nzchar(x)) x else y


# =============================================================================
# Wrapper: Giotto → Vitessce site (all 3 stages)
# =============================================================================

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

    invisible(file.path(output_dir, repo_name))
}
