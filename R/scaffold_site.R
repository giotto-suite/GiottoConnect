# Stage 3: Scaffold the Vite/React GitHub Pages site
#
# Generates a ready-to-push repository containing a Vite + React app that
# embeds Vitessce. Template files live in inst/templates/ and use
# {{TOKEN}} placeholders replaced by gsub() — NOT glue/sprintf, because
# JSX uses {expression} syntax everywhere which would conflict.

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
        # Fallback for development (devtools::load_all())
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
    # file: template source name, dest: destination path relative to repo_dir
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
    gitignore_lines <- c("node_modules/", "dist/", ".DS_Store", "*.log")
    writeLines(gitignore_lines, file.path(repo_dir, ".gitignore"))

    # -- Print post-scaffold instructions -------------------------------------
    if (verbose) {
        .print_deploy_instructions(repo_dir, repo_name, github_username, gh_pages_url)
    }

    invisible(repo_dir)
}

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' Fill a template file by replacing {{TOKEN}} placeholders
#'
#' Uses gsub() with regex patterns. Does NOT use glue/sprintf because JSX
#' source files contain {expression} syntax throughout.
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

#' NULL-coalescing operator (internal)
#' @noRd
`%||%` <- function(x, y) if (!is.null(x) && nzchar(x)) x else y
