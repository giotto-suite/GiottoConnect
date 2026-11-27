# Scrublet (Wolock, Lopez and Klein 2019) connection. Scrublet is a Python
# package, run in the R session through reticulate in an environment that
# reticulate::py_require() builds on first use; this file hands it the counts
# of a giotto object and stores its doublet scores and calls as cell metadata
# (adr/0008 in GiottoClass).

#' @title Scrublet parameters
#' @name scrubletParam
#' @description
#' Parameters for [GiottoClass::analyzeData()] that run Scrublet (Wolock,
#' Lopez and Klein 2019), which detects doublets by simulating them from pairs
#' of observed cells and scoring each cell by the share of simulated doublets
#' among its nearest neighbors. It is meant for single-cell data, not for spots
#' that hold several cells. Scrublet is a Python package: the first run
#' installs it with reticulate (reticulate >= 1.41, no conda).
#' **If you use it, please cite the original paper** (see References).
#' @param expected_doublet_rate expected fraction of doublets (0.06 is from the
#'   10x Chromium guide)
#' @param min_counts,min_cells Scrublet's feature filter before its PCA: a
#'   feature needs at least `min_counts` counts in at least `min_cells` cells
#' @param min_gene_variability_pctl percentile of feature variability above
#'   which features are used
#' @param n_prin_comps number of principal components
#' @param threshold `NULL` (Scrublet sets the doublet score threshold
#'   itself) or a score above which cells are called doublets. Scrublet
#'   cannot always set one, e.g. on spot data, and then calls no doublets.
#' @param seed integer random state of Scrublet
#' @returns A `scrubletParam` object.
#' @references Wolock SL, Lopez R, Klein AM. Scrublet: computational
#'   identification of cell doublets in single-cell transcriptomic data.
#'   *Cell Systems* 8, 281-291.e9 (2019). \doi{10.1016/j.cels.2018.11.005}
#' @seealso [runScrublet()]
#' @examples
#' p <- scrubletParam(expected_doublet_rate = 0.08)
#' p$expected_doublet_rate
#' @md
#' @export
scrubletParam <- function(expected_doublet_rate = 0.06, min_counts = 1L,
    min_cells = 1L, min_gene_variability_pctl = 85, n_prin_comps = 30L,
    threshold = NULL, seed = 1234L) {
    checkmate::assert_number(expected_doublet_rate, lower = 0, upper = 1)
    checkmate::assert_count(min_counts)
    checkmate::assert_count(min_cells)
    checkmate::assert_number(min_gene_variability_pctl, lower = 0, upper = 100)
    checkmate::assert_count(n_prin_comps, positive = TRUE)
    checkmate::assert_number(threshold, lower = 0, upper = 1, null.ok = TRUE)
    checkmate::assert_int(seed)
    p <- methods::new("scrubletParam")
    p$expected_doublet_rate <- expected_doublet_rate
    p$min_counts <- as.integer(min_counts)
    p$min_cells <- as.integer(min_cells)
    p$min_gene_variability_pctl <- min_gene_variability_pctl
    p$n_prin_comps <- as.integer(n_prin_comps)
    p$threshold <- threshold
    p$seed <- as.integer(seed)
    p
}

#' @rdname scrubletParam
#' @exportClass scrubletParam
setClass("scrubletParam", contains = "analyzeParam")

# the paper, shown once per session and listed by objCitations()
setMethod("methodCitation", signature("scrubletParam"), function(x, ...) {
    .cite_article(
        author = c("Samuel L. Wolock", "Romain Lopez", "Allon M. Klein"),
        title = paste("Scrublet: computational identification of cell",
            "doublets in single-cell transcriptomic data"),
        journal = "Cell Systems", year = 2019, volume = 8,
        pages = "281-291.e9", doi = "10.1016/j.cels.2018.11.005"
    )
})

#' @rdname scrubletParam
#' @param x count matrix, features x cells
#' @param param a `scrubletParam`
#' @param verbose show Scrublet's progress messages
#' @details `analyzeData(x, param)` returns a `data.table` with one row per
#'   cell: `cell_ID`, `doublet_score` and `doublet` (Scrublet's call, `NA`
#'   when there is no `threshold` and Scrublet could not set one).
#' @export
setMethod("analyzeData", signature(x = "allMatrix", param = "scrubletParam"),
    function(x, param, verbose = FALSE, ...) {
        if (is.null(colnames(x))) colnames(x) <- paste0("cell_", seq_len(ncol(x)))
        # Scrublet takes cells as rows; it runs in a worker R session with
        # its own Python environment (GiottoUtils::py_run())
        res <- GiottoUtils::py_run(.scrublet_worker,
            args = list(counts = Matrix::t(methods::as(x, "CsparseMatrix")),
                p = as.list(param@param)),
            requirements = "scrublet", name = "scrublet",
            verbose = isTRUE(verbose))
        out <- data.table::data.table(cell_ID = colnames(x),
            doublet_score = res$score, doublet = res$call)
        attr(out, "py_versions") <- attr(res, "py_versions")
        out
    }
)

# Runs in the py_run() worker: only its arguments and `::` calls are
# available there.
.scrublet_worker <- function(counts, p) {
    scr <- reticulate::import("scrublet")
    s <- scr$Scrublet(counts_matrix = reticulate::r_to_py(counts),
        expected_doublet_rate = p$expected_doublet_rate,
        random_state = p$seed)
    res <- s$scrub_doublets(min_counts = p$min_counts,
        min_cells = p$min_cells,
        min_gene_variability_pctl = p$min_gene_variability_pctl,
        n_prin_comps = p$n_prin_comps, verbose = FALSE)
    calls <- res[[2]]
    if (!is.null(p$threshold)) {
        calls <- s$call_doublets(threshold = p$threshold, verbose = FALSE)
    }
    score <- as.numeric(res[[1]])
    list(score = score,
        call = if (is.null(calls)) rep(NA, length(score)) else as.logical(calls))
}

#' @title Run Scrublet doublet detection
#' @name runScrublet
#' @description
#' Runs Scrublet (Wolock, Lopez and Klein 2019) on the counts of a giotto
#' object and stores each cell's doublet score and call as cell metadata, in
#' columns `<name>_score` and `<name>_doublet`. When Scrublet cannot set a
#' doublet score threshold and the param gives none, only the scores are
#' stored, with a warning. Scrublet is meant for
#' single-cell data, not for spots that hold several cells. **If you use it,
#' please cite the original paper** (see [scrubletParam()]).
#' @param gobject giotto object
#' @param spat_unit,feat_type spatial unit and feature type
#' @param expression_values expression values to use; Scrublet expects counts
#' @param param a [scrubletParam()]; its defaults when `NULL`
#' @param name prefix of the cell metadata columns
#' @param verbose be verbose
#' @returns The giotto object with Scrublet's results in its cell metadata.
#' @inherit scrubletParam references
#' @seealso [scrubletParam()]
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- runScrublet(g)
#' pDataDT(g)[, c("cell_ID", "scrublet_score", "scrublet_doublet")]
#' }
#' @md
#' @export
runScrublet <- function(gobject, spat_unit = NULL, feat_type = NULL,
    expression_values = "raw", param = NULL, name = "scrublet",
    verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_string(name)
    param <- param %null% scrubletParam()
    checkmate::assert_class(param, "scrubletParam")
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject, spat_unit = spat_unit,
        feat_type = feat_type)

    x <- getExpression(gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = expression_values, output = "matrix")
    vmsg(.v = verbose, "runScrublet: running Scrublet on", nrow(x),
        "features and", ncol(x), "cells")
    res <- analyzeData(x, param, verbose = verbose)

    out <- data.table::data.table(cell_ID = res$cell_ID)
    data.table::set(out, j = paste0(name, "_score"), value = res$doublet_score)
    if (all(is.na(res$doublet))) {
        warning(wrap_txt("[runScrublet] Scrublet could not set a doublet score",
            "threshold, so only scores are stored; give one with",
            "scrubletParam(threshold =)"), call. = FALSE)
    } else {
        data.table::set(out, j = paste0(name, "_doublet"),
            value = res$doublet)
    }
    gobject <- addCellMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, new_metadata = out, by_column = TRUE,
        column_cell_ID = "cell_ID")
    update_giotto_params(gobject, description = "_scrublet",
        attachments = list(param = param,
            py_versions = attr(res, "py_versions")))
}
