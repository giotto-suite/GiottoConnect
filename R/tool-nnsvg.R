# nnSVG (Weber et al. 2023) connection. nnSVG runs its own R code (a
# nearest-neighbor Gaussian process per gene, fitted with BRISC); this file
# builds its input from a giotto object (normalized expression and spatial
# locations, optionally covariates from the cell metadata) and stores its
# per-gene statistics as feature metadata (adr/0008 in GiottoClass).

#' @title nnSVG parameters
#' @name nnsvgParam
#' @description
#' Parameters for [GiottoClass::analyzeData()] that run nnSVG (Weber et al.
#' 2023), which finds spatially variable genes by fitting a nearest-neighbor
#' Gaussian process to each gene and testing it against a model without
#' spatial terms (likelihood ratio test). nnSVG is a Bioconductor package,
#' install it with `BiocManager::install("nnSVG")`.
#' **If you use it, please cite the original paper** (see References).
#' @param n_neighbors integer, nearest neighbors of the Gaussian process
#'   (nnSVG's default 10; 5 is faster with similar results)
#' @param order ordering of the spots for the nearest-neighbor
#'   approximation: `"AMMD"` (nnSVG's default) or `"Sum_coords"` (for fewer
#'   than 65 spots)
#' @param n_threads integer, genes fitted in parallel
#' @param seed integer seed set before nnSVG runs (its fits draw random
#'   numbers, so results vary slightly between runs without one); `NULL`
#'   leaves the random number generator alone. The caller's random state is
#'   restored afterwards
#' @returns A `nnsvgParam` object.
#' @references Weber LM, Saha A, Datta A, Hansen KD, Hicks SC. nnSVG for the
#'   scalable identification of spatially variable genes using
#'   nearest-neighbor Gaussian processes. *Nature Communications* 14, 4059
#'   (2023). \doi{10.1038/s41467-023-39748-z}
#' @seealso [runNnSVG()]
#' @examples
#' p <- nnsvgParam(n_neighbors = 5)
#' p$n_neighbors
#' @md
#' @export
nnsvgParam <- function(n_neighbors = 10L, order = c("AMMD", "Sum_coords"),
    n_threads = 1L, seed = 1L) {
    checkmate::assert_count(n_neighbors, positive = TRUE)
    order <- match.arg(order)
    checkmate::assert_count(n_threads, positive = TRUE)
    checkmate::assert_int(seed, null.ok = TRUE)
    p <- methods::new("nnsvgParam")
    p$n_neighbors <- as.integer(n_neighbors)
    p$order <- order
    p$n_threads <- as.integer(n_threads)
    p$seed <- if (!is.null(seed)) as.integer(seed)
    p
}

#' @rdname nnsvgParam
#' @exportClass nnsvgParam
setClass("nnsvgParam", contains = "analyzeParam")

#' @rdname nnsvgParam
#' @param x normalized (log) expression matrix, features x cells
#' @param param a `nnsvgParam`
#' @param spatial_locs two-column matrix of coordinates, rows in the column
#'   order of `x`
#' @param X `NULL` or a covariate model matrix (cells x covariates, with an
#'   intercept column), e.g. `model.matrix(~ cell_type)`
#' @details `analyzeData(x, param, spatial_locs, X)` returns a `data.table`
#'   with one row per feature (`feat_ID`) and nnSVG's columns: the fitted
#'   `sigma.sq` (spatial variance), `tau.sq` (nonspatial variance), `phi`
#'   (length scale), `loglik`, `runtime`, `mean`, `var`, `spcov` (spatial
#'   coefficient of variation), `prop_sv` (proportion of spatial variance),
#'   `loglik_lm`, `LR_stat`, `rank` (by `LR_stat`), `pval` and `padj` (BH).
#'   `mean`, `var` and `spcov` are computed from `x`, as nnSVG does for the
#'   `logcounts` of a `SpatialExperiment`.
#' @export
setMethod("analyzeData", signature(x = "allMatrix", param = "nnsvgParam"),
    function(x, param, spatial_locs, X = NULL, ...) {
        package_check("nnSVG", repository = "Bioc")
        spatial_locs <- as.matrix(spatial_locs)[, 1:2, drop = FALSE]
        if (nrow(spatial_locs) != ncol(x)) {
            stop("[nnsvgParam] `spatial_locs` needs one row per column of `x`",
                call. = FALSE)
        }
        if (!is.null(X) && nrow(X) != ncol(x)) {
            stop("[nnsvgParam] `X` needs one row per column of `x`",
                call. = FALSE)
        }
        if (is.null(rownames(x))) rownames(x) <- paste0("feat_", seq_len(nrow(x)))
        if (!is.null(param$seed)) {
            genv <- globalenv()
            old <- genv[[".Random.seed"]]
            on.exit({
                if (is.null(old)) {
                    suppressWarnings(rm(".Random.seed", envir = genv))
                } else {
                    assign(".Random.seed", old, envir = genv)
                }
            }, add = TRUE)
            set.seed(param$seed)
        }
        res <- nnSVG::nnSVG(x, spatial_coords = spatial_locs, X = X,
            n_neighbors = param$n_neighbors, order = param$order,
            n_threads = param$n_threads)
        res <- .nnsvg_add_moments(res, x)
        dt <- data.table::as.data.table(res)
        dt <- cbind(data.table::data.table(feat_ID = rownames(res)), dt)
        dt
    }
)

# mean, var and spcov of each gene, as nnSVG computes them for logcounts
.nnsvg_add_moments <- function(res, x) {
    n <- ncol(x)
    m <- Matrix::rowMeans(x)
    v <- (Matrix::rowSums(x^2) - n * m^2) / (n - 1)
    res[, "mean"] <- m
    res[, "var"] <- v
    res[, "spcov"] <- sqrt(res[, "sigma.sq"]) / m
    res
}

#' @title Run nnSVG spatially variable genes
#' @name runNnSVG
#' @description
#' Runs nnSVG (Weber et al. 2023) on the normalized expression and spatial
#' locations of a giotto object and stores its per-gene statistics as
#' feature metadata, in columns named `<name>_<statistic>` (e.g.
#' `nnSVG_LR_stat`, `nnSVG_rank`, `nnSVG_padj`, `nnSVG_prop_sv`; see
#' [nnsvgParam()] for the statistics). Features that are not tested (not in
#' `feats_to_use`, or constant across cells) get `NA`; ranks and adjusted
#' p-values are within the tested features. nnSVG fits one model per gene,
#' so testing a subset of features (e.g. those expressed in enough cells) is
#' much faster. **If you use it, please cite the original paper** (see
#' [nnsvgParam()]).
#' @param gobject giotto object
#' @param spat_unit,feat_type spatial unit and feature type
#' @param expression_values normalized expression values to use
#' @param spat_loc_name name of the spatial locations
#' @param feats_to_use features to test: `NULL` (all), a feature metadata
#'   column (features with `"yes"` or `TRUE`, e.g. `"hvf"`), or feature IDs
#' @param covariates `NULL` or cell metadata columns used as covariates
#'   (e.g. a cell type column), passed to nnSVG as
#'   `model.matrix(~ covariates)`
#' @param param a [nnsvgParam()]; its defaults when `NULL`
#' @param name prefix of the feature metadata columns
#' @param verbose be verbose
#' @returns The giotto object with nnSVG's statistics in its feature
#'   metadata.
#' @inherit nnsvgParam references
#' @seealso [nnsvgParam()]
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- runNnSVG(g, feats_to_use = "hvf")
#' fm <- fDataDT(g)
#' head(fm[order(fm$nnSVG_rank), ])
#' }
#' @md
#' @export
runNnSVG <- function(gobject, spat_unit = NULL, feat_type = NULL,
    expression_values = "normalized", spat_loc_name = "raw",
    feats_to_use = NULL, covariates = NULL, param = NULL, name = "nnSVG",
    verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_character(feats_to_use, null.ok = TRUE)
    checkmate::assert_character(covariates, null.ok = TRUE)
    checkmate::assert_string(name)
    param <- param %null% nnsvgParam()
    checkmate::assert_class(param, "nnsvgParam")
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject, spat_unit = spat_unit,
        feat_type = feat_type)

    x <- getExpression(gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = expression_values, output = "matrix")
    fm <- fDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
    feats <- .nnsvg_feats(feats_to_use, fm, rownames(x))
    x <- x[feats, , drop = FALSE]
    n <- ncol(x)
    m <- Matrix::rowMeans(x)
    keep <- (Matrix::rowSums(x^2) - n * m^2) > 1e-12 * n
    if (any(!keep)) {
        vmsg(.v = verbose, "runNnSVG:", sum(!keep),
            "constant features are not tested")
        x <- x[keep, , drop = FALSE]
    }
    if (!nrow(x)) stop("[runNnSVG] no features to test", call. = FALSE)

    sl <- getSpatialLocations(gobject, spat_unit = spat_unit,
        name = spat_loc_name, output = "data.table")
    sl <- sl[match(colnames(x), sl$cell_ID), ]
    if (anyNA(sl$cell_ID)) {
        stop("[runNnSVG] spatial locations miss some cells", call. = FALSE)
    }
    locs <- as.matrix(sl[, c("sdimx", "sdimy"), with = FALSE])

    X <- NULL
    if (!is.null(covariates)) {
        meta <- pDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
        miss <- setdiff(covariates, names(meta))
        if (length(miss)) {
            stop("[runNnSVG] no cell metadata column(s) ",
                paste(miss, collapse = ", "), call. = FALSE)
        }
        cov <- as.data.frame(meta[match(colnames(x), meta$cell_ID),
            covariates, with = FALSE])
        X <- stats::model.matrix(~., data = cov)
    }

    vmsg(.v = verbose, "runNnSVG: running nnSVG",
        as.character(utils::packageVersion("nnSVG")), "on", nrow(x),
        "features and", ncol(x), "cells")
    res <- analyzeData(x, param, spatial_locs = locs, X = X)

    stats <- setdiff(names(res), c("feat_ID", "runtime"))
    out <- data.table::data.table(feat_ID = fm$feat_ID)
    idx <- match(out$feat_ID, res$feat_ID)
    for (s in stats) {
        col <- paste0(name, "_", gsub(".", "_", s, fixed = TRUE))
        data.table::set(out, j = col, value = res[[s]][idx])
    }
    gobject <- addFeatMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, new_metadata = out, by_column = TRUE,
        column_feat_ID = "feat_ID")
    update_giotto_params(gobject, description = "_nnSVG",
        attachments = list(param = param))
}

# features to test from `feats_to_use` (NULL, a feature metadata column or IDs)
.nnsvg_feats <- function(feats_to_use, fm, all_feats) {
    if (is.null(feats_to_use)) return(all_feats)
    if (length(feats_to_use) == 1L && feats_to_use %in% names(fm)) {
        v <- fm[[feats_to_use]]
        sel <- fm$feat_ID[v %in% c("yes", TRUE)]
    } else {
        sel <- feats_to_use
    }
    sel <- intersect(all_feats, sel)
    if (!length(sel)) {
        stop("[runNnSVG] `feats_to_use` selects no features", call. = FALSE)
    }
    sel
}
