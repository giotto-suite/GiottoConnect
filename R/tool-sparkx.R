# SPARK-X (Zhu, Sun and Zhou 2021) connection. SPARK-X runs its own R code (a
# covariance test between expression and spatial kernels); this file builds
# its input from a giotto object (counts and spatial locations, optionally
# covariates from the cell metadata) and stores its per-gene p-values as
# feature metadata (adr/0008 in GiottoClass).

#' @title SPARK-X parameters
#' @name sparkxParam
#' @description
#' Parameters for [GiottoClass::analyzeData()] that run SPARK-X (Zhu, Sun and
#' Zhou 2021), which finds spatially variable genes with a non-parametric
#' covariance test between each gene's counts and spatial kernels built from
#' the locations. It scales to hundreds of thousands of cells. SPARK-X is part
#' of the SPARK package on GitHub, install it with
#' `remotes::install_github("xzhoulab/SPARK")`.
#' **If you use it, please cite the original paper** (see References).
#' @param option `"mixture"` (SPARK-X's default: the projection kernel, five
#'   Gaussian and five cosine kernels, combined with the Cauchy combination
#'   rule) or `"single"` (projection kernel only)
#' @param n_cores integer, cores SPARK-X uses
#' @returns A `sparkxParam` object.
#' @references Zhu J, Sun S, Zhou X. SPARK-X: non-parametric modeling enables
#'   scalable and robust detection of spatial expression patterns for large
#'   spatial transcriptomic studies. *Genome Biology* 22, 184 (2021).
#'   \doi{10.1186/s13059-021-02404-0}
#' @seealso [runSparkX()]
#' @examples
#' p <- sparkxParam(option = "single")
#' p$option
#' @md
#' @export
sparkxParam <- function(option = c("mixture", "single"), n_cores = 1L) {
    option <- match.arg(option)
    checkmate::assert_count(n_cores, positive = TRUE)
    p <- methods::new("sparkxParam")
    p$option <- option
    p$n_cores <- as.integer(n_cores)
    p
}

#' @rdname sparkxParam
#' @exportClass sparkxParam
setClass("sparkxParam", contains = "analyzeParam")

# the paper, shown once per session and listed by objCitations()
setMethod("methodCitation", signature("sparkxParam"), function(x, ...) {
    .cite_article(
        author = c("Jiaqiang Zhu", "Shiquan Sun", "Xiang Zhou"),
        title = paste("SPARK-X: non-parametric modeling enables scalable and",
            "robust detection of spatial expression patterns for large",
            "spatial transcriptomic studies"),
        journal = "Genome Biology", year = 2021, volume = 22,
        pages = "184", doi = "10.1186/s13059-021-02404-0"
    )
})

#' @rdname sparkxParam
#' @param x count matrix, features x cells
#' @param param a `sparkxParam`
#' @param spatial_locs two-column matrix of coordinates, rows in the column
#'   order of `x`
#' @param X `NULL` or a covariate model matrix (cells x covariates, with an
#'   intercept column), e.g. `model.matrix(~ cell_type)`
#' @param verbose show SPARK-X's progress messages
#' @details `analyzeData(x, param, spatial_locs, X)` returns a `data.table`
#'   with one row per tested feature (`feat_ID`): `pval`, the combined p-value
#'   over the kernels, and `padj`, adjusted with Benjamini-Yekutieli, as
#'   SPARK-X reports them. SPARK-X leaves out cells and features without
#'   counts.
#' @export
setMethod("analyzeData", signature(x = "allMatrix", param = "sparkxParam"),
    function(x, param, spatial_locs, X = NULL, verbose = FALSE, ...) {
        package_check("SPARK", repository = "github:xzhoulab/SPARK")
        spatial_locs <- as.matrix(spatial_locs)[, 1:2, drop = FALSE]
        if (nrow(spatial_locs) != ncol(x)) {
            stop("[sparkxParam] `spatial_locs` needs one row per column of `x`",
                call. = FALSE)
        }
        if (!is.null(X) && nrow(X) != ncol(x)) {
            stop("[sparkxParam] `X` needs one row per column of `x`",
                call. = FALSE)
        }
        if (is.null(rownames(x))) rownames(x) <- paste0("feat_", seq_len(nrow(x)))
        if (!methods::is(x, "sparseMatrix")) {
            x <- methods::as(methods::as(as.matrix(x), "CsparseMatrix"),
                "generalMatrix")
        }
        run <- function() {
            SPARK::sparkx(x, spatial_locs, X_in = X, numCores = param$n_cores,
                option = param$option, verbose = isTRUE(verbose))
        }
        # SPARK-X prints its progress with cat()
        if (isTRUE(verbose)) {
            res <- run()
        } else {
            utils::capture.output(res <- run())
        }
        mt <- res$res_mtest
        data.table::data.table(feat_ID = rownames(mt),
            pval = mt$combinedPval, padj = mt$adjustedPval)
    }
)

#' @title Run SPARK-X spatially variable genes
#' @name runSparkX
#' @description
#' Runs SPARK-X (Zhu, Sun and Zhou 2021) on the counts and spatial locations
#' of a giotto object and stores its per-gene p-values as feature metadata,
#' in columns `<name>_pval` (combined over the spatial kernels) and
#' `<name>_padj` (Benjamini-Yekutieli). Features that are not tested (not in
#' `feats_to_use`, or without counts) get `NA`. **If you use it, please cite
#' the original paper** (see [sparkxParam()]).
#' @param gobject giotto object
#' @param spat_unit,feat_type spatial unit and feature type
#' @param expression_values expression values to use; SPARK-X expects counts
#' @param spat_loc_name name of the spatial locations
#' @param feats_to_use features to test: `NULL` (all), a feature metadata
#'   column (features with `"yes"` or `TRUE`, e.g. `"hvf"`), or feature IDs
#' @param covariates `NULL` or cell metadata columns used as covariates
#'   (e.g. a cell type column), passed to SPARK-X as
#'   `model.matrix(~ covariates)`
#' @param param a [sparkxParam()]; its defaults when `NULL`
#' @param name prefix of the feature metadata columns
#' @param verbose be verbose
#' @returns The giotto object with SPARK-X's p-values in its feature metadata.
#' @inherit sparkxParam references
#' @seealso [sparkxParam()]
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- runSparkX(g)
#' fm <- fDataDT(g)
#' head(fm[order(fm$sparkx_padj), ])
#' }
#' @md
#' @export
runSparkX <- function(gobject, spat_unit = NULL, feat_type = NULL,
    expression_values = "raw", spat_loc_name = "raw", feats_to_use = NULL,
    covariates = NULL, param = NULL, name = "sparkx", verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_character(feats_to_use, null.ok = TRUE)
    checkmate::assert_character(covariates, null.ok = TRUE)
    checkmate::assert_string(name)
    param <- param %null% sparkxParam()
    checkmate::assert_class(param, "sparkxParam")
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject, spat_unit = spat_unit,
        feat_type = feat_type)

    x <- getExpression(gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = expression_values, output = "matrix")
    fm <- fDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
    feats <- .tool_feats(feats_to_use, fm, rownames(x), fn = "runSparkX")
    x <- x[feats, , drop = FALSE]

    sl <- getSpatialLocations(gobject, spat_unit = spat_unit,
        name = spat_loc_name, output = "data.table")
    sl <- sl[match(colnames(x), sl$cell_ID), ]
    if (anyNA(sl$cell_ID)) {
        stop("[runSparkX] spatial locations miss some cells", call. = FALSE)
    }
    locs <- as.matrix(sl[, c("sdimx", "sdimy"), with = FALSE])

    X <- NULL
    if (!is.null(covariates)) {
        meta <- pDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
        miss <- setdiff(covariates, names(meta))
        if (length(miss)) {
            stop("[runSparkX] no cell metadata column(s) ",
                paste(miss, collapse = ", "), call. = FALSE)
        }
        cov <- as.data.frame(meta[match(colnames(x), meta$cell_ID),
            covariates, with = FALSE])
        X <- stats::model.matrix(~., data = cov)
    }

    vmsg(.v = verbose, "runSparkX: running SPARK-X",
        as.character(utils::packageVersion("SPARK")), "on", nrow(x),
        "features and", ncol(x), "cells")
    res <- analyzeData(x, param, spatial_locs = locs, X = X)

    out <- data.table::data.table(feat_ID = fm$feat_ID)
    idx <- match(out$feat_ID, res$feat_ID)
    data.table::set(out, j = paste0(name, "_pval"), value = res$pval[idx])
    data.table::set(out, j = paste0(name, "_padj"), value = res$padj[idx])
    gobject <- addFeatMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, new_metadata = out, by_column = TRUE,
        column_feat_ID = "feat_ID")
    update_giotto_params(gobject, description = "_sparkx",
        attachments = list(param = param))
}
