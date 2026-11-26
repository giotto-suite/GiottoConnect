# MAST (Finak et al. 2015) connection. MAST runs its own R code (a hurdle
# model per gene: a logistic model for detection and a Gaussian model for the
# expression level); this file builds its input from a giotto object (an
# expression matrix and cell groups, optionally covariates from the cell
# metadata) and returns MAST's per-gene tests as a table (adr/0008 in
# GiottoClass).

#' @title MAST parameters
#' @name mastParam
#' @description
#' Parameters for [GiottoClass::analyzeData()] that run MAST (Finak et al.
#' 2015), which tests each feature for differential expression between two
#' groups of cells with a hurdle model: a logistic regression for whether the
#' feature is detected and a linear regression for its level where it is.
#' The likelihood ratio test combines both parts. MAST is a Bioconductor
#' package, install it with `BiocManager::install("MAST")`.
#' **If you use it, please cite the original paper** (see References).
#' @param method fitting method passed to `MAST::zlm()`: `"bayesglm"` (MAST's
#'   default), `"glm"`, `"ridge"`, `"blmer"` or `"glmer"`
#' @param ebayes logical, shrink the variances with empirical Bayes (MAST's
#'   default `TRUE`)
#' @returns A `mastParam` object.
#' @references Finak G, McDavid A, Yajima M, Deng J, Gersuk V, Shalek AK,
#'   Slichter CK, Miller HW, McElrath MJ, Prlic M, Linsley PS, Gottardo R.
#'   MAST: a flexible statistical framework for assessing transcriptional
#'   changes and characterizing heterogeneity in single-cell RNA sequencing
#'   data. *Genome Biology* 16, 278 (2015). \doi{10.1186/s13059-015-0844-5}
#' @seealso [runMAST()]
#' @examples
#' p <- mastParam(method = "glm")
#' p$method
#' @md
#' @export
mastParam <- function(method = c("bayesglm", "glm", "ridge", "blmer", "glmer"),
    ebayes = TRUE) {
    method <- match.arg(method)
    checkmate::assert_flag(ebayes)
    p <- methods::new("mastParam")
    p$method <- method
    p$ebayes <- ebayes
    p
}

#' @rdname mastParam
#' @exportClass mastParam
setClass("mastParam", contains = "analyzeParam")

# the paper, shown once per session and listed by objCitations()
setMethod("methodCitation", signature("mastParam"), function(x, ...) {
    .cite_article(
        author = c("Greg Finak", "Andrew McDavid", "Masanao Yajima",
            "Jingyuan Deng", "Vivian Gersuk", "Alex K. Shalek",
            "Chloe K. Slichter", "Hannah W. Miller", "M. Juliana McElrath",
            "Martin Prlic", "Peter S. Linsley", "Raphael Gottardo"),
        title = paste("MAST: a flexible statistical framework for assessing",
            "transcriptional changes and characterizing heterogeneity in",
            "single-cell RNA sequencing data"),
        journal = "Genome Biology", year = 2015, volume = 16,
        pages = "278", doi = "10.1186/s13059-015-0844-5"
    )
})

#' @rdname mastParam
#' @param x expression matrix, features x cells (MAST expects log-scale
#'   values; it is made dense)
#' @param param a `mastParam`
#' @param groups group label of each cell. Named by cell ID (matched to
#'   `colnames(x)`) or in the column order of `x`. Cells labelled `NA` are
#'   left out.
#' @param test,ref the labels in `groups` to compare: `test` against `ref`
#' @param covariates `NULL` or a data.frame with one row per column of `x`
#'   and the covariates to adjust for (e.g. the detection rate)
#' @param verbose show MAST's progress messages
#' @details `analyzeData(x, param, groups, test, ref, covariates)` returns a
#'   `data.table` with one row per feature, sorted by `fdr`: `feats`,
#'   `Pr(>Chisq)` (the hurdle likelihood ratio test), `coef`, `ci.hi` and
#'   `ci.lo` (the log fold change of `test` over `ref` and its confidence
#'   interval), `fdr` (Benjamini-Hochberg) and `cluster` (`"<test>_vs_<ref>"`).
#' @export
setMethod("analyzeData", signature(x = "allMatrix", param = "mastParam"),
    function(x, param, groups, test, ref, covariates = NULL, verbose = FALSE,
        ...) {
        package_check("MAST", repository = "Bioc")
        checkmate::assert_string(test)
        checkmate::assert_string(ref)
        if (identical(test, ref)) {
            stop("[mastParam] `test` and `ref` must differ", call. = FALSE)
        }
        if (is.null(names(groups))) {
            if (length(groups) != ncol(x)) {
                stop("[mastParam] unnamed `groups` needs one label per ",
                    "column of `x`", call. = FALSE)
            }
            names(groups) <- colnames(x)
        }
        if (is.null(colnames(x))) colnames(x) <- names(groups)
        lab <- as.character(groups)[match(colnames(x), names(groups))]
        keep <- !is.na(lab) & lab %in% c(test, ref)
        for (g in c(test, ref)) {
            if (!any(lab[keep] == g)) {
                stop("[mastParam] no cells in group '", g, "'", call. = FALSE)
            }
        }
        if (!is.null(covariates) && nrow(covariates) != ncol(x)) {
            stop("[mastParam] `covariates` needs one row per column of `x`",
                call. = FALSE)
        }
        if (is.null(rownames(x))) rownames(x) <- paste0("feat_", seq_len(nrow(x)))

        x <- as.matrix(x[, keep, drop = FALSE])
        cdata <- data.frame(wellKey = colnames(x), check.names = FALSE)
        cdata$group <- stats::relevel(factor(lab[keep]), ref)
        if (!is.null(covariates)) {
            cdata <- cbind(cdata,
                as.data.frame(covariates)[keep, , drop = FALSE])
        }
        fdata <- data.frame(primerid = rownames(x))
        sca <- MAST::FromMatrix(exprsArray = x, cData = cdata, fData = fdata,
            check_sanity = FALSE)
        rhs <- c("group", names(covariates))
        fml <- stats::as.formula(paste("~ 1 +", paste(rhs, collapse = " + ")))
        contrast_name <- paste0("group", test)
        fit <- function() {
            z <- MAST::zlm(fml, sca = sca, method = param$method,
                ebayes = param$ebayes, silent = !isTRUE(verbose))
            MAST::summary(z, doLRT = contrast_name)$datatable
        }
        # MAST reports its progress with message()
        sm <- if (isTRUE(verbose)) fit() else suppressMessages(fit())

        # GiottoConnect does not import data.table: index as data.frames
        sm <- as.data.frame(sm, check.names = FALSE)
        on_c <- sm$contrast == contrast_name
        h <- sm[on_c & sm$component == "H", c("primerid", "Pr(>Chisq)")]
        fc <- sm[on_c & sm$component == "logFC",
            c("primerid", "coef", "ci.hi", "ci.lo")]
        res <- merge(h, fc, by = "primerid")
        res$fdr <- stats::p.adjust(res[["Pr(>Chisq)"]], "fdr")
        res <- res[order(res$fdr), ]
        res$cluster <- paste0(test, "_vs_", ref)
        res <- data.table::as.data.table(res)
        data.table::setnames(res, "primerid", "feats")
        res
    }
)

#' @title Run MAST differential expression
#' @name runMAST
#' @description
#' Runs MAST (Finak et al. 2015) to find marker features between clusters of
#' a giotto object and returns its per-feature tests (see [mastParam()] for
#' the columns). The giotto object is not changed.
#'
#' * With `group_1` and `group_2`, one comparison: the cells of the
#'   `group_1` clusters against those of the `group_2` clusters. All features
#'   are returned.
#' * Without them, each cluster against all other clusters (those in
#'   `subset_clusters`, or all). Per cluster, the `min_feats` best features are
#'   kept and any other with `fdr < pval` and `coef > logFC`; the `ranking`
#'   column gives the rank within the cluster.
#'
#' MAST fits a model per feature on a dense matrix, so it is slow on large
#' data. **If you use it, please cite the original paper** (see
#' [mastParam()]).
#' @param gobject giotto object
#' @param cluster_column cell metadata column with the clusters
#' @param spat_unit,feat_type spatial unit and feature type
#' @param expression_values expression values to use (MAST expects
#'   log-scale values)
#' @param group_1,group_2 clusters of the two sides of a single comparison
#' @param group_1_name,group_2_name names of the two sides (default: the
#'   cluster IDs joined by `_`)
#' @param subset_clusters clusters compared one against the others (default:
#'   all)
#' @param adjust_columns cell metadata columns to adjust for (e.g. the
#'   detection rate)
#' @param pval,logFC,min_feats filters of the one-against-the-others results
#' @param param a [mastParam()]; its defaults when `NULL`
#' @param verbose be verbose
#' @returns A `data.table` of MAST's results.
#' @inherit mastParam references
#' @seealso [mastParam()]
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' # cluster 1 against clusters 2 and 3
#' runMAST(g, cluster_column = "leiden_clus", group_1 = 1, group_2 = c(2, 3))
#' # every cluster against the others
#' runMAST(g, cluster_column = "leiden_clus")
#' }
#' @md
#' @export
runMAST <- function(gobject, cluster_column, spat_unit = NULL,
    feat_type = NULL, expression_values = "normalized", group_1 = NULL,
    group_2 = NULL, group_1_name = NULL, group_2_name = NULL,
    subset_clusters = NULL, adjust_columns = NULL, pval = 0.001, logFC = 1,
    min_feats = 10L, param = NULL, verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_string(cluster_column)
    checkmate::assert_character(adjust_columns, null.ok = TRUE)
    param <- param %null% mastParam()
    checkmate::assert_class(param, "mastParam")
    if (xor(is.null(group_1), is.null(group_2))) {
        stop("[runMAST] give both `group_1` and `group_2`, or neither",
            call. = FALSE)
    }
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject, spat_unit = spat_unit,
        feat_type = feat_type)

    meta <- pDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
    miss <- setdiff(c(cluster_column, adjust_columns), names(meta))
    if (length(miss)) {
        stop("[runMAST] no cell metadata column(s) ",
            paste(miss, collapse = ", "), call. = FALSE)
    }
    x <- getExpression(gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = expression_values, output = "matrix")
    meta <- meta[match(colnames(x), meta$cell_ID), ]
    clus <- meta[[cluster_column]]
    covariates <- if (!is.null(adjust_columns)) {
        as.data.frame(meta)[, adjust_columns, drop = FALSE]
    }

    run <- function(side_1, side_2, name_1, name_2) {
        lab <- rep(NA_character_, length(clus))
        lab[clus %in% side_1] <- name_1
        lab[clus %in% side_2] <- name_2
        names(lab) <- colnames(x)
        analyzeData(x, param, groups = lab, test = name_1, ref = name_2,
            covariates = covariates, verbose = verbose)
    }

    vmsg(.v = verbose, "runMAST: running MAST",
        as.character(utils::packageVersion("MAST")), "on", nrow(x),
        "features and", ncol(x), "cells")
    if (!is.null(group_1)) {
        res <- run(group_1, group_2,
            group_1_name %null% paste0(group_1, collapse = "_"),
            group_2_name %null% paste0(group_2, collapse = "_"))
    } else {
        lv <- unique(clus[!is.na(clus)])
        if (!is.null(subset_clusters)) lv <- lv[lv %in% subset_clusters]
        lv <- mixedsort(lv)
        if (length(lv) < 2L) {
            stop("[runMAST] needs at least two clusters", call. = FALSE)
        }
        res <- data.table::rbindlist(lapply(lv, function(cl) {
            vmsg(.v = verbose, "runMAST: cluster", cl)
            run(cl, setdiff(lv, cl), as.character(cl), "others")
        }))
        res$ranking <- stats::ave(seq_len(nrow(res)), res$cluster,
            FUN = seq_along)
        keep <- res$ranking <= min_feats | (res$fdr < pval & res$coef > logFC)
        res <- res[keep %in% TRUE, ]
    }
    attr(res, "param") <- param
    citeMethod(param, what = "runMAST()")
    res
}
