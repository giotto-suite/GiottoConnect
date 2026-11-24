# BayesSpace (Zhao et al. 2021, Nature Biotechnology) connection. BayesSpace
# runs its own code; this file builds a SingleCellExperiment from a giotto
# object (counts, array grid positions, optionally a stored PCA), runs
# BayesSpace's spatialPreprocess() and spatialCluster(), and writes the spatial
# clusters to cell metadata (adr/0008 in GiottoClass). GiottoClass has no
# clustering param class, so the param extends analyzeParam for its `$`
# accessors and dispatches on clusterData().

#' @title BayesSpace parameters
#' @name bayesspaceParam
#' @description
#' Parameters for [GiottoClass::clusterData()] that run BayesSpace (Zhao et
#' al. 2021), which clusters spots of array-based spatial transcriptomics
#' (Visium, ST, Visium HD bins) with a Bayesian model whose Potts prior
#' favours the same cluster for neighbouring spots. BayesSpace is a
#' Bioconductor package (`BiocManager::install("BayesSpace")`).
#' **If you use it, please cite the original paper** (see References).
#' @param q integer, number of clusters (no default in BayesSpace; its
#'   `qTune()` helps choose one)
#' @param platform `"Visium"` (hexagonal grid, 6 neighbours), `"ST"` or
#'   `"VisiumHD"` (square grid, 4 neighbours)
#' @param d integer, number of principal components used
#' @param n_pcs,n_hvgs integer, components and highly variable genes of
#'   BayesSpace's preprocessing (`spatialPreprocess()`); not used when a
#'   stored dimension reduction is given
#' @param log_normalize logical, log-normalize counts before the PCA
#' @param init_method `"mclust"` or `"kmeans"`, initial clustering
#' @param model `"t"` or `"normal"`, error model
#' @param precision `"equal"` or `"variable"`, cluster covariance
#' @param nrep,burn_in,thin integer, MCMC iterations, burn-in and thinning
#' @param gamma numeric, smoothing strength; `NULL` uses BayesSpace's default
#'   (3 for Visium, 2 otherwise)
#' @param mu0,lambda0 prior mean and precision of the cluster means; `NULL`
#'   uses BayesSpace's defaults
#' @param alpha,beta numeric, priors of the t-model scaling
#' @param seed integer seed set before preprocessing and clustering; `NULL`
#'   leaves the random state alone
#' @returns A `bayesspaceParam` object.
#' @references Zhao E, Stone MR, Ren X, Guenthoer J, Smythe KS, Pulliam T,
#'   Williams SR, Uytingco CR, Taylor SEB, Nghiem P, Bielas JH, Gottardo R.
#'   Spatial transcriptomics at subspot resolution with BayesSpace. *Nature
#'   Biotechnology* 39, 1375-1384 (2021). \doi{10.1038/s41587-021-00935-2}
#' @seealso [runBayesSpace()]
#' @examples
#' p <- bayesspaceParam(q = 7, nrep = 10000)
#' p$q
#' @md
#' @export
bayesspaceParam <- function(q, platform = c("Visium", "ST", "VisiumHD"),
    d = 15L, n_pcs = 15L, n_hvgs = 2000L, log_normalize = TRUE,
    init_method = c("mclust", "kmeans"), model = c("t", "normal"),
    precision = c("equal", "variable"), nrep = 50000L, burn_in = 1000L,
    thin = 100L, gamma = NULL, mu0 = NULL, lambda0 = NULL, alpha = 1,
    beta = 0.01, seed = 1234L) {
    checkmate::assert_int(q, lower = 1)
    platform <- match.arg(platform)
    checkmate::assert_count(d, positive = TRUE)
    checkmate::assert_count(n_pcs, positive = TRUE)
    checkmate::assert_count(n_hvgs, positive = TRUE)
    checkmate::assert_flag(log_normalize)
    init_method <- match.arg(init_method)
    model <- match.arg(model)
    precision <- match.arg(precision)
    checkmate::assert_int(nrep, lower = 100)
    checkmate::assert_count(burn_in)
    checkmate::assert_count(thin, positive = TRUE)
    if (burn_in >= nrep) {
        stop("[bayesspaceParam] `burn_in` must be smaller than `nrep`",
            call. = FALSE)
    }
    checkmate::assert_number(gamma, null.ok = TRUE)
    checkmate::assert_number(alpha)
    checkmate::assert_number(beta)
    checkmate::assert_int(seed, null.ok = TRUE)
    p <- methods::new("bayesspaceParam", param = list())
    p$q <- as.integer(q)
    p$platform <- platform
    p$d <- as.integer(d)
    p$n_pcs <- as.integer(n_pcs)
    p$n_hvgs <- as.integer(n_hvgs)
    p$log_normalize <- log_normalize
    p$init_method <- init_method
    p$model <- model
    p$precision <- precision
    p$nrep <- as.integer(nrep)
    p$burn_in <- as.integer(burn_in)
    p$thin <- as.integer(thin)
    p$gamma <- gamma
    p$mu0 <- mu0
    p$lambda0 <- lambda0
    p$alpha <- alpha
    p$beta <- beta
    p$seed <- if (is.null(seed)) NULL else as.integer(seed)
    p
}

#' @rdname bayesspaceParam
#' @exportClass bayesspaceParam
setClass("bayesspaceParam", contains = "analyzeParam")

#' @rdname bayesspaceParam
#' @param x count matrix, features x cells (dense or sparse)
#' @param spatial_locs `data.frame` with `cell_ID` and either the array grid
#'   positions `array_row`, `array_col` or the coordinates `sdimx`, `sdimy`
#'   (the grid is then derived from the spot spacing)
#' @param pcs optional matrix, cells x dimensions (rownames are cell IDs) to
#'   cluster on instead of BayesSpace's own preprocessing
#' @details `clusterData(x, param, spatial_locs)` returns a `data.table` with
#'   `cell_ID`, `cluster` (BayesSpace's `spatial.cluster`) and `cluster_init`
#'   (the initial non-spatial clustering). Its attributes `pcs` (the
#'   components clustered) and `hvf` (BayesSpace's highly variable features,
#'   when it preprocessed) hold the intermediate results.
#' @export
setMethod("clusterData", signature(x = "allMatrix", param = "bayesspaceParam"),
    function(x, param, spatial_locs, pcs = NULL, ...) {
        package_check("BayesSpace", repository = "Bioc")
        package_check("SingleCellExperiment", repository = "Bioc")
        cells <- colnames(x)
        if (is.null(cells)) {
            stop("[bayesspaceParam] `x` needs cell IDs as column names",
                call. = FALSE)
        }
        grid <- .bayesspace_grid(spatial_locs, param$platform)
        grid <- grid[match(cells, grid$cell_ID), ]
        if (anyNA(grid$cell_ID)) {
            stop("[bayesspaceParam] `spatial_locs` lacks cells of `x`",
                call. = FALSE)
        }

        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = x),
            colData = S4Vectors::DataFrame(array_row = grid$array_row,
                array_col = grid$array_col, row.names = cells))
        if (!is.null(param$seed)) set.seed(param$seed)
        if (is.null(pcs)) {
            sce <- BayesSpace::spatialPreprocess(sce,
                platform = param$platform, n.PCs = param$n_pcs,
                n.HVGs = param$n_hvgs, log.normalize = param$log_normalize)
        } else {
            pcs <- as.matrix(pcs)
            miss <- setdiff(cells, rownames(pcs))
            if (length(miss)) {
                stop("[bayesspaceParam] `pcs` lacks ", length(miss),
                    " cell(s) of `x`", call. = FALSE)
            }
            SingleCellExperiment::reducedDim(sce, "PCA") <- pcs[cells, ,
                drop = FALSE]
            sce <- BayesSpace::spatialPreprocess(sce,
                platform = param$platform, skip.PCA = TRUE)
        }

        if (!is.null(param$seed)) set.seed(param$seed)
        sce <- BayesSpace::spatialCluster(sce, q = param$q,
            use.dimred = "PCA", d = param$d, platform = param$platform,
            init.method = param$init_method, model = param$model,
            precision = param$precision, nrep = param$nrep,
            burn.in = param$burn_in, thin = param$thin, gamma = param$gamma,
            mu0 = param$mu0, lambda0 = param$lambda0, alpha = param$alpha,
            beta = param$beta)

        res <- data.table::data.table(cell_ID = cells,
            cluster = as.integer(sce$spatial.cluster),
            cluster_init = as.integer(sce$cluster.init))
        attr(res, "pcs") <- SingleCellExperiment::reducedDim(sce, "PCA")
        rd <- SummarizedExperiment::rowData(sce)
        if ("is.HVG" %in% names(rd)) {
            attr(res, "hvf") <- rownames(sce)[rd$is.HVG]
        }
        res
    }
)

# array grid positions (array_row, array_col) of each spot. Given ones are
# kept; otherwise they come from the coordinates and the spot spacing: on the
# Visium hexagonal grid neighbours in a row are 2 columns apart and diagonal
# neighbours 1 column and 1 row apart; on the ST / Visium HD square grid every
# neighbour is 1 step apart.
.bayesspace_grid <- function(spatial_locs, platform = "Visium") {
    sl <- data.table::as.data.table(spatial_locs)
    if (all(c("cell_ID", "array_row", "array_col") %in% names(sl))) {
        return(data.table::data.table(cell_ID = as.character(sl$cell_ID),
            array_row = as.integer(sl$array_row),
            array_col = as.integer(sl$array_col)))
    }
    need <- c("cell_ID", "sdimx", "sdimy")
    miss <- setdiff(need, names(sl))
    if (length(miss)) {
        stop("[bayesspaceParam] `spatial_locs` needs array_row/array_col or ",
            "sdimx/sdimy; lacks: ", paste(miss, collapse = ", "),
            call. = FALSE)
    }
    m <- as.matrix(sl[, c("sdimx", "sdimy"), with = FALSE])
    rownames(m) <- sl$cell_ID
    nn <- createNetwork(m, kNNNetworkParam(k = 1L, include_weight = FALSE,
        output = "data.table"))
    pitch <- stats::median(nn$distance)

    fit <- function(a, b, col_step, row_step) {
        ca <- (a - min(a)) / col_step
        rb <- (b - min(b)) / row_step
        list(col = as.integer(round(ca)), row = as.integer(round(rb)),
            err = mean(abs(ca - round(ca))) + mean(abs(rb - round(rb))))
    }
    if (platform == "Visium") {
        cand <- list(
            fit(m[, 1], m[, 2], pitch / 2, pitch * sqrt(3) / 2),
            fit(m[, 2], m[, 1], pitch / 2, pitch * sqrt(3) / 2))
        best <- cand[[which.min(vapply(cand, `[[`, numeric(1), "err"))]]
        parity <- (best$row + best$col) %% 2L
        if (mean(parity == as.integer(names(which.max(table(parity))))) <
            0.95) {
            stop(wrap_txt("[bayesspaceParam] could not derive a hexagonal",
                "spot grid from the coordinates; give `array_row` and",
                "`array_col` (e.g. from the Visium tissue positions)"),
                call. = FALSE)
        }
    } else {
        best <- fit(m[, 1], m[, 2], pitch, pitch)
    }
    if (anyDuplicated(paste(best$row, best$col))) {
        stop(wrap_txt("[bayesspaceParam] spots share a derived grid position;",
            "give `array_row` and `array_col`"), call. = FALSE)
    }
    data.table::data.table(cell_ID = as.character(sl$cell_ID),
        array_row = best$row, array_col = best$col)
}

#' @title Run BayesSpace spatial clustering
#' @name runBayesSpace
#' @description
#' Runs BayesSpace (Zhao et al. 2021) spatial clustering on the counts and
#' spot positions of a giotto object and stores:
#' * cell metadata `<name>` (BayesSpace's spatial cluster) and
#'   `<name>_init` (the initial non-spatial clustering);
#' * when BayesSpace preprocesses the counts (`dim_reduction_name = NULL`),
#'   its PCA as dimension reduction `<name>_pca` (method `"pca"`), with its
#'   highly variable features in `misc$hvf`.
#'
#' The spot grid comes from cell metadata columns `array_row` and
#' `array_col` when present, otherwise from the spatial locations and the
#' spot spacing. Enhanced (subspot) resolution, `spatialEnhance()`, is not
#' wrapped. **If you use it, please cite the original paper** (see
#' [bayesspaceParam()]).
#' @param gobject giotto object
#' @param param a [bayesspaceParam()] (needs `q`)
#' @param spat_unit,feat_type spatial unit and feature type
#' @param expression_values expression values with the counts
#' @param dim_reduction_name optional name of a stored PCA (cells x
#'   dimensions) to cluster on instead of BayesSpace's own preprocessing
#' @param spat_loc_name name of the spatial locations
#' @param name name of the results
#' @param verbose be verbose
#' @returns The giotto object with the results.
#' @inherit bayesspaceParam references
#' @seealso [bayesspaceParam()]
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- runBayesSpace(g, param = bayesspaceParam(q = 7))
#' }
#' @md
#' @export
runBayesSpace <- function(gobject, param, spat_unit = NULL, feat_type = NULL,
    expression_values = "raw", dim_reduction_name = NULL,
    spat_loc_name = "raw", name = "bayesspace", verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_class(param, "bayesspaceParam")
    checkmate::assert_string(dim_reduction_name, null.ok = TRUE)
    checkmate::assert_string(name)
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject, spat_unit = spat_unit,
        feat_type = feat_type)

    x <- getExpression(gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = expression_values, output = "matrix")
    sl <- getSpatialLocations(gobject, spat_unit = spat_unit,
        name = spat_loc_name, output = "data.table")
    meta <- pDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
    if (all(c("array_row", "array_col") %in% names(meta))) {
        sl <- meta[, c("cell_ID", "array_row", "array_col"), with = FALSE]
    }
    pcs <- NULL
    if (!is.null(dim_reduction_name)) {
        pcs <- getDimReduction(gobject, spat_unit = spat_unit,
            feat_type = feat_type, reduction = "cells",
            reduction_method = "pca", name = dim_reduction_name,
            output = "matrix")
    }

    vmsg(.v = verbose, "runBayesSpace: running BayesSpace",
        as.character(utils::packageVersion("BayesSpace")))
    res <- clusterData(x, param, spatial_locs = sl, pcs = pcs)

    new_meta <- data.table::data.table(cell_ID = res$cell_ID)
    new_meta[[name]] <- factor(res$cluster, levels = seq_len(param$q))
    new_meta[[paste0(name, "_init")]] <- factor(res$cluster_init,
        levels = seq_len(param$q))
    gobject <- addCellMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, new_metadata = new_meta, by_column = TRUE,
        column_cell_ID = "cell_ID")

    if (is.null(dim_reduction_name)) {
        dobj <- createDimObj(coordinates = attr(res, "pcs"),
            name = paste0(name, "_pca"), spat_unit = spat_unit,
            feat_type = feat_type, method = "pca", reduction = "cells",
            misc = list(hvf = attr(res, "hvf"), method = "BayesSpace"))
        gobject <- setDimReduction(gobject, dobj, spat_unit = spat_unit,
            feat_type = feat_type, reduction = "cells",
            reduction_method = "pca", name = paste0(name, "_pca"))
    }
    update_giotto_params(gobject, description = "_bayesspace",
        attachments = list(param = param))
}
