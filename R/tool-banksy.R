# BANKSY (Singhal et al. 2024, Nature Genetics) connection: runs the
# Bioconductor package Banksy on a giotto object. Moved from Giotto core
# (adr/0008 in GiottoClass). The same embedding is computed natively by
# Giotto::runNeighborhoodPCA(); both are reduceData() verbs, so they can be
# compared on equal terms (same inputs, defaults and clustering).

#' @title BANKSY parameters
#' @name banksyParam
#' @description
#' Parameters for [GiottoClass::reduceData()] that run BANKSY (Singhal et al.
#' 2024) through the Bioconductor package \pkg{Banksy}: each cell's values are
#' joined by the mean (and, with `compute_agf = TRUE`, the azimuthal
#' gradient) of its spatial neighbours, weighted by `lambda`, and reduced by
#' PCA. \pkg{Banksy} is suggested and checked when this runs.
#'
#' \pkg{Banksy} takes the absolute value of the neighbour means, which is
#' harmless for expression (never negative) but discards the sign of signed
#' inputs such as PCs. Use expression as the input, or the native
#' `Giotto::runNeighborhoodPCA()`, which keeps the sign.
#'
#' **If you use it, please cite the original paper** (see References).
#' @param lambda numeric in `[0, 1]`, weight of the neighbourhood (about 0.8
#'   for spatial domains, about 0.2 for cell typing)
#' @param k_geom integer (or one per harmonic with `compute_agf = TRUE`),
#'   spatial neighbours
#' @param npcs integer, principal components
#' @param compute_agf logical, also use the azimuthal Gabor filter features
#' @param seed integer, seed of the PCA
#' @returns A `banksyParam` object.
#' @references Singhal V, Chou N, Lee J, et al. BANKSY unifies cell typing
#'   and tissue domain segmentation for scalable spatial omics data analysis.
#'   *Nature Genetics* 56, 431–441 (2024). \doi{10.1038/s41588-024-01664-3}
#' @seealso [runBanksy()]
#' @examples
#' p <- banksyParam(lambda = 0.8)
#' @md
#' @export
banksyParam <- function(lambda = 0.8, k_geom = 15L, npcs = 20L,
    compute_agf = FALSE, seed = 1234L) {
    checkmate::assert_number(lambda, lower = 0, upper = 1)
    checkmate::assert_integerish(k_geom, lower = 1, any.missing = FALSE,
        min.len = 1L, max.len = 2L)
    checkmate::assert_count(npcs, positive = TRUE)
    checkmate::assert_flag(compute_agf)
    checkmate::assert_count(seed)
    p <- methods::new("banksyParam", param = list())
    p$lambda <- lambda
    p$k_geom <- as.integer(k_geom)
    p$npcs <- as.integer(npcs)
    p$compute_agf <- compute_agf
    p$seed <- as.integer(seed)
    p
}

#' @rdname banksyParam
#' @exportClass banksyParam
setClass("banksyParam", contains = "reduceParam")

#' @rdname banksyParam
#' @param x features x cells values with cell IDs as column names
#' @param param a `banksyParam`
#' @param ... unused
#' @param spatial_locs cells x 2 coordinates with cell IDs as row names
#' @details `reduceData(x, param, spatial_locs = )` returns a list with `u`
#'   (cells x `npcs` BANKSY PCs) and `se`, the `SpatialExperiment` Banksy
#'   worked on (for `Banksy::clusterBanksy()`).
#' @export
setMethod("reduceData", signature(x = "ANY", param = "banksyParam"),
    function(x, param, ..., spatial_locs) {
        package_check(c("SpatialExperiment", "Banksy"),
            repository = c("Bioc", "Bioc"))
        x <- as.matrix(x)
        if (any(x < 0, na.rm = TRUE)) {
            warning(wrap_txt("[banksyParam] the input has negative values",
                "(e.g. PCs). Banksy takes the absolute value of the",
                "neighbour means, which discards their sign; use expression",
                "or Giotto::runNeighborhoodPCA()."), call. = FALSE)
        }
        locs <- as.matrix(spatial_locs)[colnames(x), , drop = FALSE]
        if (anyNA(locs)) {
            stop("[banksyParam] some cells have no spatial location",
                call. = FALSE)
        }
        p <- param@param
        npcs <- min(p$npcs, nrow(x) * (1L + length(p$k_geom) * p$compute_agf))
        se <- SpatialExperiment::SpatialExperiment(assays = list(values = x),
            spatialCoords = locs)
        se <- Banksy::computeBanksy(se, assay_name = "values",
            compute_agf = p$compute_agf, k_geom = p$k_geom, verbose = FALSE)
        se <- Banksy::runBanksyPCA(se, use_agf = p$compute_agf,
            lambda = p$lambda, npcs = npcs, seed = p$seed)
        rd <- getExportedValue("SingleCellExperiment", "reducedDimNames")(se)
        u <- getExportedValue("SingleCellExperiment", "reducedDim")(se,
            grep("^PCA_", rd, value = TRUE)[1L])
        dimnames(u) <- list(colnames(x), paste0("Dim.", seq_len(ncol(u))))
        list(u = u, se = se)
    }
)

#' @title Spatial domains with BANKSY
#' @name runBanksy
#' @description
#' Runs BANKSY (Singhal et al. 2024) through the Bioconductor package
#' \pkg{Banksy} (see [banksyParam()]), stores its PCs as a dimension reduction
#' (method `"banksy"`) and clusters them into domains (cell metadata `name`):
#' by default with Giotto's sNN network and Leiden, as
#' `Giotto::findSpatialDomains(method = "neighborhood")` does with the native
#' embedding, so the two can be compared directly.
#'
#' The input is expression by default, which is what \pkg{Banksy} is designed
#' for. With `input = "dim_reduction"` (e.g. PCs) \pkg{Banksy} discards the
#' sign of the neighbour means (see [banksyParam()]) and a warning is given;
#' `Giotto::runNeighborhoodPCA()` keeps it.
#'
#' **If you use it, please cite the original paper** (see References).
#' @param gobject giotto object
#' @param spat_unit,feat_type spatial unit and feature type
#' @param input `"expression"` (default) or `"dim_reduction"`
#' @param expression_values expression values with `input = "expression"`
#' @param feats_to_use features with `input = "expression"`: a character
#'   vector or the name of a logical feature metadata column (default
#'   `"hvf"`); `NULL` uses all
#' @param dim_reduction_to_use,dim_reduction_name,dimensions_to_use the
#'   embedding with `input = "dim_reduction"`
#' @param spat_loc_name name of the spatial locations
#' @param lambda,k_geom,npcs,compute_agf,seed see [banksyParam()]
#' @param cluster_method `"giotto"` (sNN network + Leiden), `"banksy"`
#'   (`Banksy::clusterBanksy()`) or `"none"`
#' @param k_neighbors neighbours of the sNN network
#' @param resolution Leiden resolution
#' @param name name of the PCs, the network and the domain column
#' @param verbose be verbose
#' @returns The giotto object with the BANKSY PCs and the domains.
#' @inherit banksyParam references
#' @seealso [banksyParam()], `Giotto::runNeighborhoodPCA()`,
#'   `Giotto::findSpatialDomains()`
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- runBanksy(g, lambda = 0.8)
#' }
#' @md
#' @export
runBanksy <- function(gobject, spat_unit = NULL, feat_type = NULL,
    input = c("expression", "dim_reduction"), expression_values = "normalized",
    feats_to_use = "hvf", dim_reduction_to_use = "pca",
    dim_reduction_name = NULL, dimensions_to_use = NULL,
    spat_loc_name = "raw", lambda = 0.8, k_geom = 15L, npcs = 20L,
    compute_agf = FALSE, seed = 1234L,
    cluster_method = c("giotto", "banksy", "none"), k_neighbors = 30L,
    resolution = 0.8, name = "banksy", verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    input <- match.arg(input)
    cluster_method <- match.arg(cluster_method)
    checkmate::assert_string(name)
    param <- banksyParam(lambda = lambda, k_geom = k_geom, npcs = npcs,
        compute_agf = compute_agf, seed = seed)
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject, spat_unit = spat_unit,
        feat_type = feat_type)

    x <- if (input == "expression") {
        x <- getExpression(gobject, spat_unit = spat_unit,
            feat_type = feat_type, values = expression_values,
            output = "matrix")
        feats <- .banksy_feats(gobject, feats_to_use, spat_unit, feat_type)
        if (is.null(feats)) x else x[intersect(feats, rownames(x)), ,
            drop = FALSE]
    } else {
        emb <- getDimReduction(gobject, spat_unit = spat_unit,
            feat_type = feat_type, reduction = "cells",
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name, output = "matrix")
        if (!is.null(dimensions_to_use)) {
            dims <- dimensions_to_use[dimensions_to_use <= ncol(emb)]
            if (!length(dims)) {
                stop("[runBanksy] no stored dimension is selected",
                    call. = FALSE)
            }
            emb <- emb[, dims, drop = FALSE]
        }
        t(emb)
    }
    sl <- getSpatialLocations(gobject, spat_unit = spat_unit,
        name = spat_loc_name, output = "data.table")
    locs <- as.matrix(sl[, c("sdimx", "sdimy")])
    rownames(locs) <- sl$cell_ID

    vmsg(.v = verbose, "runBanksy: Banksy neighbourhood and PCA")
    res <- reduceData(x, param, spatial_locs = locs)
    dim_obj <- createDimObj(coordinates = res$u, name = name,
        spat_unit = spat_unit, feat_type = feat_type, method = "banksy",
        reduction = "cells",
        misc = list(input = input, features = rownames(x)))
    gobject <- setDimReduction(gobject, x = dim_obj, spat_unit = spat_unit,
        feat_type = feat_type, reduction = "cells",
        reduction_method = "banksy", name = name, verbose = FALSE)

    if (cluster_method == "giotto") {
        package_check("Giotto", repository = "github:giotto-suite/Giotto")
        quiet <- function(code) gwith_options(
            list(giotto.update_param = FALSE), code)
        gobject <- quiet(Giotto::createNearestNetwork(gobject,
            spat_unit = spat_unit, feat_type = feat_type, type = "sNN",
            dim_reduction_to_use = "banksy", dim_reduction_name = name,
            dimensions_to_use = seq_len(ncol(res$u)), k = k_neighbors,
            name = name, verbose = FALSE))
        gobject <- quiet(Giotto::doLeidenCluster(gobject,
            spat_unit = spat_unit, feat_type = feat_type,
            nn_network_to_use = "sNN", network_name = name, name = name,
            resolution = resolution, n_iterations = -1,
            seed_number = seed))
    } else if (cluster_method == "banksy") {
        se <- Banksy::clusterBanksy(res$se, use_agf = compute_agf,
            lambda = lambda, npcs = ncol(res$u), k_neighbors = k_neighbors,
            resolution = resolution, seed = seed)
        cd <- getExportedValue("SummarizedExperiment", "colData")(se)
        new_meta <- data.frame(cell_ID = colnames(x),
            domain = as.character(cd[[grep("^clust_", colnames(cd))[1L]]]))
        names(new_meta)[2] <- name
        gobject <- addCellMetadata(gobject, spat_unit = spat_unit,
            feat_type = feat_type, new_metadata = new_meta,
            by_column = TRUE, column_cell_ID = "cell_ID")
    }
    update_giotto_params(gobject, description = "_banksy",
        attachments = list(param = param))
}

# features: NULL (all), a character vector, or a logical feature metadata
# column such as "hvf"
.banksy_feats <- function(gobject, feats_to_use, spat_unit, feat_type) {
    if (is.null(feats_to_use)) return(NULL)
    if (length(feats_to_use) == 1L) {
        fm <- fDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
        if (feats_to_use %in% names(fm)) {
            col <- fm[[feats_to_use]]
            keep <- if (is.logical(col)) col else col %in% c("yes", "TRUE")
            return(fm$feat_ID[keep])
        }
    }
    feats_to_use
}
