# AnnData target: prepared parts become an .h5ad file (written with
# anndataR, no python), in scanpy's layout:
#   expression     layer          "X", or a layer name ("counts" for raw)
#   cell_metadata  (none)         obs
#   feat_metadata  (none)         var
#   spatial_locs   obsm           obsm/<obsm> (x, y image-space)
#   dim_reduction  obsm           obsm/<obsm>; PCA loadings -> varm/PCs,
#                                 eigenvalues -> uns/pca/variance
#   nn_network     obsp, nn_type  obsp/<obsp>_distances / _connectivities +
#   spatial_network obsp            a uns/<obsp> record (scanpy keys + giotto)
# The conversion record is stored as JSON in uns/giottoConnect.

.ANNDATA_TARGET <- list(
    expression = "layer",
    cell_metadata = character(),
    feat_metadata = character(),
    spatial_locs = "obsm",
    dim_reduction = "obsm",
    nn_network = c("obsp", "nn_type"),
    spatial_network = "obsp"
)

setClass("anndataTarget", contains = "giottoConnectTarget", prototype = list(contract = "anndata"))

setMethod("targetBuild", "anndataTarget", function(target, prepared, mapping, prov, source = NULL,
    output = NULL) {
    if (!is.null(source)) stop(wrap_txt("`source` applies to Giotto targets only"), call. = FALSE)
    if (is.null(output)) stop(wrap_txt("the AnnData target needs an `output` .h5ad file"), call. = FALSE)
    .write_anndata(prepared, prov, output)
})
setMethod("targetCheckKeys", "anndataTarget", function(target, mapping, out) {
    .vm_targets_anndata(mapping, out)
})
setMethod("targetSemantic", "anndataTarget", function(target, mapping, profile, result, out) {
    .vm_semantic_anndata(mapping, profile, result, out)
})
setMethod("targetComparators", "anndataTarget", function(target) .RB_COMPARE_ANNDATA)

# cells x features, with dimnames
# (a GiottoDisk store is read into a sparse matrix)
.ad_cells_by_feats <- function(m) {
    m <- if (is.matrix(m)) m else if (methods::is(m, "Matrix")) methods::as(m, "CsparseMatrix") else
        methods::as(m, "dgCMatrix")
    Matrix::t(m)
}

# edge table -> cells x cells sparse matrix of one column (rows = from)
.ad_edges_to_matrix <- function(e, value, cells) {
    i <- match(e$from, cells)
    j <- match(e$to, cells)
    ok <- !is.na(i) & !is.na(j)
    v <- e[[value]][ok]
    v[is.na(v)] <- 0
    Matrix::sparseMatrix(i = i[ok], j = j[ok], x = as.numeric(v), dims = c(length(cells), length(cells)),
        dimnames = list(cells, cells))
}

.ad_table_rows <- function(d, id, cells) {
    d <- as.data.frame(d)
    rn <- as.character(d[[id]])
    d[[id]] <- NULL
    d <- d[match(cells, rn), , drop = FALSE]
    rownames(d) <- cells
    d
}

.write_anndata <- function(prepared, prov, output) {
    package_check("anndataR", repository = "Bioc:anndataR")
    types <- vapply(prepared, function(pp) pp$spec$part, character(1))
    ex <- prepared[types == "expression"]
    if (!length(ex)) stop(wrap_txt("AnnData needs at least one expression part"), call. = FALSE)
    mats <- lapply(ex, function(pp) .ad_cells_by_feats(pp$data))
    slots <- vapply(ex, function(pp) pp$spec$to$layer, character(1))
    # X holds normalized values when there are any (scanpy), else the counts
    if (!"X" %in% slots) slots[slots == "counts"][1L] <- "X"
    X <- mats[[which(slots == "X")[1L]]]
    cells <- rownames(X)
    feats <- colnames(X)
    layers <- mats[slots != "X"]
    names(layers) <- slots[slots != "X"]
    layers <- lapply(layers, function(m) m[cells, feats, drop = FALSE])

    obs <- data.frame(row.names = cells)
    var <- data.frame(row.names = feats)
    obsm <- list()
    varm <- list()
    obsp <- list()
    uns <- list()
    for (pp in prepared) {
        to <- pp$spec$to
        d <- pp$data
        switch(pp$spec$part,
            cell_metadata = obs <- .ad_table_rows(d, "cell_ID", cells),
            feat_metadata = var <- .ad_table_rows(d, "feat_ID", feats),
            spatial_locs = {
                d <- as.data.frame(d)
                m <- as.matrix(d[match(cells, d$cell_ID), intersect(c("sdimx", "sdimy", "sdimz"), names(d)),
                    drop = FALSE])
                dimnames(m) <- NULL
                obsm[[to$obsm]] <- m
            },
            dim_reduction = {
                coords <- if (methods::is(d, "dimObj")) d[] else as.matrix(d)
                m <- coords[cells, , drop = FALSE]
                dimnames(m) <- NULL
                obsm[[to$obsm]] <- m
                misc <- if (methods::is(d, "dimObj")) d@misc else NULL
                if (identical(to$obsm, "X_pca") && length(misc$loadings)) {
                    # scanpy keeps loadings for every var (zeros outside the PCA features)
                    load <- matrix(0, length(feats), ncol(misc$loadings))
                    hit <- match(rownames(misc$loadings), feats)
                    load[hit[!is.na(hit)], ] <- misc$loadings[!is.na(hit), , drop = FALSE]
                    varm$PCs <- load
                    if (length(misc$eigenvalues)) {
                        ev <- as.numeric(misc$eigenvalues)
                        uns$pca <- list(variance = ev, variance_ratio = ev / sum(ev))
                    }
                }
            },
            nn_network = ,
            spatial_network = {
                key <- to$obsp
                dk <- paste0(key, "_distances")
                ck <- paste0(key, "_connectivities")
                obsp[[dk]] <- .ad_edges_to_matrix(d, "distance", cells)
                obsp[[ck]] <- .ad_edges_to_matrix(d, "weight", cells)
                uns[[key]] <- list(distances_key = dk, connectivities_key = ck,
                    giotto = list(name = to$name %null% key,
                        type = if (pp$spec$part == "spatial_network") "spatial" else to$nn_type))
            },
            NULL
        )
    }
    uns$giottoConnect <- as.character(jsonlite::toJSON(prov, auto_unbox = TRUE, null = "null",
        na = "null", force = TRUE, digits = NA))
    ad <- anndataR::AnnData(X = X, obs = obs, var = var, layers = if (length(layers)) layers,
        obsm = if (length(obsm)) obsm, varm = if (length(varm)) varm,
        obsp = if (length(obsp)) obsp, uns = uns)
    if (file.exists(output)) unlink(output)
    anndataR::write_h5ad(ad, output, mode = "w")
    normalizePath(output)
}

# conversion record of an .h5ad written by the AnnData target
.anndata_provenance <- function(path) {
    h5 <- .h5ad_open(path)
    on.exit(h5$close_all())
    if (!.h5ad_exists(h5, "uns/giottoConnect")) return(NULL)
    jsonlite::fromJSON(as.character(.h5ad_value(h5, "uns/giottoConnect")), simplifyVector = TRUE,
        simplifyDataFrame = TRUE)
}

.is_h5ad_path <- function(x) is.character(x) && length(x) == 1L && grepl("\\.h5ad$", x) && file.exists(x)

# validation ####

.vm_targets_anndata <- function(mapping, out) {
    for (p in mapping$parts) {
        need <- .ANNDATA_TARGET[[p$part]]
        if (is.null(need)) {
            out$add("target", p$id, "FAIL", sprintf("AnnData has no place for %s", p$part))
            next
        }
        keys <- setdiff(names(p$to), "name")
        miss <- setdiff(need, keys)
        unknown <- setdiff(keys, need)
        lvl <- if (length(miss) || length(unknown)) "FAIL" else "OK"
        out$add("target", p$id, lvl, if (lvl == "OK") sprintf("%s -> AnnData", p$part) else
            sprintf("AnnData %s needs keys %s (missing: %s; unknown: %s)", p$part,
                paste(need, collapse = ", "), paste(miss, collapse = ", "), paste(unknown, collapse = ", ")))
    }
}

.vm_semantic_anndata <- function(mapping, profile, path, out) {
    h5 <- .h5ad_open(path)
    on.exit(h5$close_all())
    cells <- .h5ad_index(h5, "obs")
    if (!is.na(profile$n_cells)) {
        out$add("dims", NA, if (length(cells) == profile$n_cells) "OK" else "WARN",
            sprintf("%d obs in the h5ad, %d cells in source", length(cells), profile$n_cells))
    }
    counts <- if (.h5ad_exists(h5, "layers/counts")) "layers/counts" else "X"
    m <- .h5ad_matrix(h5, counts)
    v <- m[, utils::head(cells, 99L), drop = FALSE]
    v <- if (methods::is(v, "sparseMatrix")) v@x else as.vector(v)
    ok <- all(v == round(v)) && all(v >= 0)
    out$add("raw_counts", counts, if (ok) "OK" else "WARN",
        sprintf("%s is %sinteger-valued and non-negative", counts, if (ok) "" else "NOT "))
    if (.h5ad_exists(h5, "obsm")) {
        for (k in names(h5[["obsm"]])) {
            n <- rev(h5[[paste0("obsm/", k)]]$dims)[1L]
            out$add("obsm", k, if (n == length(cells)) "OK" else "FAIL", sprintf("%d rows for %d obs", n,
                length(cells)))
        }
    }
}

# readback: compare the source parts with the written file ####

.RB_COMPARE_ANNDATA <- list(
    expression = function(spec, data, path, n_cells) {
        data <- .rb_as_matrix(data)
        h5 <- .h5ad_open(path)
        on.exit(h5$close_all())
        el <- if (.h5ad_exists(h5, paste0("layers/", spec$to$layer))) paste0("layers/", spec$to$layer) else "X"
        got <- .h5ad_matrix(h5, el)
        use <- .rb_sample(colnames(data), n_cells)
        list(
            .rb_row("placement", TRUE, el),
            .rb_row("ids", setequal(colnames(got), colnames(data)) && setequal(rownames(got), rownames(data)),
                sprintf("%d obs, %d var", ncol(got), nrow(got))),
            .rb_row("values", .rb_same(as.numeric(as.matrix(got[rownames(data), use])),
                as.numeric(as.matrix(data[, use]))), sprintf("%d sampled cells", length(use)))
        )
    },
    cell_metadata = function(spec, data, path, n_cells) {
        h5 <- .h5ad_open(path)
        on.exit(h5$close_all())
        obs <- .h5ad_frame(h5, "obs")
        .rb_table(data, cbind(cell_ID = rownames(obs), obs), "cell_ID")
    },
    feat_metadata = function(spec, data, path, n_cells) {
        h5 <- .h5ad_open(path)
        on.exit(h5$close_all())
        var <- .h5ad_frame(h5, "var")
        .rb_table(data, cbind(feat_ID = rownames(var), var), "feat_ID")
    },
    spatial_locs = function(spec, data, path, n_cells) {
        h5 <- .h5ad_open(path)
        on.exit(h5$close_all())
        m <- .h5ad_obsm(h5, paste0("obsm/", spec$to$obsm))
        d <- as.data.frame(data)
        d <- d[match(rownames(m), d$cell_ID), , drop = FALSE]
        list(
            .rb_row("placement", TRUE, paste0("obsm/", spec$to$obsm)),
            .rb_row("values", .rb_same(c(m[, 1], m[, 2]), c(d$sdimx, d$sdimy)), sprintf("%d cells", nrow(m)))
        )
    },
    dim_reduction = function(spec, data, path, n_cells) {
        h5 <- .h5ad_open(path)
        on.exit(h5$close_all())
        m <- .h5ad_obsm(h5, paste0("obsm/", spec$to$obsm))
        coords <- if (methods::is(data, "dimObj")) data[] else as.matrix(data)
        list(
            .rb_row("placement", TRUE, paste0("obsm/", spec$to$obsm)),
            .rb_row("values", .rb_same(as.numeric(m), as.numeric(coords[rownames(m), , drop = FALSE])), "")
        )
    },
    nn_network = function(spec, data, path, n_cells) .rb_anndata_network(spec, data, path),
    spatial_network = function(spec, data, path, n_cells) .rb_anndata_network(spec, data, path)
)

.rb_anndata_network <- function(spec, data, path) {
    h5 <- .h5ad_open(path)
    on.exit(h5$close_all())
    m <- .h5ad_obsp(h5, paste0("obsp/", spec$to$obsp, "_distances"))
    tr <- Matrix::summary(methods::as(m, "CsparseMatrix"))
    got <- paste(rownames(m)[tr$i], colnames(m)[tr$j])
    d <- as.data.frame(data)
    exp <- paste(d$from, d$to)
    dist <- tr$x[match(exp, got)]
    list(
        .rb_row("placement", TRUE, paste0("obsp/", spec$to$obsp, "_distances")),
        .rb_row("edges", setequal(got, exp), sprintf("%d edges", length(exp))),
        .rb_row("distance", .rb_same(dist, d$distance), "", field = "distance")
    )
}
