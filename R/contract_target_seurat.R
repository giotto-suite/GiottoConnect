# Seurat as a mapping target. The engine prepares parts as for Giotto; here
# they are built with SeuratObject constructors and assembled into a Seurat
# object (expression first, it creates the object). Target keys per part:
#   expression     assay, layer          (layer "counts" creates the assay)
#   cell_metadata  -                     (table with a cell_ID column)
#   feat_metadata  assay                 (table with a feat_ID column)
#   spatial_locs   fov, boundary         (x, y, cell; centroids)
#   spatial_info   fov, boundary         (x, y, cell vertices; segmentation)
#   dim_reduction  reduction, key        (cells x dims matrix, or a dimObj with
#                                          loadings / eigenvalues in misc)
#   nn_network     graph [, nn_type]     (cells x cells sparse adjacency)
#   feature_info   fov                   (x, y, gene; FOV molecules)

.SEURAT_TARGET <- list(
    expression = c("assay", "layer"),
    cell_metadata = character(),
    feat_metadata = c("assay"),
    spatial_locs = c("fov", "boundary"),
    spatial_info = c("fov", "boundary"),
    dim_reduction = c("reduction", "key"),
    nn_network = c("graph"),
    feature_info = c("fov")
)
# keys a part may carry besides the required ones
.SEURAT_TARGET_OPTIONAL <- list(nn_network = "nn_type")

.assemble_seurat <- function(prepared) {
    package_check("SeuratObject")
    SO <- asNamespace("SeuratObject")
    by_part <- function(type) Filter(function(pp) pp$spec$part == type, prepared)

    expr <- by_part("expression")
    counts <- Filter(function(pp) identical(pp$spec$to$layer, "counts"), expr)
    if (!length(counts)) stop(wrap_txt("A Seurat target needs an expression part with layer 'counts'"),
        call. = FALSE)
    s <- NULL
    for (pp in counts) {
        a <- SO$CreateAssay5Object(counts = pp$data)
        if (is.null(s)) {
            s <- SO$CreateSeuratObject(a, assay = pp$spec$to$assay)
        } else {
            s[[pp$spec$to$assay]] <- a
        }
    }
    for (pp in Filter(function(pp) !identical(pp$spec$to$layer, "counts"), expr)) {
        s <- SO$`LayerData<-`(s, layer = pp$spec$to$layer, assay = pp$spec$to$assay, value = pp$data)
    }

    for (pp in by_part("cell_metadata")) {
        md <- as.data.frame(pp$data)
        rownames(md) <- md$cell_ID
        md$cell_ID <- NULL
        s <- SO$AddMetaData(s, metadata = md[colnames(s), , drop = FALSE])
    }
    for (pp in by_part("feat_metadata")) {
        fm <- as.data.frame(pp$data)
        rownames(fm) <- fm$feat_ID
        fm$feat_ID <- NULL
        if (ncol(fm)) {
            a <- s[[pp$spec$to$assay]]
            a <- SO$AddMetaData(a, metadata = fm[rownames(a), , drop = FALSE])
            s[[pp$spec$to$assay]] <- a
        }
    }

    # one FOV per fov name, from its centroids and/or segmentation parts and
    # transcripts (molecules)
    geo <- c(by_part("spatial_locs"), by_part("spatial_info"))
    mols <- by_part("feature_info")
    for (fov in unique(vapply(c(geo, mols), function(pp) pp$spec$to$fov, character(1)))) {
        bnd <- list()
        mol <- Filter(function(pp) pp$spec$to$fov == fov, mols)
        mol <- if (length(mol)) as.data.frame(mol[[1L]]$data)[, c("x", "y", "gene")]
        for (pp in Filter(function(pp) pp$spec$to$fov == fov, geo)) {
            df <- as.data.frame(pp$data)[, c("x", "y", "cell")]
            df <- df[df$cell %in% colnames(s), , drop = FALSE]
            bnd[[pp$spec$to$boundary]] <- if (pp$spec$part == "spatial_locs") {
                SO$CreateCentroids(df)
            } else {
                SO$CreateSegmentation(df)
            }
        }
        type <- vapply(bnd, function(b) if (methods::is(b, "Centroids")) "centroids" else "segmentation",
            character(1))
        s[[fov]] <- SO$CreateFOV(bnd, type = unname(type), molecules = mol, assay = SO$DefaultAssay(s),
            key = paste0(gsub("[^A-Za-z0-9]", "", fov), "fov_"))
    }

    for (pp in by_part("dim_reduction")) {
        d <- pp$data
        m <- if (methods::is(d, "dimObj")) d[] else as.matrix(d)
        misc <- if (methods::is(d, "dimObj")) d@misc else NULL
        # Seurat keys are alphanumeric followed by "_"
        key <- paste0(gsub("[^A-Za-z0-9]", "", sub("_$", "", pp$spec$to$key)), "_")
        colnames(m) <- paste0(key, seq_len(ncol(m)))
        load <- misc$loadings
        if (length(load)) {
            load <- as.matrix(load)[, seq_len(min(ncol(load), ncol(m))), drop = FALSE]
            rownames(load) <- .seurat_feat(rownames(load))
            colnames(load) <- colnames(m)[seq_len(ncol(load))]
        }
        ev <- misc$eigenvalues
        s[[pp$spec$to$reduction]] <- SO$CreateDimReducObject(embeddings = m[colnames(s), , drop = FALSE],
            loadings = if (length(load)) load else new("matrix"),
            stdev = if (length(ev)) sqrt(as.numeric(ev))[seq_len(ncol(m))] else numeric(),
            key = key, assay = SO$DefaultAssay(s))
    }
    for (pp in by_part("nn_network")) {
        g <- SO$as.Graph(pp$data[colnames(s), colnames(s)])
        methods::slot(g, "assay.used") <- SO$DefaultAssay(s)
        s[[pp$spec$to$graph]] <- g
    }
    s
}

.record_seurat_conversion <- function(s, prov) {
    misc <- methods::slot(s, "misc")
    misc$giottoConnect <- prov
    methods::slot(s, "misc") <- misc
    s
}


# every part carries exactly the target keys its Seurat part needs
.vm_targets_seurat <- function(mapping, out) {
    for (p in mapping$parts) {
        need <- .SEURAT_TARGET[[p$part]]
        miss <- setdiff(need, names(p$to))
        unknown <- setdiff(names(p$to), c(need, .SEURAT_TARGET_OPTIONAL[[p$part]]))
        lvl <- if (length(miss) || length(unknown)) "FAIL" else "OK"
        out$add("target", p$id, lvl, if (lvl == "OK") sprintf("%s -> Seurat", p$part) else
            sprintf("Seurat %s needs keys %s (missing: %s; unknown: %s)", p$part,
                paste(need, collapse = ", "), paste(miss, collapse = ", "), paste(unknown, collapse = ", ")))
    }
}


# semantic checks on a Seurat result ####

.vm_semantic_seurat <- function(mapping, profile, s, out) {
    SO <- asNamespace("SeuratObject")
    cells <- colnames(s)
    for (a in SO$Assays(s)) {
        if (!"counts" %in% SO$Layers(s[[a]])) next
        m <- SO$LayerData(s, layer = "counts", assay = a)
        v <- m[, utils::head(cells, 99L), drop = FALSE]
        v <- if (methods::is(v, "sparseMatrix")) v@x else as.vector(v)
        ok <- all(v == round(v)) && all(v >= 0)
        out$add("raw_counts", a, if (ok) "OK" else "FAIL",
            sprintf("counts are %sinteger-valued and non-negative", if (ok) "" else "NOT "))
        if (!is.na(profile$n_feat)) {
            out$add("dims", a, if (nrow(m) == profile$n_feat) "OK" else "WARN",
                sprintf("%d features in Seurat, %d in source", nrow(m), profile$n_feat))
        }
    }
    if (!is.na(profile$n_cells)) {
        out$add("dims", NA, if (length(cells) == profile$n_cells) "OK" else "WARN",
            sprintf("%d cells in Seurat, %d in source", length(cells), profile$n_cells))
    }
    for (im in SO$Images(s)) {
        for (b in SO$Boundaries(s[[im]])) {
            ct <- SO$GetTissueCoordinates(s[[im]][[b]])
            ov <- mean(cells %in% ct$cell)
            out$add("cell_ids", paste(im, b, sep = "/"), if (ov >= 0.99) "OK" else "FAIL",
                sprintf("%.1f%% of cells found in %s/%s", 100 * ov, im, b))
            rng <- max(diff(range(ct$x)), diff(range(ct$y)))
            out$add("coord_units", paste(im, b, sep = "/"), if (rng < 50 || rng > 1e6) "WARN" else "OK",
                sprintf("coordinate range %.3g", rng))
        }
    }
    for (r in SO$Reductions(s)) {
        e <- SO$Embeddings(s, reduction = r)
        out$add("dim_reduction", r, if (nrow(e) == length(cells)) "OK" else "FAIL",
            sprintf("%d rows for %d cells", nrow(e), length(cells)))
    }
    for (gname in SO$Graphs(s)) {
        g <- s[[gname]]
        ok <- identical(rownames(g), cells) && identical(colnames(g), cells)
        out$add("nn_network", gname, if (ok) "OK" else "FAIL",
            sprintf("%d x %d, indexed by the object's cells: %s", nrow(g), ncol(g), ok))
    }
}


# readback comparators for a Seurat result ####

# Seurat replaces "_" in feature names with "-"
.seurat_feat <- function(x) gsub("_", "-", x, fixed = TRUE)

.RB_COMPARE_SEURAT <- list(
    expression = function(spec, data, s, n_cells) {
        to <- spec$to
        got <- asNamespace("SeuratObject")$LayerData(s, layer = to$layer, assay = to$assay)
        use <- .rb_sample(colnames(data), n_cells)
        rn <- .seurat_feat(rownames(data))
        list(
            .rb_row("placement", TRUE, sprintf("%s/%s", to$assay, to$layer)),
            .rb_row("ids", setequal(colnames(got), colnames(data)) && setequal(rownames(got), rn),
                sprintf("%d cells, %d features", ncol(got), nrow(got))),
            .rb_row("values", .rb_same(as.numeric(as.matrix(got[rn, use, drop = FALSE])),
                as.numeric(as.matrix(data[, use, drop = FALSE]))), sprintf("%d sampled cells", length(use)))
        )
    },
    cell_metadata = function(spec, data, s, n_cells) {
        md <- s[[]]
        md <- cbind(cell_ID = rownames(md), md)
        .rb_table(data, md, "cell_ID")
    },
    feat_metadata = function(spec, data, s, n_cells) {
        fm <- s[[spec$to$assay]][[]]
        exp <- as.data.frame(data)
        exp$feat_ID <- .seurat_feat(exp$feat_ID)
        .rb_table(exp, cbind(feat_ID = rownames(fm), fm), "feat_ID")
    },
    spatial_locs = function(spec, data, s, n_cells) {
        ct <- asNamespace("SeuratObject")$GetTissueCoordinates(s[[spec$to$fov]][[spec$to$boundary]])
        .rb_table(as.data.frame(data)[, c("cell", "x", "y")], ct[, c("cell", "x", "y")], "cell")
    },
    spatial_info = function(spec, data, s, n_cells) {
        ct <- asNamespace("SeuratObject")$GetTissueCoordinates(s[[spec$to$fov]][[spec$to$boundary]])
        bb <- function(d) {
            r <- lapply(split(d[, c("x", "y")], d$cell), function(v) {
                c(xmin = min(v$x), xmax = max(v$x), ymin = min(v$y), ymax = max(v$y))
            })
            data.frame(cell = names(r), do.call(rbind, r), row.names = NULL)
        }
        m <- merge(bb(as.data.frame(data)), bb(ct), by = "cell")
        moved <- pmax(abs(m$xmin.x - m$xmin.y), abs(m$xmax.x - m$xmax.y),
            abs(m$ymin.x - m$ymin.y), abs(m$ymax.x - m$ymax.y)) > 1e-8
        list(
            .rb_row("placement", TRUE, sprintf("%s/%s", spec$to$fov, spec$to$boundary)),
            .rb_row("ids", setequal(unique(ct$cell), unique(data$cell)),
                sprintf("%d polygons", length(unique(ct$cell)))),
            .rb_row("geometry", !any(moved), sprintf("%d of %d polygons with a different extent",
                sum(moved), nrow(m)))
        )
    },
    dim_reduction = function(spec, data, s, n_cells) {
        SO <- asNamespace("SeuratObject")
        e <- SO$Embeddings(s, reduction = spec$to$reduction)
        misc <- if (methods::is(data, "dimObj")) data@misc else NULL
        data <- if (methods::is(data, "dimObj")) data[] else as.matrix(data)
        list(
            .rb_row("placement", TRUE, spec$to$reduction),
            .rb_row("ids", setequal(rownames(e), rownames(data)), sprintf("%d cells", nrow(e))),
            .rb_row("values", .rb_same(as.numeric(e[rownames(data), ]), as.numeric(data)), ""),
            if (length(misc$loadings)) {
                l <- SO$Loadings(s, reduction = spec$to$reduction)
                ld <- as.matrix(misc$loadings)
                .rb_row("loadings", .rb_same(as.numeric(l[.seurat_feat(rownames(ld)), seq_len(ncol(l))]),
                    as.numeric(ld[, seq_len(ncol(l))])), sprintf("%d features", nrow(ld)))
            }
        )
    },
    feature_info = function(spec, data, s, n_cells) {
        SO <- asNamespace("SeuratObject")
        fov <- s[[spec$to$fov]]
        got <- SO$GetTissueCoordinates(fov[[names(methods::slot(fov, "molecules"))[1L]]])
        a <- table(data$gene)
        b <- table(got$molecule)
        list(
            .rb_row("placement", TRUE, paste0(spec$to$fov, "/molecules")),
            .rb_row("ids", identical(sort(names(a)), sort(names(b))) && all(a == b[names(a)]),
                sprintf("%d transcripts of %d features", nrow(data), length(a)))
        )
    },
    nn_network = function(spec, data, s, n_cells) {
        g <- methods::as(s[[spec$to$graph]], "CsparseMatrix")
        a <- Matrix::summary(g)
        b <- Matrix::summary(methods::as(data, "CsparseMatrix"))
        key <- function(i, j, ids) paste(ids[i], ids[j], sep = "\r")
        ka <- key(a$i, a$j, rownames(g))
        kb <- key(b$i, b$j, rownames(data))
        list(
            .rb_row("placement", TRUE, spec$to$graph),
            .rb_row("edges", setequal(ka, kb), sprintf("%d edges", length(ka))),
            .rb_row("weight", .rb_same(a$x[order(ka)], b$x[order(kb)]), "")
        )
    }
)
