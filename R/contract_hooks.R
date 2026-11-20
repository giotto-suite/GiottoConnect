# Hooks are named R functions for conversions the op list cannot express.
# Signature: function(node, params, ctx) -> data for the part, or
# .hook_result(data, to = list(...)) to also set target keys that depend on
# the data. ctx$x is the source object; ctx$add_loss() records losses that
# depend on the data (e.g. a recovered parameter).

.HOOKS <- list(
    seurat_graph_to_nn = list(
        fn = function(node, params, ctx) .hook_graph_to_nn(node, params, ctx),
        packages = c("Matrix", "SeuratObject"),
        description = "Seurat Graph -> nearest-neighbour edges; type, name and k from the FindNeighbors record"
    ),
    seurat_visium_to_um = list(
        fn = function(node, params, ctx) .hook_visium_to_um(node, params, ctx),
        packages = c("Seurat", "SeuratObject"),
        description = "VisiumV1/V2 spot/bin centers in full-resolution pixels -> microns"
    ),
    giotto_polygons_to_vertices = list(
        fn = function(node, params, ctx) .hook_polygons_to_vertices(node),
        packages = c("GiottoClass", "terra"),
        description = "giottoPolygon (in memory or GiottoDisk) -> vertex table x, y, cell"
    ),
    giotto_nn_to_graph = list(
        fn = function(node, params, ctx) .hook_nn_to_graph(node, params, ctx),
        packages = c("GiottoClass", "Matrix"),
        description = "Giotto nearest-neighbour edges -> cells x cells adjacency for a Seurat Graph"
    )
)

#' @title Registered mapping hooks
#' @name hookRegistry
#' @description List the hooks a mapping can reference, the packages each
#' needs and whether those packages are installed.
#' @returns data.frame with columns `hook`, `packages`, `available`,
#' `description`
#' @examples
#' hookRegistry()
#' @export
hookRegistry <- function() {
    data.frame(
        hook = names(.HOOKS),
        packages = vapply(.HOOKS, function(h) paste(h$packages, collapse = ", "), character(1)),
        available = vapply(.HOOKS, function(h) {
            all(vapply(h$packages, requireNamespace, logical(1), quietly = TRUE))
        }, logical(1)),
        description = vapply(.HOOKS, `[[`, character(1), "description"),
        row.names = NULL, stringsAsFactors = FALSE
    )
}

.hook_get <- function(name) {
    h <- .HOOKS[[name]]
    if (is.null(h)) stop(wrap_txt(sprintf("Unknown hook '%s'", name)), call. = FALSE)
    for (p in h$packages) package_check(p)
    h$fn
}

.hook_result <- function(data, to = list()) {
    structure(list(data = data, to = to), class = "giottoHookResult")
}

# Edges from a (possibly self-including) sparse adjacency; self edges are
# dropped. A binary graph is a kNN, a weighted one an sNN. When Seurat's
# FindNeighbors record lists this graph, the network is named
# "<type>.<reduction>", k comes from k.param (which counts the cell itself)
# and kNN distances are recomputed in the recorded reduction, dims and metric;
# weights then follow Giotto's 1 / (1 + distance). Without a record the Seurat
# graph name is kept and distances/weights are NA.
.hook_graph_to_nn <- function(node, params, ctx) {
    m <- methods::as(node, "CsparseMatrix")
    ids <- rownames(m)
    if (is.null(ids)) stop(wrap_txt("Graph has no cell names"), call. = FALSE)
    tr <- Matrix::summary(m)
    keep <- tr$i != tr$j
    n_self <- sum(!keep)
    tr <- tr[keep, , drop = FALSE]
    binary <- all(tr$x == 1)
    nn_type <- if (binary) "kNN" else "sNN"
    edges <- data.table::data.table(from = ids[tr$i], to = ids[tr$j],
        weight = NA_real_, distance = NA_real_)

    graph <- params$nn_graph
    rec <- .seurat_neighbors_record(ctx$x, graph)
    if (!is.null(rec)) {
        name <- paste(nn_type, rec$reduction, sep = ".")
        k_note <- if (binary && !is.null(rec$k.param)) {
            sprintf("; k = %d neighbours (k.param %d includes the cell itself)",
                as.integer(rec$k.param) - 1L, as.integer(rec$k.param))
        } else ""
        ctx$add_loss("recovered_from_parameters", "name", sprintf(
            "%s built by %s on '%s'%s", graph, rec$command, rec$reduction, k_note))
    } else {
        name <- graph
        ctx$add_loss("ambiguous", "name", sprintf(
            "no FindNeighbors record for %s; its reduction is unknown, Seurat name kept", graph))
        if (binary) {
            k <- stats::median(tabulate(tr$i, nbins = nrow(m)))
            ctx$add_loss("recovered_from_parameters", "k", sprintf(
                "k = %s recovered from out-degree (%d self edges removed)", k, n_self))
        }
    }

    if (binary && !is.null(rec)) {
        d <- .neighbor_distances(ctx$x, rec, edges$from, edges$to)
        edges$distance <- d$distance
        edges$weight <- 1 / (1 + d$distance)
        ctx$add_loss("recovered_from_parameters", "distance", sprintf(
            "recomputed in %s dims %s (%s%s); weight = 1 / (1 + distance) replaces the binary adjacency",
            rec$reduction, d$dims, d$metric, if (isTRUE(d$l2)) ", L2-normalised" else ""))
    } else {
        ctx$add_loss("dropped", "distance", if (binary) {
            "Seurat graphs store adjacency only and no FindNeighbors record gives the space; distances are NA"
        } else {
            "sNN graphs carry shared-neighbour weights, not distances; distances are NA"
        })
        ctx$add_loss("dropped", "weight", "Giotto derives weights from distances; weights are NA")
    }
    .hook_result(edges, to = list(nn_type = nn_type, name = name))
}

# distances between neighbour pairs in the space FindNeighbors used
.neighbor_distances <- function(x, rec, from, to) {
    SO <- asNamespace("SeuratObject")
    emb <- SO$Embeddings(x, reduction = rec$reduction)
    dims <- rec$dims %null% seq_len(ncol(emb))
    emb <- emb[, dims, drop = FALSE]
    l2 <- isTRUE(rec$l2.norm)
    if (l2) emb <- emb / sqrt(rowSums(emb^2))
    metric <- rec$annoy.metric %null% "euclidean"
    a <- emb[from, , drop = FALSE]
    b <- emb[to, , drop = FALSE]
    distance <- switch(metric,
        euclidean = sqrt(rowSums((a - b)^2)),
        cosine = sqrt(pmax(0, 2 - 2 * rowSums(a * b) / sqrt(rowSums(a^2) * rowSums(b^2)))),
        manhattan = rowSums(abs(a - b)),
        stop(wrap_txt(sprintf("Unsupported FindNeighbors metric '%s'", metric)), call. = FALSE)
    )
    rng <- range(dims)
    list(distance = unname(distance), metric = metric, l2 = l2,
        dims = if (identical(dims, seq(rng[1], rng[2]))) sprintf("%d:%d", rng[1], rng[2]) else
            paste(dims, collapse = ","))
}

# the FindNeighbors command that produced `graph`, or NULL
.seurat_neighbors_record <- function(x, graph) {
    if (is.null(x) || is.null(graph)) return(NULL)
    SO <- asNamespace("SeuratObject")
    cmds <- grep("^FindNeighbors\\.", SO$Command(x), value = TRUE)
    for (cmd in rev(cmds)) {
        prm <- methods::slot(SO$Command(x, command = cmd), "params")
        if (graph %in% prm$graph.name) {
            return(list(command = cmd, reduction = prm$reduction, k.param = prm$k.param,
                dims = prm$dims, annoy.metric = prm$annoy.metric, l2.norm = prm$l2.norm))
        }
    }
    NULL
}

# VisiumV2 centers are full-resolution pixels. Visium HD bins: the image's
# "spot" scale factor is the bin size in pixels, and the bin size in microns
# is in the image/assay name Load10X_Spatial gives it ("slice1.016um").
# Standard Visium: spots sit on a 100 um center-to-center grid, so the pixel
# size comes from the measured nearest-neighbour spacing.
.hook_visium_to_um <- function(node, params, ctx) {
    if (methods::is(node, "VisiumV1")) {
        # VisiumV1 keeps full-resolution imagerow (y) / imagecol (x)
        v1 <- asNamespace("SeuratObject")$GetTissueCoordinates(node, scale = NULL)
        xy <- data.frame(x = v1$imagecol, y = v1$imagerow, cell = rownames(v1),
            stringsAsFactors = FALSE)
    } else {
        xy <- asNamespace("SeuratObject")$GetTissueCoordinates(node)
    }
    sf <- methods::slot(node, "scale.factors")
    bin <- regmatches(params$image, regexpr("(?<=\\.)\\d+(?=um$)", params$image, perl = TRUE))
    if (length(bin)) {
        um_per_px <- as.numeric(bin) / sf[["spot"]]
        how <- sprintf("%s um bins / spot scale factor %.4f px", bin, sf[["spot"]])
    } else {
        nn_px <- .median_nn_spacing(as.matrix(xy[, c("x", "y")]))
        um_per_px <- 100 / nn_px
        how <- sprintf("100 um spot pitch / measured spacing %.2f px", nn_px)
    }
    ctx$add_loss("recovered_from_parameters", "units", sprintf(
        "%.5f um per px (%s)", um_per_px, how))
    xy$x <- xy$x * um_per_px
    xy$y <- xy$y * um_per_px
    xy
}

.median_nn_spacing <- function(m, n = 500L) {
    if (nrow(m) > n) m <- m[as.integer(seq(1L, nrow(m), length.out = n)), , drop = FALSE]
    d <- as.matrix(stats::dist(m))
    diag(d) <- Inf
    stats::median(apply(d, 1L, min))
}

# outer ring vertices of each polygon (holes are dropped), in Giotto's
# y-up convention; the mapping's ops flip y back
.hook_polygons_to_vertices <- function(node) {
    sv <- if (inherits(node[], "SpatVector")) node[] else GiottoClass::as.terra(node[])
    geo <- as.data.frame(terra::geom(sv))
    geo <- geo[geo$hole == 0 & geo$part == 1, , drop = FALSE]
    data.frame(x = geo$x, y = geo$y, cell = sv$poly_ID[geo$geom], stringsAsFactors = FALSE)
}

# Seurat kNN graphs are binary and count each cell as its own neighbour;
# sNN graphs carry the edge weights
.hook_nn_to_graph <- function(node, params, ctx) {
    cells <- GiottoClass::spatIDs(ctx$x, spat_unit = params$spat_unit)
    e <- as.data.frame(node)
    i <- match(e$from, cells)
    j <- match(e$to, cells)
    knn <- grepl("_nn$", ctx$part_spec$to$graph %null% "")
    if (knn) {
        i <- c(i, seq_along(cells))
        j <- c(j, seq_along(cells))
        x <- rep(1, length(i))
        ctx$add_loss("estimated", "weight", "kNN graph written as binary adjacency with self edges (Seurat convention); distances dropped")
    } else {
        x <- e$weight
        ctx$add_loss("ambiguous", "weight", "Giotto sNN weights are not Seurat's Jaccard shared-neighbour weights")
    }
    m <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(length(cells), length(cells)),
        dimnames = list(cells, cells))
    m@x[m@x > 1 & knn] <- 1
    m
}
