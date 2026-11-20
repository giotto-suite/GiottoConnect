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
        description = "VisiumV2 spot/bin centers in full-resolution pixels -> microns"
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
# "<type>.<reduction>" and k comes from k.param (which counts the cell itself);
# otherwise the Seurat graph name is kept.
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
    edges <- data.table::data.table(
        from = ids[tr$i], to = ids[tr$j],
        weight = as.numeric(tr$x), distance = NA_real_
    )

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
    .hook_result(edges, to = list(nn_type = nn_type, name = name))
}

# the FindNeighbors command that produced `graph`, or NULL
.seurat_neighbors_record <- function(x, graph) {
    if (is.null(x) || is.null(graph)) return(NULL)
    SO <- asNamespace("SeuratObject")
    cmds <- grep("^FindNeighbors\\.", SO$Command(x), value = TRUE)
    for (cmd in rev(cmds)) {
        gn <- SO$Command(x, command = cmd, value = "graph.name")
        if (graph %in% gn) {
            return(list(command = cmd,
                reduction = SO$Command(x, command = cmd, value = "reduction"),
                k.param = SO$Command(x, command = cmd, value = "k.param")))
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
    xy <- asNamespace("SeuratObject")$GetTissueCoordinates(node)
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
