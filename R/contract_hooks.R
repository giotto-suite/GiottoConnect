# Hooks are named R functions for conversions the op list cannot express.
# Signature: function(node, params, ctx) -> data for the part. `ctx$add_loss()`
# records losses that depend on the data (e.g. a recovered parameter).

.HOOKS <- list(
    seurat_graph_to_nn = list(
        fn = function(node, params, ctx) .hook_graph_to_nn(node, ctx),
        packages = "Matrix",
        description = "Seurat Graph (sparse adjacency) -> nearest-neighbour edge table; recovers k"
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

# Edges from a (possibly self-including) sparse adjacency. Self edges are
# dropped; for a binary graph k is recovered as the median out-degree.
.hook_graph_to_nn <- function(node, ctx) {
    m <- methods::as(node, "CsparseMatrix")
    ids <- rownames(m)
    if (is.null(ids)) stop(wrap_txt("Graph has no cell names"), call. = FALSE)
    tr <- Matrix::summary(m)
    keep <- tr$i != tr$j
    n_self <- sum(!keep)
    tr <- tr[keep, , drop = FALSE]
    binary <- all(tr$x == 1)
    edges <- data.table::data.table(
        from = ids[tr$i], to = ids[tr$j],
        weight = as.numeric(tr$x), distance = NA_real_
    )
    if (binary) {
        k <- stats::median(tabulate(tr$i, nbins = nrow(m)))
        ctx$add_loss("recovered_from_parameters", "k", sprintf(
            "k = %s recovered from out-degree (%d self edges removed)", k, n_self
        ))
    }
    edges
}
