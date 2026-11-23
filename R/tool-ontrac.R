# ONTraC (Wang et al. 2025, Genome Biology) connection. ONTraC runs its own
# Python code; this file builds its input from a giotto object, runs its
# command line program in an environment made on first use, and writes the
# results to existing slots (adr/0008 in GiottoClass): per-cell NT score and
# niche cluster to cell metadata, per-cell niche cluster probabilities to a
# spatial enrichment object whose misc holds the niche-level NT scores and
# the niche cluster connectivity.

#' @title ONTraC parameters
#' @name ontracParam
#' @description
#' Parameters for [GiottoClass::analyzeData()] that run ONTraC (Wang et al.
#' 2025), which builds niche trajectories with a graph neural network on
#' cell-type composition of each cell's neighbourhood. ONTraC is a Python
#' package; it is installed in a cached environment the first time it runs
#' (through `reticulate::py_require()`, needs reticulate >= 1.41).
#' **If you use it, please cite the original paper** (see References).
#' @param k integer, number of niche clusters (at most 8)
#' @param n_neighbors integer, neighbours for the niche network
#' @param n_local integer, index of the local neighbour used to normalize
#'   distances
#' @param epochs,patience,min_epochs integer, training length
#' @param lr numeric, learning rate
#' @param hidden_feats,n_gcn_layers integer, model size
#' @param modularity_loss_weight,purity_loss_weight,regularization_loss_weight,beta
#'   numeric, loss weights (ONTraC's defaults)
#' @param seed integer, seed
#' @param device `"cpu"` or `"cuda"`
#' @param n_cpu integer, CPUs for building the niche network
#' @param version ONTraC version installed on first use
#' @param out_dir directory for ONTraC's files; `NULL` uses a temporary one
#' @param ... further ONTraC command line options as `name = value`, e.g.
#'   `"trajectory-construct" = "BF"`
#' @returns An `ontracParam` object.
#' @references Wang W, Zheng S, Shin SC, Chávez-Fuentes JC, Yuan GC. ONTraC
#'   characterizes spatially continuous variations of tissue microenvironment
#'   through niche trajectory analysis. *Genome Biology* 26, 117 (2025).
#'   \doi{10.1186/s13059-025-03588-5}
#' @seealso [runONTraC()]
#' @examples
#' p <- ontracParam(k = 6, epochs = 200)
#' p$k
#' @md
#' @export
ontracParam <- function(k = 6L, n_neighbors = 50L, n_local = 20L,
    epochs = 1000L, patience = 100L, min_epochs = 50L, lr = 0.03,
    hidden_feats = 4L, n_gcn_layers = 2L, modularity_loss_weight = 1,
    purity_loss_weight = 30, regularization_loss_weight = 0.1, beta = 0.3,
    seed = 42L, device = c("cpu", "cuda"), n_cpu = 4L, version = "2.1.2",
    out_dir = NULL, ...) {
    checkmate::assert_int(k, lower = 2, upper = 8)
    checkmate::assert_count(n_neighbors, positive = TRUE)
    checkmate::assert_count(n_local, positive = TRUE)
    device <- match.arg(device)
    checkmate::assert_string(version)
    checkmate::assert_string(out_dir, null.ok = TRUE)
    p <- methods::new("ontracParam", param = list(...))
    p$k <- as.integer(k)
    p$n_neighbors <- as.integer(n_neighbors)
    p$n_local <- as.integer(n_local)
    p$epochs <- as.integer(epochs)
    p$patience <- as.integer(patience)
    p$min_epochs <- as.integer(min_epochs)
    p$lr <- lr
    p$hidden_feats <- as.integer(hidden_feats)
    p$n_gcn_layers <- as.integer(n_gcn_layers)
    p$modularity_loss_weight <- modularity_loss_weight
    p$purity_loss_weight <- purity_loss_weight
    p$regularization_loss_weight <- regularization_loss_weight
    p$beta <- beta
    p$seed <- as.integer(seed)
    p$device <- device
    p$n_cpu <- as.integer(n_cpu)
    p$version <- version
    p$out_dir <- out_dir
    p
}

#' @rdname ontracParam
#' @exportClass ontracParam
setClass("ontracParam", contains = "analyzeParam")

#' @rdname ontracParam
#' @param x `data.frame` with columns `cell_ID`, `sdimx`, `sdimy`,
#'   `cell_type` and optionally `sample`
#' @param param an `ontracParam`
#' @details `analyzeData(x, param)` returns a `data.table` with `cell_ID`,
#'   `NTScore`, `niche_cluster` (the most likely cluster) and one probability
#'   column per niche cluster (`NicheCluster_0`, ...). Its attributes
#'   `niche_clusters` (each cluster's NT score), `connectivity` (cluster x
#'   cluster) and `out_dir` hold the niche-level results.
#' @export
setMethod("analyzeData", signature(x = "data.frame", param = "ontracParam"),
    function(x, param, ...) {
        need <- c("cell_ID", "sdimx", "sdimy", "cell_type")
        miss <- setdiff(need, names(x))
        if (length(miss)) {
            stop("[ontracParam] `x` lacks column(s): ",
                paste(miss, collapse = ", "), call. = FALSE)
        }
        exe <- .tool_python_cli(
            requirement = paste0("ONTraC[analysis]==", param$version),
            module = "ONTraC", cli = "ONTraC",
            python_version = ">=3.10,<3.13"
        )
        out <- param$out_dir %null% tempfile("ontrac_")
        dir.create(out, recursive = TRUE, showWarnings = FALSE)
        sample <- if ("sample" %in% names(x)) x$sample else "sample_1"
        meta <- data.table::data.table(Cell_ID = x$cell_ID,
            Sample = as.character(sample), x = x$sdimx, y = x$sdimy,
            Cell_Type = as.character(x$cell_type))
        meta_file <- file.path(out, "meta.csv")
        data.table::fwrite(meta, meta_file)

        args <- c("--meta-input", meta_file,
            "--NN-dir", file.path(out, "NN"),
            "--GNN-dir", file.path(out, "GNN"),
            "--NT-dir", file.path(out, "NT"),
            .ontrac_args(param))
        .tool_run(exe, args, file.path(out, "ontrac.log"))
        .ontrac_read(out)
    }
)

# command line options from the param
.ontrac_args <- function(p) {
    opts <- list(
        "n-cpu" = p$n_cpu, "n-neighbors" = p$n_neighbors,
        "n-local" = p$n_local, "device" = p$device, "epochs" = p$epochs,
        "patience" = p$patience, "min-epochs" = p$min_epochs,
        "seed" = p$seed, "lr" = p$lr, "hidden-feats" = p$hidden_feats,
        "n-gcn-layers" = p$n_gcn_layers, "k-clusters" = p$k,
        "modularity-loss-weight" = p$modularity_loss_weight,
        "purity-loss-weight" = p$purity_loss_weight,
        "regularization-loss-weight" = p$regularization_loss_weight,
        "beta" = p$beta
    )
    known <- c(names(formals(ontracParam)), "")
    extra <- p@param[setdiff(names(p@param), known)]
    opts <- c(opts, extra)
    unlist(Map(function(n, v) c(paste0("--", n), as.character(v)),
        names(opts), opts), use.names = FALSE)
}

# ONTraC's output files -> per-cell table with niche-level attributes
.ontrac_read <- function(out) {
    rd <- function(..., header = TRUE) {
        data.table::fread(file.path(out, ...), header = header)
    }
    nt <- rd("NT", "NTScore.csv.gz")
    prob <- rd("GNN", "cell_level_niche_cluster.csv.gz")
    nc_score <- rd("NT", "niche_cluster_score.csv.gz", header = FALSE)
    adj <- as.matrix(rd("GNN", "consolidate_out_adj.csv.gz", header = FALSE))
    clusters <- setdiff(names(prob), "Cell_ID")
    dimnames(adj) <- list(clusters, clusters)

    res <- data.table::data.table(cell_ID = as.character(prob$Cell_ID))
    res$NTScore <- nt$Cell_NTScore[match(res$cell_ID, as.character(nt$Cell_ID))]
    pm <- as.matrix(prob[, clusters, with = FALSE])
    res$niche_cluster <- clusters[max.col(pm, ties.method = "first")]
    res <- cbind(res, prob[, clusters, with = FALSE])
    attr(res, "niche_clusters") <- data.table::data.table(
        niche_cluster = clusters, NTScore = nc_score[[1]])
    attr(res, "connectivity") <- adj
    attr(res, "out_dir") <- out
    res
}

#' @title Run ONTraC niche trajectories
#' @name runONTraC
#' @description
#' Runs ONTraC (Wang et al. 2025) on the cell types and spatial locations of
#' a giotto object and stores its results:
#' * cell metadata `<name>_NTScore` (niche trajectory score of each cell)
#'   and `<name>_niche_cluster` (most likely niche cluster, ordered by the
#'   clusters' NT scores);
#' * a spatial enrichment `<name>` with each cell's niche cluster
#'   probabilities. Its `misc` holds `niche_clusters` (each cluster's NT
#'   score) and `connectivity` (cluster x cluster).
#'
#' ONTraC is a Python package, installed in a cached environment the first
#' time this runs. **If you use it, please cite the original paper** (see
#' [ontracParam()]).
#' @param gobject giotto object
#' @param cell_type cell metadata column with the cell types
#' @param spat_unit,feat_type spatial unit and feature type
#' @param spat_loc_name name of the spatial locations
#' @param sample_column optional cell metadata column with the sample of
#'   each cell; one sample when `NULL`
#' @param param an [ontracParam()]; its defaults when `NULL`
#' @param name prefix of the results
#' @param verbose be verbose
#' @returns The giotto object with the results.
#' @inherit ontracParam references
#' @seealso [ontracParam()]
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("vizgen")
#' g <- runONTraC(g, cell_type = "leiden_clus", spat_unit = "aggregate",
#'     param = ontracParam(n_neighbors = 30, n_local = 10))
#' }
#' @md
#' @export
runONTraC <- function(gobject, cell_type, spat_unit = NULL, feat_type = NULL,
    spat_loc_name = "raw", sample_column = NULL, param = NULL,
    name = "ontrac", verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_string(cell_type)
    checkmate::assert_string(name)
    param <- param %null% ontracParam()
    checkmate::assert_class(param, "ontracParam")
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject, spat_unit = spat_unit,
        feat_type = feat_type)

    meta <- pDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
    if (!cell_type %in% names(meta)) {
        stop("[runONTraC] no cell metadata column '", cell_type, "'",
            call. = FALSE)
    }
    sl <- getSpatialLocations(gobject, spat_unit = spat_unit,
        name = spat_loc_name, output = "data.table")
    x <- merge(sl[, c("cell_ID", "sdimx", "sdimy")],
        data.table::data.table(cell_ID = meta$cell_ID,
            cell_type = meta[[cell_type]],
            sample = if (is.null(sample_column)) "sample_1" else
                meta[[sample_column]]),
        by = "cell_ID")

    vmsg(.v = verbose, "runONTraC: running ONTraC", param$version)
    res <- analyzeData(x, param)
    clusters <- attr(res, "niche_clusters")
    ordered <- clusters$niche_cluster[order(clusters$NTScore)]

    new_meta <- data.table::data.table(cell_ID = res$cell_ID)
    new_meta[[paste0(name, "_NTScore")]] <- res$NTScore
    new_meta[[paste0(name, "_niche_cluster")]] <- factor(res$niche_cluster,
        levels = ordered)
    gobject <- addCellMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, new_metadata = new_meta, by_column = TRUE,
        column_cell_ID = "cell_ID")

    enr <- createSpatEnrObj(
        enrichment_data = res[, c("cell_ID", clusters$niche_cluster),
            with = FALSE],
        name = name, spat_unit = spat_unit, feat_type = feat_type,
        method = "ONTraC",
        misc = list(niche_clusters = clusters,
            connectivity = attr(res, "connectivity"),
            ontrac_version = param$version),
        verbose = FALSE
    )
    gobject <- setSpatialEnrichment(gobject, x = enr, spat_unit = spat_unit,
        feat_type = feat_type, name = name, verbose = FALSE)
    update_giotto_params(gobject, description = "_ontrac",
        attachments = list(param = param))
}
