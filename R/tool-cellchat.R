# CellChat (Jin et al. 2021; CellChat v2, Jin et al. 2025) connection.
# CellChat runs its own R code; this file builds its input from a giotto
# object (normalized expression, a cell type column and, when present, the
# spatial locations) and runs CellChat's standard workflow (adr/0008 in
# GiottoClass).
#
# Result -> slot: CellChat's result is a cell type x cell type x
# ligand-receptor table plus a CellChat object. No giotto slot holds a table
# between cell groups (every slot is per cell or per feature), inventing a
# spat_unit or feat_type to park it is not allowed, and storing the CellChat
# object would copy the expression matrix into the giotto object. So
# runCellChat() leaves the giotto object unchanged and *returns* the
# interaction table, like Giotto's own spatCellCellcom(); the pathway table,
# the CellChat object (for CellChat's plots) and the param are attributes.

#' @title CellChat parameters
#' @name cellchatParam
#' @description
#' Parameters for [GiottoClass::analyzeData()] that run CellChat (Jin et al.
#' 2021) and its spatial mode from CellChat v2 (Jin et al. 2025), which infers
#' cell-cell communication between cell groups from the expression of
#' ligands, receptors and their cofactors (CellChatDB), with spatial
#' distance between the groups as a constraint when coordinates are given.
#' CellChat is an R package from GitHub, install it with
#' `remotes::install_github("jinworks/CellChat")`.
#' **If you use it, please cite the original paper** (see References).
#' @param species CellChatDB to use: `"mouse"`, `"human"` or `"zebrafish"`
#' @param db_search CellChatDB annotations to keep (`subsetDB()`); `NULL`
#'   keeps all protein signaling ("Secreted Signaling", "ECM-Receptor",
#'   "Cell-Cell Contact"), CellChat's recommendation
#' @param non_protein logical, also keep "Non-protein Signaling"
#' @param type,trim how average expression per group is computed (CellChat's
#'   `computeCommunProb()`; `"truncatedMean"` with `trim = 0.1` is less strict
#'   for low-depth spots)
#' @param population_size logical, weigh by group sizes
#' @param variable_both logical, require both ligand and receptor to be
#'   over-expressed (`identifyOverExpressedInteractions()`)
#' @param thresh_p p-value threshold for over-expressed genes
#' @param do_fast logical, CellChat's fast Wilcoxon test, which needs the
#'   GitHub package presto; `NULL` uses it when presto is installed
#' @param nboot integer, permutations for the p-values
#' @param seed integer, seed of the permutations
#' @param min_cells integer, minimum cells per group (`filterCommunication()`)
#' @param thresh p-value threshold of the returned interactions
#' @param distance_use logical, spatial mode: use distance as a constraint
#' @param interaction_range numeric, spatial mode: maximum interaction length
#'   in microns
#' @param scale_distance numeric, spatial mode: distance scale factor
#' @param k_min integer, spatial mode: minimum number of interacting cell
#'   pairs between two groups
#' @param contact_dependent logical, spatial mode: run "Cell-Cell Contact"
#'   pairs only between touching groups
#' @param contact_range numeric, spatial mode: contact distance in microns
#'   (about the cell diameter, or the spot distance for Visium)
#' @param ratio microns per coordinate unit of the spatial locations; `NULL`
#'   derives it from `center_distance`
#' @param center_distance microns between neighbouring centres (e.g. 100 for
#'   Visium spots); with `ratio = NULL` the ratio is `center_distance` divided
#'   by the median nearest-neighbour distance of the coordinates. One of
#'   `ratio` or `center_distance` is needed for the spatial mode
#' @param tol tolerance in microns added to the ranges (CellChat uses half
#'   the spot or cell size, e.g. 32.5 for Visium); `NULL` uses half the
#'   median nearest-neighbour distance in microns
#' @param ... further options of `computeCommunProb()` as `name = value`
#' @returns A `cellchatParam` object.
#' @references Jin S, Guerrero-Juarez CF, Zhang L, Chang I, Ramos R, Kuan CH,
#'   Myung P, Plikus MV, Nie Q. Inference and analysis of cell-cell
#'   communication using CellChat. *Nature Communications* 12, 1088 (2021).
#'   \doi{10.1038/s41467-021-21246-9}
#'
#'   Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell-cell
#'   communication from single-cell transcriptomics. *Nature Protocols* 20,
#'   180-219 (2025). \doi{10.1038/s41596-024-01045-4}
#' @seealso [runCellChat()]
#' @examples
#' p <- cellchatParam(species = "mouse", center_distance = 100, tol = 32.5)
#' p$species
#' @md
#' @export
cellchatParam <- function(species = c("mouse", "human", "zebrafish"),
    db_search = NULL, non_protein = FALSE,
    type = c("triMean", "truncatedMean", "thresholdedMean", "median"),
    trim = 0.1, population_size = FALSE, variable_both = TRUE,
    thresh_p = 0.05, do_fast = NULL, nboot = 100L, seed = 1L, min_cells = 10L, thresh = 0.05,
    distance_use = TRUE, interaction_range = 250, scale_distance = 0.01,
    k_min = 10L, contact_dependent = TRUE, contact_range = 100,
    ratio = NULL, center_distance = NULL, tol = NULL, ...) {
    species <- match.arg(species)
    type <- match.arg(type)
    checkmate::assert_character(db_search, null.ok = TRUE)
    checkmate::assert_flag(non_protein)
    checkmate::assert_number(trim, lower = 0, upper = 0.5)
    checkmate::assert_flag(population_size)
    checkmate::assert_flag(variable_both)
    checkmate::assert_flag(do_fast, null.ok = TRUE)
    checkmate::assert_count(nboot, positive = TRUE)
    checkmate::assert_int(seed)
    checkmate::assert_count(min_cells)
    checkmate::assert_number(thresh, lower = 0, upper = 1)
    checkmate::assert_flag(distance_use)
    checkmate::assert_number(interaction_range, lower = 0)
    checkmate::assert_number(contact_range, lower = 0, null.ok = TRUE)
    checkmate::assert_number(ratio, lower = 0, null.ok = TRUE)
    checkmate::assert_number(center_distance, lower = 0, null.ok = TRUE)
    checkmate::assert_number(tol, lower = 0, null.ok = TRUE)
    p <- methods::new("cellchatParam", param = list(...))
    p$species <- species
    p$db_search <- db_search
    p$non_protein <- non_protein
    p$type <- type
    p$trim <- trim
    p$population_size <- population_size
    p$variable_both <- variable_both
    p$thresh_p <- thresh_p
    p$do_fast <- do_fast
    p$nboot <- as.integer(nboot)
    p$seed <- as.integer(seed)
    p$min_cells <- as.integer(min_cells)
    p$thresh <- thresh
    p$distance_use <- distance_use
    p$interaction_range <- interaction_range
    p$scale_distance <- scale_distance
    p$k_min <- as.integer(k_min)
    p$contact_dependent <- contact_dependent
    p$contact_range <- contact_range
    p$ratio <- ratio
    p$center_distance <- center_distance
    p$tol <- tol
    p
}

#' @rdname cellchatParam
#' @exportClass cellchatParam
setClass("cellchatParam", contains = "analyzeParam")

# the paper, shown once per session and listed by objCitations()
setMethod("methodCitation", signature("cellchatParam"), function(x, ...) {
    c(
        .cite_article(
            author = c("Suoqin Jin", "Christian F. Guerrero-Juarez",
                "Lihua Zhang", "others"),
            title = paste("Inference and analysis of cell-cell communication",
                "using CellChat"),
            journal = "Nature Communications", year = 2021, volume = 12,
            pages = "1088", doi = "10.1038/s41467-021-21246-9"
        ),
        .cite_article(
            author = c("Suoqin Jin", "Maksim V. Plikus", "Qing Nie"),
            title = paste("CellChat for systematic analysis of cell-cell",
                "communication from single-cell transcriptomics"),
            journal = "Nature Protocols", year = 2025, volume = 20,
            pages = "180--219", doi = "10.1038/s41596-024-01045-4"
        )
    )
})

#' @rdname cellchatParam
#' @param x normalized (log) expression matrix, features x cells, feature
#'   names as gene symbols of `species`
#' @param param a `cellchatParam`
#' @param cell_type cell group of each column of `x`
#' @param spatial_locs `NULL` (non-spatial mode) or a two-column matrix of
#'   coordinates, rows in the column order of `x`
#' @details `analyzeData(x, param, cell_type, spatial_locs)` returns a
#'   `data.table` with one row per significant interaction (`source`,
#'   `target`, `ligand`, `receptor`, `prob`, `pval`, `interaction_name`,
#'   `pathway_name`, `annotation`, ...; CellChat's `subsetCommunication()`).
#'   Its attributes are `pathways` (the pathway-level table), `cellchat`
#'   (the CellChat object, for CellChat's plotting functions) and
#'   `spatial_factors` (the `ratio` and `tol` used; `NULL` when non-spatial).
#' @export
setMethod("analyzeData", signature(x = "allMatrix", param = "cellchatParam"),
    function(x, param, cell_type, spatial_locs = NULL, ...) {
        package_check("CellChat", repository = "github:jinworks/CellChat")
        if (length(cell_type) != ncol(x)) {
            stop("[cellchatParam] `cell_type` needs one value per column of `x`",
                call. = FALSE)
        }
        if (is.null(colnames(x))) colnames(x) <- paste0("cell_", seq_len(ncol(x)))
        meta <- data.frame(labels = droplevels(as.factor(cell_type)),
            samples = factor("sample1"), row.names = colnames(x))
        sf <- NULL
        if (is.null(spatial_locs)) {
            cc <- CellChat::createCellChat(x, meta = meta, group.by = "labels",
                datatype = "RNA")
        } else {
            spatial_locs <- as.matrix(spatial_locs)[, 1:2, drop = FALSE]
            dimnames(spatial_locs) <- list(colnames(x), c("x", "y"))
            sf <- .cellchat_spatial_factors(spatial_locs, param)
            cc <- CellChat::createCellChat(x, meta = meta, group.by = "labels",
                datatype = "spatial", coordinates = spatial_locs,
                spatial.factors = sf)
        }
        cc <- .cellchat_workflow(cc, param, spatial = !is.null(sf))
        res <- data.table::as.data.table(CellChat::subsetCommunication(cc,
            slot.name = "net", thresh = param$thresh))
        attr(res, "pathways") <- data.table::as.data.table(
            CellChat::subsetCommunication(cc, slot.name = "netP",
                thresh = param$thresh))
        attr(res, "cellchat") <- cc
        attr(res, "spatial_factors") <- sf
        res
    }
)

# CellChat's ratio (microns per coordinate unit) and tol (microns)
.cellchat_spatial_factors <- function(locs, p) {
    nn <- NULL
    nn_dist <- function() {
        if (is.null(nn)) {
            net <- createNetwork(locs, kNNNetworkParam(k = 1L,
                include_weight = FALSE, output = "data.table"))
            nn <<- stats::median(net$distance)
        }
        nn
    }
    ratio <- p$ratio
    if (is.null(ratio)) {
        if (is.null(p$center_distance)) {
            stop(wrap_txt("[cellchatParam] the spatial mode needs the micron",
                "scale of the coordinates: set `ratio` (microns per unit, 1",
                "when the coordinates are in microns) or `center_distance`",
                "(e.g. 100 for Visium spots)"), call. = FALSE)
        }
        ratio <- p$center_distance / nn_dist()
    }
    tol <- p$tol %null% (ratio * nn_dist() / 2)
    data.frame(ratio = ratio, tol = tol)
}

# CellChat's standard workflow (CellChat v2 tutorials) on a CellChat object
.cellchat_workflow <- function(cc, p, spatial) {
    db <- .cellchat_db(p$species)
    cc@DB <- CellChat::subsetDB(db, search = p$db_search,
        non_protein = p$non_protein)
    cc <- CellChat::subsetData(cc)
    do_fast <- p$do_fast %null% requireNamespace("presto", quietly = TRUE)
    cc <- CellChat::identifyOverExpressedGenes(cc, thresh.p = p$thresh_p,
        do.fast = do_fast)
    cc <- CellChat::identifyOverExpressedInteractions(cc,
        variable.both = p$variable_both)
    args <- list(cc, type = p$type, trim = p$trim,
        population.size = p$population_size, nboot = p$nboot,
        seed.use = p$seed)
    if (spatial) {
        args <- c(args, list(distance.use = p$distance_use,
            interaction.range = p$interaction_range,
            scale.distance = p$scale_distance, k.min = p$k_min,
            contact.dependent = p$contact_dependent,
            contact.range = p$contact_range))
    }
    known <- c(names(formals(cellchatParam)), "")
    args <- c(args, p@param[setdiff(names(p@param), known)])
    cc <- do.call(CellChat::computeCommunProb, args)
    cc <- CellChat::filterCommunication(cc, min.cells = p$min_cells)
    cc <- CellChat::computeCommunProbPathway(cc, thresh = p$thresh)
    CellChat::aggregateNet(cc, thresh = p$thresh)
}

.cellchat_db <- function(species) {
    e <- new.env()
    nm <- paste0("CellChatDB.", species)
    utils::data(list = nm, package = "CellChat", envir = e)
    e[[nm]]
}

#' @title Run CellChat cell-cell communication
#' @name runCellChat
#' @description
#' Runs CellChat (Jin et al. 2021) on the normalized expression and a cell
#' type column of a giotto object, in CellChat v2's spatial mode (Jin et al.
#' 2025) when `spatial = TRUE`, and returns the significant ligand-receptor
#' interactions between cell types.
#'
#' The giotto object is not changed: the result is a table between cell
#' groups, which no giotto slot holds (as with Giotto's `spatCellCellcom()`
#' it is returned). The spatial mode needs the micron scale of the spatial
#' locations: give `ratio` or `center_distance` in the param (see
#' [cellchatParam()]). **If you use it, please cite the original paper** (see
#' [cellchatParam()]).
#' @param gobject giotto object
#' @param cell_type cell metadata column with the cell types
#' @param spat_unit,feat_type spatial unit and feature type
#' @param expression_values normalized expression values to use
#' @param spat_loc_name name of the spatial locations
#' @param spatial logical, run the spatial mode (needs spatial locations)
#' @param param a [cellchatParam()]; its defaults when `NULL`
#' @param verbose be verbose
#' @returns A `data.table` of significant interactions (see
#'   [cellchatParam()]) with attributes `pathways`, `cellchat` (the CellChat
#'   object, for CellChat's plots), `spatial_factors` and `param`.
#' @inherit cellchatParam references
#' @seealso [cellchatParam()]
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' res <- runCellChat(g, cell_type = "leiden_clus",
#'     param = cellchatParam(species = "mouse", center_distance = 100,
#'         tol = 32.5))
#' CellChat::netVisual_circle(attr(res, "cellchat")@net$count)
#' }
#' @md
#' @export
runCellChat <- function(gobject, cell_type, spat_unit = NULL, feat_type = NULL,
    expression_values = "normalized", spat_loc_name = "raw", spatial = TRUE,
    param = NULL, verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_string(cell_type)
    checkmate::assert_flag(spatial)
    param <- param %null% cellchatParam()
    checkmate::assert_class(param, "cellchatParam")
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject, spat_unit = spat_unit,
        feat_type = feat_type)

    meta <- pDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
    if (!cell_type %in% names(meta)) {
        stop("[runCellChat] no cell metadata column '", cell_type, "'",
            call. = FALSE)
    }
    x <- getExpression(gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = expression_values, output = "matrix")
    labels <- meta[[cell_type]][match(colnames(x), meta$cell_ID)]
    locs <- NULL
    if (spatial) {
        sl <- getSpatialLocations(gobject, spat_unit = spat_unit,
            name = spat_loc_name, output = "data.table")
        sl <- sl[match(colnames(x), sl$cell_ID), ]
        if (anyNA(sl$cell_ID)) {
            stop("[runCellChat] spatial locations miss some cells",
                call. = FALSE)
        }
        locs <- as.matrix(sl[, c("sdimx", "sdimy")])
    }
    vmsg(.v = verbose, "runCellChat: running CellChat",
        as.character(utils::packageVersion("CellChat")),
        if (spatial) "(spatial mode)" else "")
    res <- analyzeData(x, param, cell_type = labels, spatial_locs = locs)
    attr(res, "param") <- param
    citeMethod(param, what = "runCellChat()")
    res
}
