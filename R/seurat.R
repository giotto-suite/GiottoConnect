#' @import GiottoClass
#' @import GiottoUtils
NULL




## Seurat object ####



#' @title Deprecated
#' @name giottoToSeurat
#' @description Deprecated. Please use either [giottoToSeuratV4] or
#' [giottoToSeuratV5]
#' @param gobject Giotto object
#' @param obj_use Giotto object (deprecated, use gobject)
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param ... additional params to pass to \code{\link{getSpatialLocations}}
#' @returns Seurat object
#' @export
giottoToSeurat <- function(
        gobject,
        spat_unit = NULL,
        obj_use = NULL,
        ...) {
    stop(wrap_txt(
        "Deprecated. Please use either giottoToSeuratV4() or giottoToSeuratV5()"
    ))
}




#' @title Convert Giotto to Seurat V4
#' @name giottoToSeuratV4
#' @description Converts Giotto object into a Seurat object. This functions
#' extracts specific sets of data belonging to specified spatial unit.
#' The default values are 'cell' and 'rna' respectively.
#' @param gobject Giotto object
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param ... additional params to pass to \code{\link{getSpatialLocations}}
#' @returns Seurat object
#' @keywords seurat interoperability
#' @export
giottoToSeuratV4 <- function(
        gobject,
        spat_unit = NULL,
        ...) {
    # data.table vars
    feat_type <- name <- dim_type <- nn_type <- NULL
    # set default spat_unit and feat_type to be extracted as a Seurat assay
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    # verify if optional package is installed
    package_check(pkg_name = "Seurat", repository = "CRAN")
    requireNamespace("Seurat")
    # check whether any raw data exist -- required for Seurat
    avail_expr <- list_expression(gobject = gobject, spat_unit = spat_unit)
    raw_exist <- avail_expr[, "raw" %in% name, by = feat_type]
    # raw_exist <- sapply(gobject@expression_feat,function(x)
    #   'raw' %in% names(gobject@expression[[x]]))
    if (nrow(raw_exist) > 0) {
        assays_all <- raw_exist[, feat_type]
        # assays_all <- names(raw_exist[1])
        # assays_all <- union(assays_all,gobject@expression_feat)
    } else {
        stop("Raw count data not found. Required for Seurat object.")
    }
    # create Seurat object when at least one raw data is available
    for (i in seq_along(assays_all)) {
        assay_use <- assays_all[i]
        expr_use <- lapply(
            avail_expr[feat_type == assay_use, name],
            function(x) {
                getExpression(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = assay_use,
                    values = x,
                    output = "exprObj"
                )
            }
        )
        # expr_use <- gobject@expression[[assay_use]]
        names(expr_use) <- unlist(lapply(expr_use, objName))
        slot_use <- names(expr_use)
        if (i == 1) {
            data_raw <- expr_use[["raw"]][]
            sobj <- Seurat::CreateSeuratObject(
                counts = data_raw,
                assay = assay_use
            )
            if ("normalized" %in% slot_use) {
                sobj <- Seurat::SetAssayData(
                    sobj,
                    layer = "data",
                    new.data = expr_use[["normalized"]][],
                    assay = assay_use
                )
            }
            if ("scaled" %in% slot_use) {
                sobj <- Seurat::SetAssayData(
                    sobj,
                    layer = "scale.data",
                    # does not accept 'dgeMatrix'
                    new.data = as.matrix(expr_use[["scaled"]][]),
                    assay = assay_use
                )
            }
        } else {
            if ("raw" %in% slot_use) {
                data_raw <- expr_use[["raw"]][]
                flag_raw <- 1
            } else {
                flag_raw <- 0
            }
            if ("normalized" %in% slot_use) {
                data_norm <- expr_use[["normalized"]][]
                flag_norm <- 1
            } else {
                flag_norm <- 0
            }
            if (flag_raw == 1) {
                assay_obj <- Seurat::CreateAssayObject(counts = data_raw)
            } else if (flag_raw == 0 & flag_norm == 1) {
                assay_obj <- Seurat::CreateAssayObject(data = data_norm)
            } else {
                stop(paste0(
                    "Raw and normalized data not found for assay ",
                    assay_use
                ))
            }
            sobj[[assay_use]] <- assay_obj
            if ("scaled" %in% slot_use) {
                data_scale <- as.matrix(expr_use[["scaled"]][])
                sobj <- Seurat::SetAssayData(
                    sobj,
                    layer = "scale.data",
                    new.data = data_scale,
                    assay = assay_use
                )
            }
        }
        # add cell metadata
        meta_cells <- data.table::setDF(
            getCellMetadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = assay_use,
                output = "data.table",
                copy_obj = TRUE
            )
        )
        rownames(meta_cells) <- meta_cells$cell_ID
        meta_cells <- meta_cells[, -which(colnames(meta_cells) == "cell_ID")]
        if (ncol(meta_cells) > 0) {
            colnames(meta_cells) <- paste0(assay_use, "_", colnames(meta_cells))
        }
        sobj <- Seurat::AddMetaData(sobj,
            metadata = meta_cells[Seurat::Cells(sobj), ]
        )
        # add feature metadata
        meta_genes <- data.table::setDF(
            getFeatureMetadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = assay_use,
                output = "data.table",
                copy_obj = TRUE
            )
        )
        rownames(meta_genes) <- meta_genes$feat_ID
        sobj[[assay_use]]@meta.features <- cbind(
            sobj[[assay_use]]@meta.features, meta_genes
        )
        # dim reduction
        # note: Seurat requires assay name specification for each dim reduc
        avail_dr <- list_dim_reductions(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = assay_use,
        )
        if (!is.null(avail_dr)) {
            if (nrow(avail_dr) > 0) {
                dr_use <- avail_dr[, name]
                for (i in seq(nrow(avail_dr))) {
                    dr_name <- avail_dr[i, name]
                    dr_type <- avail_dr[i, dim_type]
                    dr_obj <- get_dimReduction(
                        gobject = gobject,
                        output = "dimObj",
                        spat_unit = spat_unit,
                        feat_type = assay_use,
                        reduction_method = dr_type,
                        name = dr_name
                    )
                    emb_use <- dr_obj[][Seurat::Cells(sobj), ]
                    if (sum(c("loadings", "eigenvalues") %in%
                        names(slot(dr_obj, "misc"))) == 2) {
                        loadings_use <- slot(dr_obj, "misc")$loadings
                        stdev_use <- slot(dr_obj, "misc")$eigenvalues
                        sobj[[dr_name]] <- Seurat::CreateDimReducObject(
                            embeddings = as.matrix(emb_use),
                            loadings = loadings_use,
                            key = paste0(dr_name, "_"),
                            stdev = stdev_use,
                            assay = assay_use
                        )
                    } else {
                        sobj[[dr_name]] <- Seurat::CreateDimReducObject(
                            embeddings = as.matrix(emb_use),
                            key = paste0(dr_name, "_"),
                            assay = assay_use
                        )
                    }
                }
            }
        }
        # network objects
        # expression network
        avail_nn <- list_nearest_networks(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = assay_use
        )
        if (!is.null(avail_nn)) {
            if (nrow(avail_nn) > 0) {
                for (i in seq(nrow(avail_nn))) {
                    nn_name <- avail_nn[i, name]
                    nn_type <- avail_nn[i, nn_type]
                    nn_use <- get_NearestNetwork(
                        gobject = gobject,
                        spat_unit = spat_unit,
                        feat_type = assay_use,
                        nn_network_to_use = nn_type,
                        network_name = nn_name,
                        output = "data.table"
                    )
                    idx1 <- match(nn_use$from, Seurat::Cells(sobj))
                    idx2 <- match(nn_use$to, Seurat::Cells(sobj))
                    edge_weight <- nn_use$weight
                    edge_dist <- nn_use$distance
                    nn_mtx <- Matrix::sparseMatrix(
                        i = idx1, j = idx2, x = edge_dist,
                        dims = c(ncol(sobj), ncol(sobj))
                    )
                    rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
                    nn_name <- paste0(nn_name)
                    sGraph <- Seurat::as.Graph(nn_mtx)
                    sGraph@assay.used <- assay_use
                    sobj[[nn_name]] <- sGraph
                }
            }
        }
    }
    # spatial coordinates
    loc_use <- data.table::setDF(
        get_spatial_locations(
            gobject = gobject,
            spat_unit = spat_unit,
            output = "data.table",
            copy_obj = TRUE,
            ... # allow setting of spat_loc_name through additional params
        )
    )
    rownames(loc_use) <- loc_use$cell_ID
    sobj <- Seurat::AddMetaData(sobj, metadata = loc_use)
    # add spatial coordinates as new dim reduct object
    loc_2 <- loc_use[, c("sdimx", "sdimy")]
    colnames(loc_2) <- c("spatial_1", "spatial_2")
    sobj[["spatial"]] <- Seurat::CreateDimReducObject(
        embeddings = as.matrix(loc_2),
        assay = names(sobj@assays)[1],
        key = "spatial_"
    )
    # spatial network
    avail_sn <- list_spatial_networks(gobject = gobject, spat_unit = spat_unit)
    if (!is.null(avail_sn)) {
        if (nrow(avail_sn) > 0) {
            sn_all <- avail_sn[, name]
            for (i in sn_all) {
                snt_use <- get_spatialNetwork(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    name = i,
                    output = "networkDT"
                )
                idx1 <- match(snt_use$from, Seurat::Cells(sobj))
                idx2 <- match(snt_use$to, Seurat::Cells(sobj))
                edge_weight <- snt_use$weight
                edge_dist <- snt_use$distance
                nn_mtx <- Matrix::sparseMatrix(
                    i = idx1, j = idx2, x = edge_dist,
                    dims = c(ncol(sobj), ncol(sobj))
                )
                rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
                nn_name <- paste0(i)
                spatGraph <- Seurat::as.Graph(nn_mtx)
                spatGraph@assay.used <- names(sobj@assays)[1]
                sobj[[nn_name]] <- spatGraph
            }
        }
    }
    return(sobj)
}


#' @title Convert Giotto to Seurat V5
#' @name giottoToSeuratV5
#' @description Converts Giotto object into a Seurat object. This functions
#' extracts specific sets of data belonging to specified spatial unit.
#' The default values are 'cell' and 'rna' respectively.
#' @param gobject Giotto object
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param tech technology the dataset is using (e.g. "Visium","Xenium", "Slide-seq")
#' @param res_type type of 10x image output resolution
#' @param ... additional params to pass to \code{\link{getSpatialLocations}}
#' @returns Seurat object
#' @keywords seurat interoperability
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' giottoToSeuratV5(g, tech = "Visium")
#' 
#' @export
giottoToSeuratV5 <- function(gobject,
    spat_unit = NULL,
    tech = c("Visium", "Xenium", "Slide-seq"),
    res_type = c("hires", "lowres", "fullres"),
    ...) {
    # data.table vars
    feat_type <- name <- dim_type <- nn_type <- NULL
    res_type <- match.arg(res_type,
        choices = c("hires", "lowres", "fullres"),
        several.ok = FALSE
    )
    tech <- match.arg(tech,
        choices = c("Visium", "Xenium", "Slide-seq")
    )
    # set default spat_unit and feat_type to be extracted as a Seurat assay
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    assay_names <- names(gobject@expression[[spat_unit]])

    # Identify assays with spaces and replace with underscores
    new_assay_names <- gsub(" ", "_", assay_names)

    # Apply the new names to the gobject expression slot
    names(gobject@expression[[spat_unit]]) <- new_assay_names

    # verify if optional package is installed
    package_check(pkg_name = "Seurat", repository = "CRAN")
    loadNamespace("Seurat")

    # check whether any raw data exist -- required for Seurat
    avail_expr <- list_expression(gobject = gobject, spat_unit = spat_unit)
    raw_exist <- avail_expr[, "raw" %in% name, by = feat_type]
    # raw_exist <- sapply(gobject@expression_feat,function(x)
    # 'raw' %in% names(gobject@expression[[x]]))
    if (nrow(raw_exist) > 0) {
        assays_all <- raw_exist[, feat_type]
        # assays_all <- names(raw_exist[1])
        # assays_all <- union(assays_all,gobject@expression_feat)
    } else {
        stop("Raw count data not found. Required for Seurat object.")
    }

    # create Seurat object when at least one raw data is available
    for (i in seq_along(assays_all)) {
        assay_use <- assays_all[i]
        expr_use <- lapply(
            avail_expr[feat_type == assay_use, name],
            function(x) {
                getExpression(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = assay_use,
                    values = x,
                    output = "exprObj"
                )
            }
        )
        # expr_use <- gobject@expression[[assay_use]]
        names(expr_use) <- unlist(lapply(expr_use, objName))
        slot_use <- names(expr_use)
        if (i == 1) {
            data_raw <- expr_use[["raw"]][]
            sobj <- Seurat::CreateSeuratObject(
                counts = data_raw,
                assay = assay_use
            )
            if ("normalized" %in% slot_use) {
                sobj <- Seurat::SetAssayData(
                    sobj,
                    layer = "data",
                    new.data = expr_use[["normalized"]][],
                    assay = assay_use
                )
            }
            if ("scaled" %in% slot_use) {
                sobj <- Seurat::SetAssayData(
                    sobj,
                    layer = "scale.data",
                    # does not accept 'dgeMatrix'
                    new.data = as.matrix(expr_use[["scaled"]][]),
                    assay = assay_use
                )
            }
        } else {
            if ("raw" %in% slot_use) {
                data_raw <- expr_use[["raw"]][]
                flag_raw <- 1
            } else {
                flag_raw <- 0
            }
            if ("normalized" %in% slot_use) {
                data_norm <- expr_use[["normalized"]][]
                flag_norm <- 1
            } else {
                flag_norm <- 0
            }
            if (flag_raw == 1) {
                assay_obj <- Seurat::CreateAssayObject(counts = data_raw)
            } else if (flag_raw == 0 & flag_norm == 1) {
                assay_obj <- Seurat::CreateAssayObject(data = data_norm)
            } else {
                stop(paste0(
                    "Raw and normalized data not found for assay ",
                    assay_use
                ))
            }
            sobj[[assay_use]] <- assay_obj
            if ("scaled" %in% slot_use) {
                data_scale <- as.matrix(expr_use[["scaled"]][])
                sobj <- Seurat::SetAssayData(
                    sobj,
                    layer = "scale.data",
                    new.data = data_scale,
                    assay = assay_use
                )
            }
        }

        # add cell metadata
        names(gobject@cell_metadata$cell) <- gsub(
            " ", "_", names(gobject@cell_metadata$cell))
        meta_cells <- data.table::setDF(
            getCellMetadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = assay_use,
                output = "data.table",
                copy_obj = TRUE
            )
        )
        rownames(meta_cells) <- meta_cells$cell_ID
        meta_cells <- meta_cells[
            , -which(colnames(meta_cells) == "cell_ID"), drop = FALSE]
        if (ncol(meta_cells) > 0) {
            colnames(meta_cells) <- paste0(
                assay_use, "_",
                colnames(meta_cells)
            )
            sobj <- Seurat::AddMetaData(sobj,
                metadata = meta_cells[Seurat::Cells(sobj), ],
                col.name = names(meta_cells)
            )
        }

        # add feature metadata
        names(gobject@feat_metadata$cell) <- gsub(
            " ", "_", names(gobject@feat_metadata$cell))
        meta_genes <- data.table::setDF(
            getFeatureMetadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = assay_use,
                output = "data.table",
                copy_obj = TRUE
            )
        )
        rownames(meta_genes) <- meta_genes$feat_ID
        if ("meta.data" %in% slotNames(sobj@assays[[assay_use]])) {
            sobj@assays[[assay_use]]@meta.data <- meta_genes
            message(paste("Meta data updated for assay:", assay_use))
        } else if ("meta.features" %in% slotNames(sobj@assays[[assay_use]])) {
            sobj@assays[[assay_use]]@meta.features <- meta_genes
            message(paste("Meta features updated for assay:", assay_use))
        } else {
            warning(paste("No suitable metadata slot found for assay", assay_use))
        }

        # dim reduction
        # note: Seurat requires assay name specification for each dim reduc
        avail_dr <- list_dim_reductions(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = assay_use,
        )
        if (!is.null(avail_dr)) {
            if (nrow(avail_dr) > 0) {
                dr_use <- avail_dr[, name]
                for (i in seq(nrow(avail_dr))) {
                    dr_name <- avail_dr[i, name]
                    dr_type <- avail_dr[i, dim_type]
                    dr_obj <- get_dimReduction(
                        gobject = gobject,
                        output = "dimObj",
                        spat_unit = spat_unit,
                        feat_type = assay_use,
                        reduction_method = dr_type,
                        name = dr_name
                    )
                    emb_use <- dr_obj[][Seurat::Cells(sobj), ]
                    if (sum(c("loadings", "eigenvalues") %in%
                        names(slot(dr_obj, "misc"))) == 2) {
                        loadings_use <- slot(dr_obj, "misc")$loadings
                        stdev_use <- slot(dr_obj, "misc")$eigenvalues
                        sobj[[dr_name]] <- Seurat::CreateDimReducObject(
                            embeddings = as.matrix(emb_use),
                            loadings = loadings_use,
                            key = paste0(dr_name, "_"),
                            stdev = stdev_use,
                            assay = assay_use
                        )
                    } else {
                        sobj[[dr_name]] <- Seurat::CreateDimReducObject(
                            embeddings = as.matrix(emb_use),
                            key = paste0(dr_name, "_"),
                            assay = assay_use
                        )
                    }
                }
            }
        }



        # network objects
        # expression network
        avail_nn <- list_nearest_networks(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = assay_use
        )
        if (!is.null(avail_nn)) {
            if (nrow(avail_nn) > 0) {
                for (i in seq(nrow(avail_nn))) {
                    nn_name <- avail_nn[i, name]
                    nn_type <- avail_nn[i, nn_type]
                    nn_use <- get_NearestNetwork(
                        gobject = gobject,
                        spat_unit = spat_unit,
                        feat_type = assay_use,
                        nn_network_to_use = nn_type,
                        network_name = nn_name,
                        output = "data.table"
                    )
                    idx1 <- match(nn_use$from, Seurat::Cells(sobj))
                    idx2 <- match(nn_use$to, Seurat::Cells(sobj))
                    edge_weight <- nn_use$weight
                    edge_dist <- nn_use$distance
                    nn_mtx <- Matrix::sparseMatrix(
                        i = idx1, j = idx2, x = edge_dist,
                        dims = c(ncol(sobj), ncol(sobj))
                    )
                    rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
                    nn_name <- paste0(nn_name)
                    sGraph <- Seurat::as.Graph(nn_mtx)
                    sGraph@assay.used <- assay_use
                    sobj[[nn_name]] <- sGraph
                }
            }
        }
    }

    # spatial coordinates
    loc_use <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        output = "spatLocsObj",
        copy_obj = TRUE,
        ... # allow setting of spat_loc_name through additional params
    )

    # flip y vals

    loc_use <- flip(loc_use)[] %>%
        data.table::setDF()

    rownames(loc_use) <- loc_use$cell_ID
    sobj <- Seurat::AddMetaData(sobj, metadata = loc_use)
    # add spatial coordinates as new dim reduct object
    loc_2 <- loc_use[, c("sdimx", "sdimy")]
    colnames(loc_2) <- c("spatial_1", "spatial_2")
    sobj[["spatial"]] <- Seurat::CreateDimReducObject(
        embeddings = as.matrix(loc_2),
        assay = names(sobj@assays)[1],
        key = "spatial_"
    )



    # spatial network
    avail_sn <- list_spatial_networks(gobject = gobject, spat_unit = spat_unit)
    if (!is.null(avail_sn)) {
        if (nrow(avail_sn) > 0) {
            sn_all <- avail_sn[, name]
            for (i in sn_all) {
                snt_use <- get_spatialNetwork(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    name = i,
                    output = "networkDT"
                )
                idx1 <- match(snt_use$from, Seurat::Cells(sobj))
                idx2 <- match(snt_use$to, Seurat::Cells(sobj))
                edge_weight <- snt_use$weight
                edge_dist <- snt_use$distance
                nn_mtx <- Matrix::sparseMatrix(
                    i = idx1, j = idx2, x = edge_dist,
                    dims = c(ncol(sobj), ncol(sobj))
                )
                rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
                nn_name <- paste0(i)
                spatGraph <- Seurat::as.Graph(nn_mtx)
                spatGraph@assay.used <- names(sobj@assays)[1]
                sobj[[nn_name]] <- spatGraph
            }
        }
    }

    all_x <- NULL
    all_y <- NULL

    gimgs <- getGiottoImage(gobject, name = ":all:")
    if (length(gimgs) > 0) {
        for (i in seq_along(gimgs)) {
            gimg <- gimgs[[i]]
            key <- objName(gimg)
            imagerow <- loc_use$sdimy
            imagecol <- loc_use$sdimx
            img_array <- as(gimg, "array")
            img_array <- img_array / 255
            coord <- data.frame(
                imagerow = imagerow, imagecol = imagecol,
                row.names = loc_use$cell_ID
            )
            scalef <- .estimate_scalefactors(
                gimg,
                res_type = res_type,
                spatlocs = loc_use
            )
            # There does not seem to be a way to tell seurat which image type
            # you are using. The lowres scalefactor seems to be the important
            # one in mapping the image
            scalefactors <- Seurat::scalefactors(
                spot = scalef$spot,
                fiducial = scalef$fiducial,
                hires = scalef$hires,
                lowres = scalef[res_type] # this looks like the main one
                # so instead of strictly supplying lowres scalef, we use
                # the scalef belonging to whichever image was used in Giotto
                # since we allow use non-lowres images
            )
            # see https://github.com/satijalab/seurat/issues/3595
            if (tech != 0) {
                if (tech == "Xenium") {
                    coord1 <- coord
                    coord$cell_id <- rownames(coord)
                    coord <- coord[, c("cell_id", "imagerow", "imagecol")]
                    segmentations.data <- list(
                        "centroids" = SeuratObject::CreateCentroids(coord1),
                        "segmentation" = SeuratObject::CreateSegmentation(coord)
                    )
                    coords <- SeuratObject::CreateFOV(
                        coords = segmentations.data,
                        type = c("segmentation", "centroids"),
                        assay = "rna"
                    )
                    fov <- "default_fov"
                    sobj[[fov]] <- coords
                } else {
                    newV1 <- new(
                        Class = "VisiumV1",
                        image = img_array,
                        scale.factors = scalefactors,
                        coordinates = coord,
                        spot.radius =
                            scalef$fiducial * scalef$lowres / max(dim(img_array)),
                        key = paste0(key, "_")
                    )

                    sobj@images[[key]] <- newV1
                }
            }
        }
    }



    return(sobj)
}

#' @param x image object
#' @param res_type type of 10x image output resolution
#' @param spatlocs a data.frame of spatial locations coordinates
#' @noRd

.estimate_scalefactors <- function(x,
    res_type = c("hires", "lowres", "fullres"),
    spatlocs) {
    res_type <- match.arg(res_type,
        choices = c("hires", "lowres", "fullres")
    )
    pxdims <- dim(x)[seq(2)]
    edims <- range(ext(x))
    scalef <- mean(pxdims / edims)
    # assume that lowres and hires follow a general ratio
    # may not be that important since the scalefactor should theoretically
    # only matter for the image res that we are using
    # this ratio is roughly 3.333334 based on Visium BreastCancerA1 dataset
    res_ratio <- 3.333334
    # fullres should have a scalef of roughly 1.
    # No way to guess hires or lowres scalefs so use arbitrary values.
    hres_scalef <- switch(res_type,
        "hires" = scalef,
        "lowres" = scalef * res_ratio,
        "fullres" = 0.08250825 # arbitrary
    )
    lres_scalef <- switch(res_type,
        "hires" = scalef / res_ratio,
        "lowres" = scalef,
        "fullres" = 0.02475247 # arbitrary
    )
    # spot diameter and fid diameter are variable based on how spatial info was
    # mapped to the image. Estimate this by getting the center to center
    # px distance vs fullsize px dims ratio.
    # ! fullsize px dims is the same as edims !
    coords <- data.table::as.data.table(spatlocs)
    # create a delaunay
    dnet <- createNetwork(
        as.matrix(coords[, c("sdimx", "sdimy")]),
        type = "delaunay",
        method = "geometry",
        include_distance = TRUE,
        as.igraph = FALSE,
        include_weight = TRUE,
        verbose = FALSE
    )

    # expect center to center be most common edge distance
    # this gives CC dist as fullres px distance
    distances <- sort(unique(dnet$distance))
    cc_px <- distances[which.max(table(dnet$distance))]
    # assume constant ratios between diameters and cc_px
    fid_cc_ratio <- 1.045909
    fid_diam <- cc_px * fid_cc_ratio
    spot_cc_ratio <- 0.6474675
    spot_diam <- cc_px * spot_cc_ratio
    scalef_list <- list(
        spot = spot_diam,
        fiducial = fid_diam,
        hires = hres_scalef,
        lowres = lres_scalef
    )
    return(scalef_list)
}


#' @title Convert a Seurat object to a Giotto object
#' @name seuratToGiotto
#' @description Convert a Seurat object with the mapping that fits its
#' structure. The object is profiled, the matching mapping is selected from
#' those bundled with GiottoConnect (or the one given), and the conversion runs
#' through [convertWithMapping()]. Object versions are handled by choosing a
#' mapping, not by separate functions.
#' @param sobject Seurat object
#' @param mapping optional `giottoMapping`, mapping file or bundled mapping id.
#' By default the bundled mapping that applies to `sobject` is used.
#' @param source `NULL` for an in-memory giotto object, or a directory for a
#' GiottoDisk-backed object
#' @param verbose logical
#' @returns giotto object; see [mappingProvenance()] for what was converted
#' and which losses apply
#' @seealso [findMapping()], [convertWithMapping()], [validateMapping()]
#' @export
seuratToGiotto <- function(sobject, mapping = NULL, source = NULL, verbose = FALSE) {
    profile <- profileObject(sobject)
    if (is.null(mapping)) {
        found <- findMapping(profile, ecosystem = "seurat")
        if (found$status == "none") {
            stop(wrap_txt("No bundled mapping applies to this Seurat object:",
                paste(found$reasons, collapse = "; ")), call. = FALSE)
        }
        if (found$status == "repair_needed") {
            warning(wrap_txt(sprintf(
                "Structure differs from mapping '%s' (fingerprint %s, expected %s); converting with re-bound params",
                found$mapping$mapping$id, found$fingerprint,
                found$mapping$mapping$fingerprint %null% "unset"
            )), call. = FALSE)
        }
        mapping <- found$mapping
    }
    convertWithMapping(sobject, mapping, profile = profile, source = source,
        verbose = verbose)
}
