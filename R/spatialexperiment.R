#' @import GiottoClass
#' @import GiottoUtils
NULL


#' Utility function to convert a Giotto object to a SpatialExperiment object.
#'
#' @param gobject Input Giotto object to convert to a SpatialExperiment object
#' @param verbose A boolean value specifying if progress messages should be
#' displayed or not. Default \code{TRUE}.
#' @param giottoObj deprecated
#'
#' @returns A SpatialExperiment object that contains data from the input Giotto
#' object.
#' @examples
#' \dontrun{
#' mini_gobject <- GiottoData::loadGiottoMini("vizgen")
#' giottoToSpatialExperiment(mini_gobject)
#' }
#' @export
giottoToSpatialExperiment <- function(gobject,
    verbose = TRUE,
    giottoObj = deprecated()) {
    # deprecation
    gobject <- deprecate_param(giottoObj, gobject,
        fun = "giottoToSpatialExperiment",
        when = "0.4.6"
    )

    spat_unit <- NULL

    # Load required packages
    # package_check(pkg_name = "SummarizedExperiment", repository = 'Bioc')
    # SP should load this?
    # package_check(pkg_name = "SingleCellExperiment", repository = 'Bioc')
    # SP should load this?
    package_check(pkg_name = "SpatialExperiment", repository = "Bioc")
    package_check(pkg_name = "S4Vectors", repository = "Bioc")

    # SpatialExperiment objects list (one object for each spatial unit)
    speList <- list()

    # Expression Matrices
    giottoExpr <- list_expression(gobject)

    # Iterate over spatial units
    spatialUnits <- unique(giottoExpr$spat_unit) # a function to get spat units?
    for (su in seq(spatialUnits)) {
        if (verbose) {
            message("Processing spatial unit: '", spatialUnits[su], "'")
        }
        # Check if expression matrices exist in input object
        if (!is.null(giottoExpr)) {
            if (verbose) {
                message(
                    "Copying expression matrix: '", giottoExpr[1]$name,
                    "' for spatial unit: '", spatialUnits[su], "'"
                )
            }
            exprMat <- getExpression(
                gobject = gobject,
                spat_unit = spatialUnits[su],
                feat_type = giottoExpr[1]$feat_type,
                values = giottoExpr[1]$name,
                output = "matrix"
            )
            names(rownames(exprMat)) <- NULL
            names(colnames(exprMat)) <- NULL
            exprMat <- list(exprMat)
            names(exprMat)[1] <- giottoExpr[1]$name
            # Creating SPE object with first expression matrix
            spe <- SpatialExperiment::SpatialExperiment(assays = exprMat)
            SummarizedExperiment::assayNames(spe) <- paste0(
                giottoExpr[1]$name, "_",
                giottoExpr[1]$feat_type, "_",
                spatialUnits[su]
            )
            giottoExpr <- giottoExpr[-1, ]
        } else {
            stop("The input Giotto object must contain atleast one expression
                matrix.")
        }

        # Copying remaining expression matrices if they exist
        if (nrow(giottoExpr[spat_unit == spatialUnits[su]]) > 0) {
            for (i in seq(nrow(giottoExpr))) {
                if (verbose) {
                    message(
                        "Copying expression matrix: '",
                        giottoExpr[i]$name, "' for spatial unit: '",
                        spatialUnits[su], "'"
                    )
                }
                # SPE does not have specific slots for different units,
                # instead joining multiple unit names to identify them
                SummarizedExperiment::assay(
                    spe,
                    paste0(
                        giottoExpr[i]$name, "_",
                        giottoExpr[i]$feat_type, "_",
                        spatialUnits[su]
                    ),
                    withDimnames = FALSE
                ) <- getExpression(
                    gobject = gobject,
                    spat_unit = spatialUnits[su],
                    feat_type = giottoExpr[i]$feat_type,
                    values = giottoExpr[i]$name,
                    output = "matrix"
                )
            }
        }

        # Cell Metadata to ColData
        pData <- pDataDT(gobject = gobject, spat_unit = spatialUnits[su])
        if (nrow(pData) > 0) {
            if (verbose) {
                message(
                    "Copying phenotype data for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            SummarizedExperiment::colData(spe) <- S4Vectors::DataFrame(
                pData,
                row.names = pData$cell_ID
            )
        } else {
            if (verbose) {
                message("No phenotype data found in input Giotto object")
            }
        }

        # Feature Metadata to RowData
        fData <- fDataDT(gobject = gobject, spat_unit = spatialUnits[su])
        if (nrow(fData) > 0) {
            if (verbose) {
                message(
                    "Copying feature metadata for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            SummarizedExperiment::rowData(spe) <- fData
        } else {
            if (verbose) {
                message("No feature metadata found in input Giotto object")
            }
        }

        # Spatial Locations to Spatial Coordinates
        spatialLocs <- get_spatial_locations(
            gobject = gobject,
            spat_unit = spatialUnits[su],
            output = "data.table"
        )
        if (!is.null(spatialLocs)) {
            if (verbose) {
                message(
                    "Copying spatial locations for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            if (all(colnames(spatialLocs[, seq_along(2)]) == c("sdimx", "sdimy"))) {
                spatialLocs <- spatialLocs[, c("sdimx", "sdimy"), drop = FALSE]
            } else {
                # Rename the first two columns to sdimx and sdimy
                colnames(spatialLocs)[seq_along(2)] <- c("sdimx", "sdimy")
                spatialLocs <- spatialLocs[, c("sdimx", "sdimy"), drop = FALSE]
            }
            SpatialExperiment::spatialCoords(spe) <- data.matrix(spatialLocs)
        } else {
            if (verbose) {
                message("No spatial locations found in the input Giotto object")
            }
        }

        # DimReductions
        giottoReductions <- list_dim_reductions(
            gobject = gobject,
            spat_unit = spatialUnits[su]
        )
        if (!is.null(giottoReductions)) {
            if (verbose) {
                message(
                    "Copying reduced dimensions for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            for (i in seq(nrow(giottoReductions))) {
                SingleCellExperiment::reducedDim(
                    spe,
                    giottoReductions[i]$name
                ) <- get_dimReduction(
                    gobject = gobject,
                    reduction = "cells",
                    spat_unit = spatialUnits[su],
                    feat_type = giottoReductions[i]$feat_type,
                    reduction_method = giottoReductions[i]$dim_type,
                    name = giottoReductions[i]$name,
                    output = "data.table"
                )
            }
        } else {
            if (verbose) {
                message("No reduced dimensions found in the input Giotto
                        object")
            }
        }


        # NN Graph
        giottoNearestNetworks <- list_nearest_networks(
            gobject = gobject,
            spat_unit = spatialUnits[su]
        )
        if (!is.null(giottoNearestNetworks)) {
            if (verbose) {
                message(
                    "Copying nearest networks for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            for (i in seq(nrow(giottoNearestNetworks))) {
                nn_network <- get_NearestNetwork(
                    gobject = gobject,
                    spat_unit = spatialUnits[su],
                    nn_network_to_use = giottoNearestNetworks[i]$type,
                    network_name = giottoNearestNetworks[i]$name,
                    output = "data.table"
                )

                # SPE stores in colpairs, with col indices instead of colnames
                cell1 <- match(nn_network$from, pData$cell_ID)
                cell2 <- match(nn_network$to, pData$cell_ID)

                SingleCellExperiment::colPair(
                    spe,
                    giottoNearestNetworks[i]$name
                ) <- S4Vectors::SelfHits(
                    cell1, cell2,
                    nnode = ncol(spe),
                    nn_network[, seq(-1, -2)] # removing from and to
                )
            }
        } else {
            if (verbose) {
                message("No nearest networks found in the input Giotto object")
            }
        }

        # Spatial Networks
        giottoSpatialNetworks <- list_spatial_networks(
            gobject = gobject,
            spat_unit = spatialUnits[su]
        )
        if (!is.null(giottoSpatialNetworks)) {
            if (verbose) {
                message(
                    "Copying spatial networks for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            for (i in seq(nrow(giottoSpatialNetworks))) {
                sp_network <- get_spatialNetwork(
                    gobject = gobject,
                    spat_unit = spatialUnits[su],
                    name = giottoSpatialNetworks[i]$name,
                    output = "networkDT"
                )

                # spe stores in colpairs, with col indices instead of colnames
                cell1 <- match(sp_network$from, pData$cell_ID)
                cell2 <- match(sp_network$to, pData$cell_ID)

                SingleCellExperiment::colPair(
                    spe,
                    giottoSpatialNetworks[i]$name
                ) <- S4Vectors::SelfHits(
                    cell1, cell2,
                    nnode = ncol(spe),
                    sp_network[, seq(-1, -2)] # removing from and to
                )
            }
        } else {
            if (verbose) {
                message("No spatial networks found in the input Giotto object")
            }
        }

        # SpatialImages
        giottoImages <- list_images(gobject = gobject)
        if (!is.null(giottoImages)) {
            for (i in seq(nrow(giottoImages))) {
                img <- get_giottoImage(
                    gobject = gobject,
                    image_type = giottoImages[i]$img_type,
                    name = giottoImages[i]$name
                )

                if (!is.null(img@file_path)) {
                    if (verbose) message("Copying spatial image: ", img@name)

                    tryCatch(
                        expr = {
                            spe <- SpatialExperiment::addImg(spe,
                                sample_id = spe$sample_id[i],
                                image_id = img@name,
                                imageSource = img@file_path,
                                scaleFactor = mean(img@scale_factor),
                                load = TRUE
                            )
                        },
                        error = function(e) {
                            message(
                                "Error copying spatial image: ",
                                img@name,
                                ". Please check if the image path is
                                    correct and the image exists at that path."
                            )
                        }
                    )
                } else {
                    if (verbose) {
                        message(
                            "\t - Skipping image with NULL file path: ",
                            img@name
                        )
                    }
                }
                S4Vectors::metadata(spe)[[img@name]] <- img
            }
        } else {
            if (verbose) {
                message("No spatial images found in the input Giotto object")
            }
        }

        if (verbose) message("")

        # Add spe for current spatial unit to speList
        speList[[su]] <- spe
    }

    # return list of spe objects
    return(speList)
}


#' Utility function to convert a SpatialExperiment object to a Giotto object
#'
#' @param spe Input SpatialExperiment object to convert to a Giotto object.
#' @param python_path Specify the path to python.
#' @param nn_network Specify the name of the nearest neighbour network(s)
#' in the input SpatialExperiment object. Default \code{NULL} will use
#' all existing networks.
#' @param sp_network Specify the name of the spatial network(s) in the input
#' SpatialExperiment object. Default \code{NULL} will use all existing
#' networks. This can be a vector of multiple network names.
#' @param verbose A boolean value specifying if progress messages should
#' be displayed or not. Default \code{TRUE}.
#' @returns Giotto object
#' @examples
#' \dontrun{
#' library(SpatialExperiment)
#' example(read10xVisium, echo = FALSE)
#' spatialExperimentToGiotto(spe)
#' }
#' @export
spatialExperimentToGiotto <- function(
        spe,
        python_path,
        nn_network = NULL,
        sp_network = NULL,
        verbose = TRUE) {
    # Create giotto instructions and set python path
    instrs <- createGiottoInstructions(python_path = python_path)

    # Create Giotto object with first matrix
    exprMats <- SummarizedExperiment::assays(spe)
    exprMatsNames <- SummarizedExperiment::assayNames(spe)
    firstMatrix <- exprMats[[1]]

    # check unique colnames
    if (length(unique(colnames(firstMatrix))) !=
        length(colnames(firstMatrix))) {
        colnames(firstMatrix) <- make.names(colnames(firstMatrix),
            unique = TRUE
        )
    }

    if (verbose) {
        message("Creating Giotto object with ", exprMatsNames[1], " matrix")
    }
    suppressWarnings(suppressMessages(giottoObj <- createGiottoObject(
        expression = firstMatrix,
        instructions = instrs
    )))
    exprMats[[1]] <- NULL
    exprMatsNames <- exprMatsNames[-1]

    # Copying remaining matrices
    if (length(exprMats) > 0) {
        for (i in seq(exprMats)) {
            if (verbose) {
                message("Copying expression matrix: ", exprMatsNames[i])
            }
            exprObj <- create_expr_obj(
                name = exprMatsNames[i],
                exprMat = exprMats[[i]]
            )
            giottoObj <- set_expression_values(
                gobject = giottoObj,
                values = exprObj
            )
        }
    }

    # Phenotype Data
    pData <- SummarizedExperiment::colData(spe)

    # To handle cell_ID duplication issues
    if ("cell_ID" %in% colnames(pData)) {
        pData$`cell_ID` <- NULL
    }
    if (nrow(pData) > 0) {
        if (verbose) message("Copying phenotype data")
        giottoObj <- addCellMetadata(
            gobject = giottoObj,
            new_metadata = as.data.table(pData)
        )
    }

    # Feature Metadata
    fData <- SummarizedExperiment::rowData(spe)
    if (nrow(fData) > 0) {
        if (verbose) message("Copying feature metadata")
        giottoObj <- addFeatMetadata(
            gobject = giottoObj,
            new_metadata = as.data.table(fData)
        )
    }

    # Reduced Dimensions
    redDims <- SingleCellExperiment::reducedDims(spe)
    redDimsNames <- SingleCellExperiment::reducedDimNames(spe)
    if (length(redDims) > 0) {
        for (i in seq(length(redDims))) {
            if (verbose) message("Copying reduced dimensions")
            dimRedObj <- create_dim_obj(
                name = redDimsNames[i],
                coordinates = redDims[[i]],
                reduction_method = redDimsNames[i]
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            giottoObj <- set_dimReduction(
                gobject = giottoObj,
                dimObject = dimRedObj
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }
    }

    # Spatial Locations
    spatialLocs <- SpatialExperiment::spatialCoords(spe)
    if (ncol(spatialLocs) > 0) {
        if (verbose) message("Copying spatial locations")
        spatialLocsDT <- data.table::data.table(
            sdimx = spatialLocs[, 1],
            sdimy = spatialLocs[, 2],
            cell_ID = colnames(spe)
        )
        spatLocsObj <- create_spat_locs_obj(
            name = "spatLocs",
            coordinates = spatialLocsDT
        )
        giottoObj <- set_spatial_locations(
            gobject = giottoObj,
            spatlocs = spatLocsObj
        )
    }

    # Spatial Images
    spatialImages <- SpatialExperiment::imgData(spe)
    if (nrow(spatialImages) > 0) {
        for (i in seq(nrow(spatialImages))) {
            if (verbose) message("Copying spatial images")
            spImg <- SpatialExperiment::getImg(
                spe,
                spatialImages[i, "sample_id"],
                spatialImages[i, "image_id"]
            )
            mObject <- magick::image_read(grDevices::as.raster(spImg))
            giottoImage <- createGiottoImage(
                gobject = giottoObj,
                mg_object = mObject,
                scale_factor = spatialImages[i, "scaleFactor"]
            )
            giottoObj <- addGiottoImage(
                gobject = giottoObj,
                images = list(giottoImage)
            )
        }
    }

    # Networks
    networks <- SingleCellExperiment::colPairs(spe)
    # Spatial Networks
    if (!is.null(sp_network)) {
        if (all(sp_network %in% names(networks))) {
            for (i in seq(sp_network)) {
                if (verbose) message("Copying spatial networks")
                networkDT <- as.data.table(networks[[sp_network[i]]])
                networkDT$to <- colnames(spe)[networkDT$to]
                networkDT$from <- colnames(spe)[networkDT$from]
                spatNetObj <- create_spat_net_obj(
                    networkDT = networkDT
                )
                giottoObj <- set_spatialNetwork(
                    gobject = giottoObj,
                    spatial_network = spatNetObj,
                    name = sp_network[i]
                )
                networks[[sp_network[i]]] <- NULL
            }
        }
    }

    # Nearest Neighbour Networks
    if (!is.null(nn_network)) {
        if (nn_network %in% names(networks)) {
            for (i in seq(nn_network)) {
                if (verbose) message("Copying nearest neighbour networks")
                nnNetObj <- create_nn_net_obj(
                    name = nn_network[i],
                    igraph = networks[[nn_network[i]]]
                )
                giottoObj <- set_NearestNetwork(
                    gobject = giottoObj,
                    nn_network = nnNetObj
                )
                networks[[nn_network[i]]] <- NULL
            }
        }
    }

    # if not specified, storing remaining as NN
    if (length(networks) > 0) {
        for (i in seq(networks)) {
            if (verbose) message("Copying additional networks")
            nnNetObj <- create_nn_net_obj(
                name = names(networks)[i],
                igraph = networks[[i]]
            )
            giottoObj <- set_NearestNetwork(
                gobject = giottoObj,
                nn_network = nnNetObj
            )
        }
    }

    return(giottoObj)
}





# SpatialExperiment methods ####

# register methods if SpatialExperiment is present
if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
    setMethod(
        SingleCellExperiment::colData,
        signature("giotto"),
        function(x, spat_unit = NULL, feat_type = NULL, ...) {
            pDataDT(
                gobject = x,
                spat_unit = spat_unit,
                feat_type = feat_type, ...
            ) %>%
                S4Vectors::DataFrame()
        }
    )

    setMethod(
        # TODO generate DataFrame of ENSG and gene symbols
        SingleCellExperiment::rowData,
        signature("giotto"),
        function(x, spat_unit = NULL, feat_type = NULL, ...) {
            fDataDT(
                gobject = x,
                spat_unit = spat_unit,
                feat_type = feat_type, ...
            ) %>%
                S4Vectors::DataFrame()
        }
    )
}
