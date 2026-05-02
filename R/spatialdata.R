#' @import GiottoClass
#' @import GiottoUtils
NULL




## SpatialData object to Giotto ####

#' @title Convert SpatialData to Giotto
#' @name spatialdataToGiotto
#' @description Converts a saved SpatialData object into a Giotto object
#'
#' @param spatialdata_path path to SpatialData object
#' @param n_key_added equivalent of "key_added" argument from scanpy.pp.neighbors().
#' If multiple spatial networks are in the anndata object, a list of key_added
#' terms may be provided.
#' If converting an anndata object from giottoToAnnData, the keys are saved in `.uns['NN_keys']`
#' and all keys are used in conversion unless specified in the function call.
#' Cannot be "spatial". This becomes the name of the nearest network in the gobject.
#' @param spatial_n_key_added
#' equivalent of "key_added" argument from squidpy.gr.spatial_neighbors.
#' If multiple spatial networks are in the anndata object, a list of key_added
#' terms may be provided.
#' If converting a SpatialData object from giottoToSpatialData, the keys are saved in `.uns['SN_keys']`
#' and all keys are used in conversion unless specified in the function call.
#' Cannot be the same as n_key_added.
#' @param delaunay_spat_net binary parameter for spatial network. If TRUE,
#' the spatial network is a delaunay network.
#' @param spat_enrich_key_added
#' list of names of spatial enrichment annotations added to the SpatialData object.
#' If converting an anndata object from giottoToAnnData and the original Giotto object had
#' spatial enrichment annotations, the keys are saved in `.uns['SE_keys']`
#' and all keys are used in conversion unless specified in the function call.
#' @param spat_unit desired spatial unit for conversion, default NULL
#' @param feat_type desired feature type for conversion, default NULL
#' @param python_path path to python executable within a conda/miniconda environment
#' @param env_name name of environment containing python_path executable
#'
#' @return Giotto object
#' @details Function in beta. Converts a structured SpatialData file into a Giotto object.
#'    The returned Giotto Object will take default insructions with the
#'    exception of the python path, which may be customized.
#'    See \code{\link{changeGiottoInstructions}} to modify instructions after creation.
#' @export

spatialdataToGiotto <- function(spatialdata_path = NULL,
    n_key_added = NULL,
    spatial_n_key_added = NULL,
    delaunay_spat_net = TRUE,
    spat_enrich_key_added = NULL,
    spat_unit = NULL,
    feat_type = NULL,
    python_path = NULL,
    env_name = NULL) {
    # File check
    if (is.null(spatialdata_path)) {
        stop("Please provide a path to SpatialData object for conversion.\n")
    }
    if (!file.exists(spatialdata_path)) {
        stop("The provided path to SpatialData object does not exist.\n")
    }

    # Initialize reticulate
    instrs <- createGiottoInstructions(
        python_path = python_path,
        show_plot = FALSE,
        return_plot = FALSE,
        save_plot = TRUE,
        save_dir = NULL,
        plot_format = NULL,
        dpi = NULL,
        units = NULL,
        height = NULL,
        width = NULL,
        is_docker = FALSE,
        plot_count = 0,
        fiji_path = NULL,
        no_python_warn = FALSE
    )

    # Check spatialdata dependencies
    package_check(
        pkg_name = "spatialdata",
        repository = "pip:spatialdata"
    )

    # Import sd2g, a python module for parsing SpatialData
    sd2g_path <- system.file("python", "sd2g.py", package = "GiottoConnect")
    reticulate::source_python(sd2g_path)
    sdata <- read_spatialdata_from_path(spatialdata_path)

    if (is.null(n_key_added)) {
        n_key_added <- list.files(
            file.path(spatialdata_path, "giotto_meta"),
            pattern = "nn_network", full.names = TRUE
        )
        n_key_added <- n_key_added %none% NULL
    }
    if (is.null(spatial_n_key_added)) {
        spatial_n_key_added <- list.files(
            file.path(spatialdata_path, "giotto_meta"),
            pattern = "spatial_network", full.names = TRUE
        )
        spatial_n_key_added <- spatial_n_key_added %none% NULL
    }

    # Create baseline Giotto object
    gobject <- createGiottoObject(
        instructions = instrs
    )

    # Extract expression matrices
    expr_df_dict <- extract_expression(sdata)

    su_list_initial <- list()
    ft_list_initial <- list()
    su_list <- list()
    ft_list <- list()

    # Set expression data
    if (!is.null(spat_unit)) {
        su_list_initial <- spat_unit
    }
    if (!is.null(feat_type)) {
        ft_list_initial <- feat_type
    }

    for (key in names(expr_df_dict)) {
        # If table name is used to extract su and ft
        if (length(su_list_initial) == 0 && length(ft_list_initial) == 0 && length(strsplit(key, "_")[[1]]) >= 2) {
            parts <- strsplit(key, "_")[[1]]
            su <- parts[1]
            su_list <- c(su_list, su)
            if (length(parts) == 3) {
                ft <- paste0(parts[2], "_", parts[3])
            } else {
                ft <- parts[2]
            }
            ft_list <- c(ft_list, ft)
            cat("Spatial unit and feature type extracted from table name.\n")
        } else if (length(su_list_initial) > 0 && length(ft_list_initial) == 0) {
            su <- spat_unit
            ft <- "rna"
            ft_list <- c(ft_list, ft)
            cat("Only the spatial unit has been specified in function call. Default feature type has been set to rna.\n")
        } else if (length(su_list_initial) == 0 && length(ft_list_initial) > 0) {
            su <- "cell"
            ft <- feat_type
            su_list <- c(su_list, su)
            cat("Only the feature type has been specified in function call. Default spatial unit has been set to cell.\n")
        } else if (length(su_list_initial) > 0 && length(ft_list_initial) > 0) {
            su <- spat_unit[[1]]
            ft <- feat_type[[1]]
            cat("Using spat_unit and feat_type specified in function call.\n")
        } else {
            su <- "cell"
            ft <- "rna"
            su_list <- c(su_list, su)
            ft_list <- c(ft_list, ft)
            cat("Default spatial unit and feature type have been set to [cell][rna]. If you wish to set a specific\nspatial unit and feature type, please specify it in the function call.\n")
        }

        gobject <- setExpression(
            gobject,
            x = createExprObj(expr_df_dict[[key]], name = "raw"),
            spat_unit = su,
            feat_type = ft
        )
    }

    # ID spat_unit and feat_type if not already provided.
    if (is.null(spat_unit) && is.null(feat_type)) {
        spat_unit <- su_list
        feat_type <- ft_list
    } else if (is.null(spat_unit) && !is.null(feat_type)) {
        spat_unit <- su_list
    } else if (!is.null(spat_unit) && is.null(feat_type)) {
        feat_type <- ft_list
    }

    for (su in spat_unit) {
        wrap_msg("Spatial unit(s)", su, "will be used in conversion.")
    }
    for (ft in feat_type) {
        wrap_msg("Feature type(s)", ft, "will be used in conversion.")
    }

    cID_dict <- extract_cell_IDs(sdata)
    fID_dict <- extract_feat_IDs(sdata)

    ## Layers
    lay_dict <- extract_layer_names(sdata)
    if (!is.null(lay_dict)) {
        for (key in names(lay_dict)) {
            parts <- strsplit(key, "_")[[1]]
            layer_su <- parts[1]
            if (length(parts) == 3) {
                layer_ft <- paste0(parts[2], "_", parts[3])
            } else {
                layer_ft <- parts[2]
            }
            for (l_n in lay_dict[[key]]) {
                lay <- extract_layered_data(sdata, key, layer_name = l_n)
                if ("data.frame" %in% class(lay)) {
                    names(lay) <- fID_dict[[key]]
                    row.names(lay) <- cID_dict[[key]]
                } else {
                    lay@Dimnames[[1]] <- fID_dict[[key]]
                    lay@Dimnames[[2]] <- cID_dict[[key]]
                }
                parts <- strsplit(l_n, "_")[[1]]
                if (length(parts) == 4) {
                    name <- parts[4]
                } else {
                    name <- parts[3]
                }
                layExprObj <- createExprObj(lay, name = name)
                gobject <- set_expression_values(
                    gobject = gobject,
                    spat_unit = layer_su,
                    feat_type = layer_ft,
                    values = layExprObj
                )
            }
        }
    }

    # Extract polygons
    polygon_dict <- extract_polygons(sdata)
    for (su in names(polygon_dict)) {
        gpoly_dt <- as.data.table(polygon_dict[[su]])
        gpoly <- createGiottoPolygon(
            gpoly_dt
        )
        gobject <- setPolygonInfo(gobject, gpoly, name = su)
    }

    gobject <- update_giotto_params(
        gobject = gobject,
        description = "_AnnData_Conversion"
    )
    # Extract spatial locations
    spatial_dict <- extract_spatial(sdata)

    # Set spatial locations
    for (key in names(spatial_dict)) {
        if (length(strsplit(key, "_")[[1]]) >= 2) {
            parts <- strsplit(key, "_")[[1]]
            spat_loc_su <- parts[1]
            gobject <- setSpatialLocations(
                gobject,
                x = createSpatLocsObj(spatial_dict[[key]], name = "raw"),
                spat_unit = spat_loc_su
            )
        } else {
            gobject <- setSpatialLocations(
                gobject,
                x = createSpatLocsObj(spatial_dict[[key]], name = "raw"),
                spat_unit = key
            )
        }
    }
    sp_dict <- parse_obsm_for_spat_locs(sdata)

    # Spatial Enrichment
    spat_enrich_all <- extract_spat_enrich(
        sdata,
        key_added = spat_enrich_key_added
    )
    num_spat_enrich <- length(spat_enrich_all)

    if (is.null(spat_enrich_key_added) && !is.null(spat_enrich_all)) {
        num_spat_enrich <- length(spat_enrich_all)
    }

    for (i in seq_len(num_spat_enrich)) {
        if (inherits(spat_enrich_key_added, "list")) {
            spat_enrich_key_added_it <- spat_enrich_key_added[[i]]
        } else {
            sk <- names(spat_enrich_all)[[i]]
            split_key <- strsplit(gsub("[()']", "", sk), ", ")[[1]]
            tn <- split_key[1]
            tn_parts <- strsplit(tn, "_")[[1]]
            se_su <- tn_parts[1]
            se_ft <- tn_parts[2]
            se_name <- split_key[2]
            spat_enrich_key_added_it <- se_name
        }
        se <- spat_enrich_all[[i]]
        spatEnrObj <- createSpatEnrObj(
            se,
            name = se_name,
            spat_unit = se_su,
            feat_type = se_ft
        )
        gobject <- setSpatialEnrichment(
            gobject = gobject,
            spatEnrObj
        )
    }

    ## Spatial Network
    s_weights_sd <- NULL
    num_SN_nets <- length(spatial_n_key_added)
    # Check for the case where NULL is provided, since the
    # anndata object takes the default value for SN

    if (is.null(spatial_n_key_added) && !is.null(extract_SN_connectivities(sdata, key_added = spatial_n_key_added))) {
        num_SN_nets <- length(extract_SN_connectivities(sdata, key_added = spatial_n_key_added))
    }

    for (i in seq_len(num_SN_nets)) {
        if (inherits(spatial_n_key_added, "list")) {
            spatial_n_key_added_it <- spatial_n_key_added[[i]]
        } else {
            spatial_n_key_added_it <- spatial_n_key_added
        }

        s_weights_list <- extract_SN_connectivities(
            sdata,
            key_added = spatial_n_key_added_it
        )
        if (!is.null(s_weights_list)) {
            for (connectivity_key in names(s_weights_list)) {
                s_weights_sd <- s_weights_list[[connectivity_key]]
                split_key <- strsplit(connectivity_key, ", ")[[1]]
                tn <- gsub("[()']", "", split_key[1])
                sk <- gsub("[()']", "", split_key[2])

                s_distances_sd <- extract_SN_distances(
                    sdata,
                    key_added = spatial_n_key_added_it,
                    tn = tn, sn_key_list = sk
                )
                ij_matrix <- methods::as(s_distances_sd, "TsparseMatrix")
                from_idx <- ij_matrix@i + 1 # zero index!!!
                to_idx <- ij_matrix@j + 1 # zero index!!!

                # pre-allocate DT variables
                from <- to <- weight <- distance <- from_cell_ID <- to_cell_ID <- uniq_ID <- NULL
                sn_dt <- data.table::data.table(
                    from = from_idx,
                    to = to_idx,
                    weight = s_weights_sd@x,
                    distance = s_distances_sd@x
                )
                cID <- cID_dict[[tn]]
                sn_dt[, from_cell_ID := cID[from]]
                sn_dt[, to_cell_ID := cID[to]]

                sdimx <- "sdimx"
                sdimy <- "sdimy"
                xbegin_name <- paste0(sdimx, "_begin")
                ybegin_name <- paste0(sdimy, "_begin")
                xend_name <- paste0(sdimx, "_end")
                yend_name <- paste0(sdimy, "_end")

                network_DT <- data.table::data.table(
                    from = unlist(sn_dt$from_cell_ID),
                    to = unlist(sn_dt$to_cell_ID),
                    xbegin_name = sp_dict[[tn]][sn_dt$from, sdimx],
                    ybegin_name = sp_dict[[tn]][sn_dt$from, sdimy],
                    xend_name = sp_dict[[tn]][sn_dt$to, sdimx],
                    yend_name = sp_dict[[tn]][sn_dt$to, sdimy],
                    weight = s_weights_sd@x,
                    distance = s_distances_sd@x
                )
                data.table::setnames(network_DT,
                    old = c("xbegin_name", "ybegin_name", "xend_name", "yend_name"),
                    new = c(xbegin_name, ybegin_name, xend_name, yend_name)
                )
                data.table::setorder(network_DT, from, to)
                dist_mean <- get_distance(network_DT, method = "mean")
                dist_median <- get_distance(network_DT, method = "median")
                cellShapeObj <- list(
                    "meanCellDistance" = dist_mean,
                    "medianCellDistance" = dist_median
                )

                # TODO filter network?
                # TODO 3D handling?
                net_name <- strsplit(sk, "_")[[1]][1]
                if (net_name == "Delaunay") {
                    spatObj <- createSpatNetObj(
                        network = network_DT,
                        name = "Delaunay_network",
                        method = "delaunay",
                        networkDT_before_filter = network_DT,
                        cellShapeObj = cellShapeObj
                    )
                } else {
                    net_name <- sub("_connectivities", "", sk)
                    spatObj <- createSpatNetObj(
                        network = network_DT,
                        name = net_name,
                        method = "non-delaunay",
                        networkDT_before_filter = network_DT,
                        cellShapeObj = cellShapeObj
                    )
                }
                gobject <- set_spatialNetwork(
                    gobject = gobject,
                    spatial_network = spatObj,
                    spat_unit = strsplit(tn, "_")[[1]][1]
                )
            }
        }
    }

    # Dimension reduction
    ## Add PCA
    p_dict <- extract_pca(sdata)
    if (!is.null(p_dict)) {
        for (tn in names(p_dict)) {
            for (i in seq_along(p_dict[[tn]][[1]])) {
                p <- p_dict[[tn]][[1]][[i]]
                if (!is.null(p)) {
                    pca_name <- names(p_dict[[tn]][[1]])[[i]]
                    pca <- p[["pca"]]
                    evs <- if (!is.null(p[["eigenvalues"]])) p[["eigenvalues"]] else NULL
                    loads <- if (!is.null(p[["loadings"]])) p[["loadings"]] else NULL
                    parts <- strsplit(tn, "_")[[1]]
                    pca_su <- parts[1]
                    pca_ft <- parts[2]
                    rownames_vec <- if (!is.null(expr_df_dict[[tn]])) colnames(expr_df_dict[[tn]]) else NULL
                    dobj <- createDimObj(
                        coordinates = pca,
                        name = pca_name,
                        spat_unit = pca_su,
                        feat_type = pca_ft,
                        method = "pca",
                        reduction = "cells",
                        provenance = NULL,
                        misc = list(
                            eigenvalues = evs,
                            loadings = loads
                        ),
                        my_rownames = rownames_vec
                    )
                    gobject <- set_dimReduction(
                        gobject = gobject,
                        dimObject = dobj
                    )
                }
            }
        }
    }

    # Add UMAP
    u_dict <- extract_umap(sdata)
    if (!is.null(u_dict)) {
        for (tn in names(u_dict)) {
            for (i in seq_along(u_dict[[tn]][[1]])) {
                u <- u_dict[[tn]][[1]][[i]]
                if (!is.null(u)) {
                    umap_name <- names(u_dict[[tn]][[1]])[[i]]
                    parts <- strsplit(tn, "_")[[1]]
                    umap_su <- parts[1]
                    umap_ft <- parts[2]
                    rownames_vec <- if (!is.null(expr_df_dict[[tn]])) colnames(expr_df_dict[[tn]]) else NULL
                    dobj <- createDimObj(
                        coordinates = u,
                        name = umap_name,
                        spat_unit = umap_su,
                        feat_type = umap_ft,
                        method = "umap",
                        reduction = "cells",
                        provenance = NULL,
                        misc = NULL,
                        my_rownames = rownames_vec
                    )
                    gobject <- set_dimReduction(
                        gobject = gobject,
                        dimObject = dobj
                    )
                }
            }
        }
    }

    # Add tSNE
    t_dict <- extract_tsne(sdata)
    if (!is.null(t_dict)) {
        for (tn in names(t_dict)) {
            for (i in seq_along(t_dict[[tn]][[1]])) {
                t <- t_dict[[tn]][[1]][[i]]
                if (!is.null(t)) {
                    tsne_name <- names(t_dict[[tn]][[1]])[[i]]
                    parts <- strsplit(tn, "_")[[1]]
                    tsne_su <- parts[1]
                    tsne_ft <- parts[2]
                    rownames_vec <- if (!is.null(expr_df_dict[[tn]])) colnames(expr_df_dict[[tn]]) else NULL
                    dobj <- createDimObj(
                        coordinates = t,
                        name = tsne_name,
                        spat_unit = tsne_su,
                        feat_type = tsne_ft,
                        method = "tsne",
                        reduction = "cells",
                        provenance = NULL,
                        misc = NULL,
                        my_rownames = rownames_vec
                    )
                    gobject <- set_dimReduction(
                        gobject = gobject, dimObject = dobj
                    )
                }
            }
        }
    }

    ## Nearest Network
    weights_sd_all <- extract_NN_connectivities(sdata, key_added = n_key_added)
    num_NN_nets <- length(n_key_added)

    if (is.null(n_key_added) && !is.null(weights_sd_all)) {
        num_NN_nets <- length(weights_sd_all)
    }

    for (i in seq_len(num_NN_nets)) {
        if (inherits(n_key_added, "list")) {
            n_key_added_it <- n_key_added[[i]]
        } else {
            n_key_added_it <- n_key_added
        }

        weights_list <- extract_NN_connectivities(sdata, key_added = n_key_added_it)
        if (!is.null(weights_list)) {
            for (connectivity_key in names(weights_list)) {
                weights_sd <- weights_list[[connectivity_key]]
                split_key <- strsplit(connectivity_key, ", ")[[1]]
                tn <- gsub("[()']", "", split_key[1])
                nk <- gsub("[()']", "", split_key[2])

                distances_sd <- extract_NN_distances(
                    sdata,
                    key_added = n_key_added_it,
                    tn = tn, nn_key_list = nk
                )

                nn_dt <- align_network_data(
                    distances = weights_sd, weights = distances_sd
                )

                # pre-allocate DT variables
                from <- to <- weight <- distance <- from_cell_ID <- to_cell_ID <- uniq_ID <- NULL
                nn_dt <- data.table::data.table(nn_dt)
                cID <- cID_dict[[tn]]
                nn_dt[, from_cell_ID := cID[from]]
                nn_dt[, to_cell_ID := cID[to]]
                nn_dt[, from_cell_ID := as.character(from_cell_ID)]
                nn_dt[, to_cell_ID := as.character(to_cell_ID)]
                nn_dt[, uniq_ID := paste0(from, to)]
                nn_dt[order(uniq_ID)]
                nn_dt[, uniq_ID := NULL]
                vert <- unique(x = c(nn_dt$from_cell_ID, nn_dt$to_cell_ID))
                nn_network_igraph <- igraph::graph_from_data_frame(
                    nn_dt[, .(from_cell_ID, to_cell_ID, weight, distance)],
                    directed = TRUE, vertices = vert
                )

                nn_info <- extract_NN_info(
                    sdata = sdata, key_added = n_key_added_it, tn = tn
                )

                net_type <- "kNN" # anndata default
                if (("sNN" %in% n_key_added_it) & !is.null(n_key_added_it)) {
                    net_type <- "sNN"
                    net_name <- paste0(n_key_added_it, ".", nn_info["method"])
                } else if (!("sNN" %in% n_key_added_it) & !is.null(nk)) {
                    net_name <- sub("_connectivities", "", nk)
                } else {
                    net_name <- paste0(net_type, ".", nn_info["method"])
                }
                parts <- strsplit(tn, "_")[[1]]
                nn_su <- parts[1]
                nn_ft <- parts[2]
                netObj <- createNearestNetObj(
                    name = net_name,
                    network = nn_network_igraph,
                    spat_unit = nn_su,
                    feat_type = nn_ft
                )

                gobject <- set_NearestNetwork(
                    gobject = gobject,
                    nn_network = netObj,
                    spat_unit = nn_su,
                    feat_type = nn_ft,
                    nn_network_to_use = net_type,
                    network_name = net_name,
                    set_defaults = FALSE
                )
            }
        }
    }

    # Attach hires image
    extracted_images <- extract_image(sdata)
    extract_image_names <- extract_image_names(sdata)

    raster_image_list <- lapply(extracted_images, terra::rast)
    large_image_list <- createGiottoLargeImageList(
        raster_image_list,
        names = extract_image_names
    )
    gobject <- addGiottoLargeImage(
        gobject = gobject, largeImages = large_image_list
    )

    # Extract points
    point_dict <- extract_points(sdata)
    for (ft in names(point_dict)) {
        gpoint_dt <- as.data.table(point_dict[[ft]])
        gpoint <- createGiottoPoints(
            gpoint_dt,
            feat_type = ft
        )
        gobject <- setGiotto(
            gobject,
            gpoint
        )
    }
    return(gobject)
}



##  Giotto to SpatialData####

#' @title Convert Giotto to SpatialData
#' @name giottoToSpatialData
#' @description Converts a Giotto object to SpatialData object
#'
#' @param gobject giotto object to be converted
#' @param spat_unit spatial unit which will be used in conversion
#' @param feat_type feature type which will be used in conversion
#' @param spot_radius radius of the spots
#' @param python_path path to python executable within a conda/miniconda
#' environment
#' @param env_name name of environment containing python_path executable
#' @param save_directory directory in which the SpatialData object will be saved
#'
#' @return SpatialData object saved on disk.
#' @details Function in beta. Converts and saves a Giotto object in SpatialData
#' format on disk.
#' @export

giottoToSpatialData <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        spot_radius = NULL,
        python_path = NULL,
        env_name = NULL,
        save_directory = NULL) {
    # Initialize reticulate
    instrs <- createGiottoInstructions(python_path = python_path)

    if (dir.exists(save_directory)) {
        stop("save_directory must be a non-existing directory")
    }

    # Check spatialdata dependencies
    package_check(
        pkg_name = "spatialdata",
        repository = "pip:spatialdata"
    )

    # Import sd2g, a python module for parsing SpatialData
    g2sd_path <- system.file("python", "g2sd.py", package = "GiottoConnect")
    reticulate::source_python(g2sd_path)

    spat_unit <- spat_unit %null% spatUnit(gobject)
    feat_type <- feat_type %null% featType(gobject)

    for (su in spat_unit) {
        wrap_msg("Spatial unit(s)", su, "will be used in conversion.")
    }
    for (ft in feat_type) {
        wrap_msg("Feature type(s)", ft, "will be used in conversion.")
    }

    # Create a temporary folder to hold anndata
    temp <- "temp_conversion_files/"
    # Delete temporary files and folders when done
    on.exit(
        {
            if (dir.exists(temp)) unlink(temp, recursive = TRUE)
        },
        add = TRUE
    )


    # First, convert Giotto object to AnnData using an existing function
    giottoToAnnData(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        python_path = python_path,
        env_name = env_name,
        save_directory = temp
    )

    # Extract GiottoImage only if an image exists
    images_exist <- NULL
    if (length(slot(gobject, "images")) > 0) {
        images_exist <- TRUE
        gimg_list <- slot(gobject, "images")
        for (i in seq_along(gimg_list)) {
            img_name <- slot(gimg_list[[i]], "name")
            writeGiottoLargeImage(
                giottoLargeImage = gimg_list[[i]],
                gobject = gobject,
                largeImage_name = img_name,
                filename = paste0(temp, img_name, ".png"),
                dataType = NULL,
                max_intensity = NULL,
                overwrite = TRUE,
                verbose = TRUE
            )
        }
    }

    # Extract polygons
    if (!is.null(list_spatial_info(gobject))) {
        dir.create(paste0(temp, "shapes"))
        for (su in spat_unit) {
            gpoly <- getPolygonInfo(gobject, polygon_name = su)
            gpoly_sf <- as.sf(gpoly)
            sf::st_write(gpoly_sf, paste0(temp, "shapes/", su, ".geojson"),
                delete_dsn = TRUE
            )
        }
    }

    # Extract points
    if (!is.null(list_feature_info(gobject))) {
        dir.create(paste0(temp, "points"))
        for (ft in feat_type) {
            gpoint <- getFeatureInfo(gobject, feat_type = ft)
            gpoint_dt <- as.data.table(gpoint, geom = "XY")
            fwrite(gpoint_dt, paste0(temp, "points/", ft, ".csv"),
                sep = ",", row.names = FALSE
            )
        }
    }

    # Create SpatialData object
    createSpatialData(temp, save_directory, images_exist)

    # copy over giotto metadata
    dir.create(file.path(save_directory, "giotto_meta"))
    file.copy(
        from = file.path(temp, "giotto_meta"),
        to = save_directory,
        recursive = TRUE
    )

    # Successful Conversion
    cat(
        "Giotto object has been converted and saved to SpatialData object at: ",
        save_directory, "\n"
    )
}
