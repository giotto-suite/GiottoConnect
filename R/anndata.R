#' @import GiottoClass
#' @import GiottoUtils
NULL




## anndata object ####


#' @title Check Scanpy Installation
#' @name check_py_for_scanpy
#' @returns character
#' @description checks current python environment for scanpy 1.9.0
#' @keywords internal
check_py_for_scanpy <- function() {
    # test if scanpy is found
    module_test <- reticulate::py_module_available("scanpy")
    py_path <- reticulate::py_config()$python
    genv_in_use <- grepl(pattern = "giotto_env", x = py_path)

    if (module_test == FALSE && !genv_in_use) {
        stop(wrap_txt("scanpy python module is not installed:
            install in the environment or python path with:

            'pip install scanpy==1.9.0'

            Alternatively, install in the active python
            environment with reticulate:

            reticulate::py_install(packages = 'scanpy==1.9.0',
                                    pip = TRUE)
            \n
            ", errWidth = TRUE))
    } else if (module_test == FALSE && genv_in_use) {
        cat("Python module scanpy is required for conversion.
        Installing scanpy now in the Giotto Miniconda Environment.\n")

        conda_path <- reticulate::miniconda_path()
        py_ver <- reticulate::py_config()$version_string
        py_ver <- strsplit(py_ver, "|", fixed = TRUE)[[1]][1]
        py_ver <- gsub(" ", "", py_ver, fixed = TRUE)
        conda_full_path <- paste0(conda_path, "/", "bin/conda")
        full_envname <- paste0(conda_path, "/envs/giotto_env")

        reticulate::py_install(
            packages = "scanpy==1.9.0",
            envname = full_envname,
            method = "conda",
            conda = conda_full_path,
            pip = TRUE,
            python_version = py_ver
        )
    } else {
        cat("Required Python module scanpy has been previously installed.
            Proceeding with conversion.\n")
    }
}


#' @title Convert an AnnData .h5ad file to a Giotto object
#' @name anndataToGiotto
#' @description Convert an AnnData (.h5ad) file with the mapping that fits
#' its structure. The file is read directly (hdf5r, no python): it is
#' profiled, the matching bundled mapping is selected (or the one given) and
#' the conversion runs through [convertWithMapping()].
#' @param anndata_path path to the .h5ad file
#' @param mapping optional `giottoMapping`, mapping file or bundled mapping id.
#' By default the bundled mapping that applies to the file is used.
#' @param source `NULL` for an in-memory giotto object, or a directory for a
#' GiottoDisk-backed object
#' @param verbose logical
#' @returns giotto object; see [mappingProvenance()] for what was converted
#' and which losses apply
#' @seealso [findMapping()], [convertWithMapping()], [readbackCheck()]
#' @export
anndataToGiotto <- function(anndata_path, mapping = NULL, source = NULL, verbose = FALSE) {
    profile <- adapterProfile(.adapter_for("anndata"), anndata_path)
    if (is.null(mapping)) {
        found <- findMapping(profile, ecosystem = "anndata")
        if (found$status == "none") {
            stop(wrap_txt("No bundled mapping applies to this AnnData file:",
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
    convertWithMapping(anndata_path, mapping, profile = profile, source = source,
        verbose = verbose)
}


#' @title Convert Giotto to anndata
#' @name giottoToAnnData
#' @description Converts a Giotto object to a spatial anndata (e.g. scanpy)
#' .h5ad file
#' @param gobject giotto object to be converted
#' @param spat_unit spatial unit which will be used in conversion.
#' @param feat_type feature type which will be used in conversion.
#' @param python_path path to python executable within a conda/miniconda
#' environment
#' @param env_name name of environment containing python_path executable
#' @param save_directory directory in which the file will be saved.
#' @details Function in beta. Converts a Giotto object into .h5ad file(s).
#'
#' If there are multiple spatial units and/or feature types, but only
#' one spatial unit and/or feature type is specified, then only the
#' specified spatial unit and/or feature type will be used. If NULL,
#' by default, all spatial units will be used in conversion.
#'
#' If multiple spatial units or feature types are specified, multiple
#' AnnData object will be created and returned.
#'
#' This function will create .txt files which will record any `key_added`
#' parameters for networks. They are named after the corresponding spatial unit
#' and feature type pair.
#'
#' The save_directory will be created if it does not already exist.
#' The default save_directory is the working directory.
#' @returns vector containing .h5ad file path(s)
#' @export
giottoToAnnData <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        python_path = NULL,
        env_name = "giotto_env",
        save_directory = NULL) {
    # Check gobject
    invalid_obj <- !inherits(gobject, "giotto")
    if (is.null(gobject) || invalid_obj) {
        stop(wrap_msg("Please provide a valid Giotto Object for conversion."))
    }

    package_check(
        pkg_name = c("anndata", "scanpy"),
        repository = c("pip:anndata", "pip:scanpy")
    )

    # Python module import
    g2ad_path <- system.file("python", "g2ad.py", package = "GiottoConnect")
    reticulate::source_python(g2ad_path)
    if (!is.null(save_directory)) {
        save_directory <- paste0(save_directory, "/")
    }

    # Check directory, make it if it doesn't exist
    if (is.null(save_directory)) {
        save_directory <- paste0(getwd(), "/")
    } else if (!dir.exists(save_directory)) {
        warning(wrap_msg("Provided save directory not found. Creating save
                        directory at location:"))
        cat(save_directory)
        dir.create(save_directory, recursive = TRUE)
        if (dir.exists(save_directory)) {
            cat("Created directory", save_directory)
        } else {
            stop(wrap_msg("Unable to create directory. Please change the
                        provided path and try again."))
        }
    } else {
        wrap_msg(
            "Directory", save_directory,
            "found. The converted Giotto object will be saved here as an
                .h5ad file."
        )
    }

    # set up giotto meta folder
    gmdir <- file.path(save_directory, "giotto_meta")
    if (!dir.exists(gmdir)) dir.create(gmdir)

    # Expresion
    expr_dt <- list_expression(gobject)
    # ID spat_unit and feat_type if not already provided.
    if (is.null(spat_unit) && is.null(feat_type)) {
        spat_unit <- unique(expr_dt$spat_unit)
        feat_type <- unique(expr_dt$feat_type)
    } else if (is.null(spat_unit) && !is.null(feat_type)) {
        spat_unit <- unique(expr_dt$spat_unit)
    } else if (!is.null(spat_unit) && is.null(feat_type)) {
        feat_type <- unique(expr_dt$feat_type)
    }

    for (su in spat_unit) {
        wrap_msg("Spatial unit(s)", su, "will be used in conversion.")
    }
    for (ft in feat_type) {
        wrap_msg("Feature type(s)", ft, "will be used in conversion.")
    }

    # Iterate through spat_unit and feat_type to pull out expression data.
    # By default, the raw expression in the slot of the first spatial unit
    # and first feature type (if multiple are present) will be transferred to
    # the AnnData.anndata.X slot
    # Any other expression data will be inserted into AnnData.anndata.layers
    # By default, layer names are formed by
    # "'spatial_unit'_'feature_type'_'value'"

    adata <- NULL # scope
    su_ft_length <- 0

    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            su_ft_length <- su_ft_length + 1
        }
    }


    adata_list <- lapply(seq_len(su_ft_length), function(i) adata)
    adata_pos <- 1

    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            expr_names <- list_expression_names(
                gobject = gobject,
                spat_unit = su,
                feat_type = ft
            )

            for (en in expr_names) {
                if (en == "raw") {
                    raw_x <- getExpression(
                        gobject = gobject,
                        values = en,
                        spat_unit = su,
                        feat_type = ft,
                        output = "matrix"
                    )

                    adata <- ad_obj(x = raw_x)
                } else {
                    ad_layer_name <- paste0(su, "_", ft, "_", en)

                    x <- getExpression(
                        gobject = gobject,
                        values = en,
                        spat_unit = su,
                        feat_type = ft,
                        output = "matrix"
                    )

                    if ("dgeMatrix" %in% class(x)) x <- methods::as(x, "array")

                    adata <- set_adg_layer_data(
                        adata = adata,
                        lay = x,
                        lay_name = ad_layer_name
                    )
                }
            }
            adata_list[[adata_pos]] <- adata
            adata_pos <- adata_pos + 1
            adata <- NULL
        }
    }
    # Reset indexing variable
    adata_pos <- 1

    # Spatial Locations
    for (su in spat_unit) {
        for (ft_ in names(gobject@expression[[su]])) {
            sl <- get_spatial_locations(
                gobject = gobject,
                output = "data.table",
                spat_unit = su
            )
            n_col_sl <- dim(sl)[2]

            # preallocate data.table params
            sdimx <- sdimy <- sdimz <- NULL

            if (n_col_sl == 3) {
                sl <- sl[, .(sdimx, sdimy)]
            } else {
                sl <- sl[, .(sdimx, sdimy, sdimz)]
            }
            adata <- adata_list[[adata_pos]]
            adata <- set_adg_spat_locs(
                adata = adata,
                spat_locs = sl
            )
            adata_pos <- adata_pos + 1
        }
    }
    # Reset indexing variable
    adata_pos <- 1

    # Spatial Info

    # Cell Metadata
    # Feat Metadata
    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            cmeta <- getCellMetadata(
                gobject = gobject,
                spat_unit = su,
                feat_type = ft,
                output = "data.table",
                set_defaults = FALSE
            )

            fm <- getFeatureMetadata(
                gobject = gobject,
                spat_unit = su,
                feat_type = ft,
                output = "data.table",
                set_defaults = FALSE
            )

            adata_list[[adata_pos]] <- set_adg_metadata(
                adata = adata_list[[adata_pos]],
                cell_meta = cmeta,
                feat_meta = fm
            )

            adata_pos <- adata_pos + 1
        }
    }
    # Reset indexing variable
    adata_pos <- 1

    # Dimension Reductions

    # error hanldling wrapper to get_dimReduction
    try_get_dimReduction <- function(
        gobject,
        spat_unit,
        feat_type,
        reduction,
        reduction_method,
        name,
        output,
        set_defaults) {
        tryCatch(
            {
                dim_red <- get_dimReduction(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    reduction = reduction,
                    reduction_method = reduction_method,
                    name = name,
                    output = output,
                    set_defaults = set_defaults
                )
                return(dim_red)
            },
            error = function(e) {
                return(NULL)
            }
        )
    }

    ## PCA

    # pca on feats not supported by anndata because of dimensionality
    # agreement reqs
    reduction_list <- list_giotto_data(gobject, "dimension_reduction")
    dim_red <- NULL
    if (!is.null(reduction_list) && is.data.frame(reduction_list)) {
        for (i in seq_len(nrow(reduction_list))) {
            ro <- reduction_list[[i, 1]]
            su <- reduction_list[[i, 2]]
            ft <- reduction_list[[i, 3]]
            method <- reduction_list[[i, 4]]
            name <- reduction_list[[i, 5]]

            if (ro != "cells") {
                warning("AnnData does not support storing PCA by features.
                        Skipping PCA data conversion.")
                break
            }

            if (ft != "rna") name <- paste0(ft, ".pca")
            dim_red <- try_get_dimReduction(
                gobject = gobject,
                spat_unit = su,
                feat_type = ft,
                reduction = ro,
                reduction_method = method,
                name = name,
                output = "dimObj",
                set_defaults = FALSE
            )

            if (is.null(dim_red)) {
                adata_pos <- adata_pos + 1
                next
            }

            if (method == "pca") {
                pca_coord <- dim_red[]
                pca_loadings <- data.table(dim_red@misc$loadings)
                feats_used <- dimnames(dim_red@misc$loadings)[[1]]
                evs <- dim_red@misc$eigenvalues

                adata_list[[adata_pos]] <- set_adg_pca(
                    adata = adata_list[[adata_pos]],
                    pca_coord = pca_coord,
                    loadings = pca_loadings,
                    eigenv = evs,
                    feats_used = feats_used,
                    pca_name = name
                )
                wrap_msg(sprintf("\n%s converted.\n", name))
                adata_pos <- adata_pos + 1
            }

            if (method == "umap") {
                umap_data <- dim_red[]
                adata_list[[adata_pos]] <- set_adg_umap(
                    adata = adata_list[[adata_pos]],
                    umap_data = umap_data,
                    umap_name = name
                )
                wrap_msg(sprintf("\n%s converted.\n", name))
                adata_pos <- adata_pos + 1
            }

            if (method == "tsne") {
                tsne_data <- dim_red[]
                adata_list[[adata_pos]] <- set_adg_tsne(
                    adata = adata_list[[adata_pos]],
                    tsne_data = tsne_data,
                    tsne_name = name
                )
                wrap_msg(sprintf("\n%s converted.\n", name))
                adata_pos <- adata_pos + 1
            }
            adata_pos <- 1
        }
    }


    # Nearest Neighbor Network

    # error hanldling wrapper to get_NearestNetwork
    try_get_NN <- function(
        gobject,
        spat_unit,
        feat_type,
        nn_network_to_use,
        network_name,
        output,
        set_defaults) {
        tryCatch(
            {
                nearest_net <- get_NearestNetwork(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    nn_network_to_use = nn_network_to_use,
                    network_name = network_name,
                    output = output,
                    set_defaults = set_defaults
                )
                return(nearest_net)
            },
            error = function(e) {
                return(NULL)
            }
        )
    }

    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            nn_network_to_use <- c("sNN", "kNN")
            for (nn_net_tu in nn_network_to_use) {
                network_name <- list_nearest_networks_names(
                    gobject = gobject,
                    spat_unit = su,
                    feat_type = ft,
                    nn_type = nn_net_tu
                )
                if (is.null(network_name)) {
                    next
                }
                for (n_name in network_name) {
                    gob_NN <- try_get_NN(
                        gobject = gobject,
                        spat_unit = su,
                        feat_type = ft,
                        nn_network_to_use = nn_net_tu,
                        network_name = n_name,
                        output = "nnNetObj",
                        set_defaults = FALSE
                    )

                    pidx <- grep("nn_network", names(gobject@parameters))
                    for (p in pidx) {
                        if (gobject@parameters[[p]]["type"] == nn_net_tu) {
                            kval <- gobject@parameters[[p]]["k"]
                            dim_red_used <- gobject@parameters[[p]][
                                "dim_red_to_use"
                            ]
                        }
                    }

                    df_gob_NN <- igraph::as_data_frame(gob_NN[])

                    adata_list[[adata_pos]] <- set_adg_nn(
                        adata = adata_list[[adata_pos]],
                        df_NN = df_gob_NN,
                        net_name = n_name,
                        n_neighbors = kval,
                        dim_red_used = dim_red_used
                    )
                }

                # Save NN keys to .uns['NN_keys']
                network_name <- network_name[!grepl("kNN.", network_name)]
                save_NN_keys(adata = adata_list[[adata_pos]], network_name = network_name)
            }
            adata_pos <- adata_pos + 1
        }
    }

    # Reset indexing variable
    adata_pos <- 1

    try_get_SN <- function(
        gobject,
        spat_unit,
        name,
        output,
        set_defaults,
        verbose) {
        tryCatch(
            {
                spatial_net <- get_spatialNetwork(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    name = name,
                    output = output,
                    set_defaults = set_defaults,
                    verbose = verbose
                )
                return(spatial_net)
            },
            error = function(e) {
                return(NULL)
            }
        )
    }

    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            # Spatial networks do not have a feature type slot.
            # Iterate through anyways to properly assign to anndata objects
            network_name <- list_spatial_networks_names(
                gobject = gobject,
                spat_unit = su
            )
            if (is.null(network_name)) {
                next
            }
            for (sn_name in network_name) {
                gob_SN <- try_get_SN(
                    gobject = gobject,
                    spat_unit = su,
                    name = sn_name,
                    output = "networkDT",
                    set_defaults = FALSE,
                    verbose = TRUE
                )

                pidx <- grep("spatial_network", names(gobject@parameters))
                for (p in pidx) {
                    if (gobject@parameters[[p]]["name of spatial network"] == sn_name) {
                        current_param <- gobject@parameters[[p]]
                        kval <- current_param["k neighbours"] # null in the case of Delaunay
                        maxdist <- current_param["maximum distance threshold"] # "auto" for Delaunay, "400 (string)" for sNN
                        dimused <- current_param["dimensions used"]
                    }
                }

                adata_list[[adata_pos]] <- set_adg_sn(
                    adata = adata_list[[adata_pos]],
                    df_SN = gob_SN,
                    net_name = sn_name,
                    n_neighbors = kval,
                    max_distance = maxdist,
                    dim_used = dimused
                )
            }

            # Save SN keys to .uns['SN_keys']
            if (length(network_name) != 0) {
                save_SN_keys(
                    adata = adata_list[[adata_pos]],
                    network_name = network_name
                )
            }
        }

        adata_pos <- adata_pos + 1
    }

    # Reset indexing variable
    adata_pos <- 1

    # Spatial Enrichment
    spat_enrich_list <- list_giotto_data(
        gobject = gobject,
        slot = "spatial_enrichment"
    )

    if (!is.null(spat_enrich_list) && is.data.frame(spat_enrich_list)) {
        for (i in seq_len(nrow(spat_enrich_list))) {
            se_su <- spat_enrich_list[i]$spat_unit
            se_ft <- spat_enrich_list[i]$feat_type
            se_name <- spat_enrich_list[i]$name

            se <- getSpatialEnrichment(
                gobject = gobject,
                spat_unit = se_su,
                feat_type = se_ft,
                name = se_name,
                output = "data.table"
            )
            adata_list[[adata_pos]] <- set_adg_spat_enrich(
                adata = adata_list[[adata_pos]],
                enrichment = se,
                name = se_name
            )
        }
        save_SE_keys(
            adata = adata_list[[adata_pos]],
            enrichment_name = spat_enrich_list$name
        )
    }

    # Write AnnData object to .h5ad file
    # Verify it exists, and return upon success
    fname_list <- lapply(seq_len(su_ft_length), function(i) NULL)
    wrap_msg("\n")
    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            adata <- adata_list[[adata_pos]]
            path_adata <- write_ad_h5ad(
                adata = adata,
                save_directory = save_directory,
                spat_unit = su,
                feat_type = ft
            )
            if (!is.null(path_adata)) {
                wrap_msg(
                    "Spatial unit", su, "and feature type",
                    ft, "converted to:"
                )
                wrap_msg(path_adata)
                fname_list[[adata_pos]] <- path_adata
                adata_pos <- adata_pos + 1
            } else {
                wrap_msg(
                    "Unable to convert spatial unit feature type pair",
                    su, ft
                )
                stop(wrap_msg("File writing error. Please try again."))
            }
        }
    }

    wrap_msg("\nGiotto object successfully converted to .h5ad file(s)\n")

    return(fname_list)
}
