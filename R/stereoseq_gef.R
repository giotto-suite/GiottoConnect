#' @import GiottoClass
#' @import GiottoUtils
NULL

## sp classes ####



## gef object ####
#' @title Convert gef to Giotto
#' @name gefToGiotto
#' @description Converts .gef file (output stereo-seq pipeline) into
#' giotto subcellular object
#'
#' @param gef_file path to .gef file
#' @param bin_size bin size to select from .gef file
#' @param gene_column (optional) character. Which column contains gene names
#' within the geneExp information.
#' @param h5_file name to create and on-disk HDF5 file
#' @param verbose be verbose
#'
#' @details Function in beta. Converts .gef object to Giotto object.
#'
#' There are six possible choices for bin_size:
#' bin1, bin10, bin20, bin50, bin100, bin200.
#'
#' See SAW pipeline for additional information about the gef file.
#' @returns giotto object
#' @export

gefToGiotto <- function(gef_file,
    bin_size = "bin100",
    gene_column = NULL,
    verbose = FALSE,
    h5_file = NULL) {
    # data.table vars
    genes <- gene_idx <- x <- y <- sdimx <- sdimy <- cell_ID <- bin_ID <-
        count <- i.bin_ID <- NULL

    # package check
    package_check(pkg_name = "rhdf5", repository = "Bioc")
    if (!file.exists(gef_file)) stop("File path to .gef file does not exist")

    # check if proper bin_size is selected. These are determined in SAW pipeline
    wrap_msg("1. gefToGiotto() begin... \n")
    bin_size_options <- c("bin1", "bin10", "bin20", "bin50", "bin100", "bin200")
    if (!(bin_size %in% bin_size_options)) {
        stop("Please select valid bin size, see ?gefToGiotto for details.")
    }

    # 1. read .gef file at specific bin size
    geneExpData <- rhdf5::h5read(
        file = gef_file,
        name = paste0("geneExp/", bin_size)
    )

    exprDT <- data.table::as.data.table(geneExpData[["expression"]])
    exprDT[, count := lapply(.SD, as.integer), .SDcols = "count"]
    data.table::setorder(exprDT, x, y) # sort by x, y coords (ascending)
    geneDT <- data.table::as.data.table(geneExpData[["gene"]])

    # gene_column selection
    name_to_replace <- gene_column %null% "gene"
    if (name_to_replace %in% colnames(geneDT)) {
        data.table::setnames(geneDT, old = name_to_replace, new = "geneName")
    }

    # process duplicated gene symbol
    if (any(duplicated(geneDT$geneName))) {
        duplicated_genes <- unique(geneDT$geneName[duplicated(geneDT$geneName)])
        cat("Ops!!! Duplicated_genes,processing:sum(count),mean(offset)")
        # merge
        for (gene in duplicated_genes) {
            # indices
            idx <- which(geneDT$geneName == gene)

            #
            cat("Processing gene:", gene, "\n")
            cat("Original count values for", gene, ":", geneDT$count[idx], "\n")

            # update
            geneDT$count[idx[1]] <- sum(geneDT$count[idx]) # 对重复的 count 求和
            geneDT$offset[idx[1]] <- mean(geneDT$offset[idx]) # 对重复的 offset 求平均

            #
            cat("Updated count for", gene, ":", geneDT$count[idx[1]], "\n")
            cat("Updated offset for", gene, ":", geneDT$offset[idx[1]], "\n")

            #
            geneDT <- geneDT[-idx[-1], ]

            #
            cat("Deleted duplicate entries for gene:", gene, "\n\n")
        }
    }

    if (isTRUE(verbose)) wrap_msg("finished reading in .gef", bin_size, "\n")

    # 2. create spatial locations
    if (isTRUE(verbose)) wrap_msg("2. create spatial_locations... \n")
    cell_locations <- unique(exprDT[, c("x", "y")], by = c("x", "y"))
    cell_locations[, bin_ID := as.factor(seq_len(nrow(cell_locations)))]
    cell_locations[, cell_ID := paste0("cell_", bin_ID)]
    data.table::setcolorder(cell_locations, c("x", "y", "cell_ID", "bin_ID"))
    # ensure first non-numerical col is cell_ID
    if (isTRUE(verbose)) wrap_msg(nrow(cell_locations), " bins in total \n")
    if (isTRUE(verbose)) wrap_msg("finished spatial_locations \n")

    # 3. create expression matrix
    if (isTRUE(verbose)) wrap_msg("3. create expression matrix... \n")
    exprDT[, genes := as.character(rep(x = geneDT$geneName, geneDT$count))]
    exprDT[, gene_idx := as.integer(factor(exprDT$genes,
        levels = unique(exprDT$genes)
    ))]

    # merge on x,y and populate based on bin_ID values in cell_locations
    exprDT[cell_locations, cell_ID := i.bin_ID, on = .(x, y)]
    exprDT$cell_ID <- as.integer(exprDT$cell_ID)

    expMatrix <- Matrix::sparseMatrix(
        i = exprDT$gene_idx,
        j = exprDT$cell_ID,
        x = exprDT$count
    )

    colnames(expMatrix) <- cell_locations$cell_ID
    rownames(expMatrix) <- geneDT$geneName
    rm(exprDT)
    if (isTRUE(verbose)) wrap_msg("finished expression matrix")

    # 4. create minimal giotto object
    if (isTRUE(verbose)) wrap_msg("4. create giotto object... \n")
    stereo <- createGiottoObject(
        expression = expMatrix,
        spatial_locs = cell_locations,
        verbose = FALSE,
        h5_file = h5_file
    )
    if (isTRUE(verbose)) wrap_msg("finished giotto object... \n")

    wrap_msg("gefToGiotto() finished \n")
    return(stereo)
}
