library(testthat)

giotto_obj <- GiottoData::loadGiottoMini("visium")
seurat_obj <- giottoToSeurat(giotto_obj)
giotto_obj_roundtrip <- seuratToGiotto(seurat_obj)

test_that("giottoToSeurat handles basic conversion", {
    expect_s4_class(seurat_obj, "Seurat")
    expect_true(!is.null(seurat_obj@meta.data))
    expect_equal(
        nrow(seurat_obj@meta.data),
        nrow(giotto_obj@cell_metadata$cell$rna@metaDT)
    )
    # spatial locations and spot polygons become an FOV
    expect_s4_class(seurat_obj[["fov"]], "FOV")
})

test_that("Seurat to Giotto Conversion Works", {
    expect_s4_class(giotto_obj_roundtrip, "giotto")
})

test_that("Assay names are converted correctly (no spaces)", {
    assay_names <- names(seurat_obj@assays)
    expect_false(any(grepl(" ", assay_names)))
})

test_that("Data is consistent after roundtrip conversion", {
    expect_equal(
        dim(giotto_obj@expression$cell$rna$raw[]),
        dim(giotto_obj_roundtrip@expression$cell$rna$raw[])
    )
    expect_equal(
        nrow(giotto_obj@cell_metadata$cell$rna@metaDT),
        nrow(giotto_obj_roundtrip@cell_metadata$cell$rna@metaDT)
    )
})

test_that("Feature names are transferred correctly", {
    expect_setequal(
        rownames(seurat_obj[["rna"]]),
        gsub("_", "-", GiottoClass::featIDs(giotto_obj))
    )
})

