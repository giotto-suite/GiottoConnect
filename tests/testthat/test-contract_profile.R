node <- function(p, path) p$nodes[p$nodes$path == path, , drop = FALSE]

test_that("profileObject describes generic R structures", {
    x <- list(
        counts = matrix(c(1L, 0L, 2L, 3L), 2),
        meta = data.frame(a = 1:2, b = c(0.5, 1.5), c = c("u", "v")),
        polys = list(list(x = 1), list(x = 2), list(x = 3)),
        env = list2env(list(k = "v")),
        f = function() NULL
    )
    p <- profileObject(x)
    expect_s3_class(p, "giottoProfile")
    expect_equal(node(p, "counts")$kind, "dense_matrix")
    expect_equal(node(p, "counts")$shape, "2x2")
    expect_match(node(p, "counts")$role_hints, "integer_valued")
    expect_equal(node(p, "meta")$kind, "data.frame")
    expect_match(node(p, "meta")$role_hints, "two_numeric_cols")
    expect_equal(node(p, "meta/c")$dtype, "character")
    # unnamed lists are represented by their first element
    expect_true("polys/*/x" %in% p$nodes$path)
    expect_false(any(grepl("^polys/2", p$nodes$path)))
    expect_equal(node(p, "env/k")$kind, "vector")
    expect_equal(node(p, "f")$kind, "function")
})

test_that("cells are inferred from the shallowest matching table", {
    # more features than cells, with feature metadata nested deeper
    x <- list(
        assays = list(rna = list(
            counts = matrix(rep(0:2, 100L), 30, 10),
            meta.data = data.frame(f = seq_len(30))
        )),
        meta.data = data.frame(c = seq_len(10))
    )
    p <- profileObject(x)
    expect_equal(p$n_cells, 10L)
    expect_equal(p$n_feat, 30L)
    expect_equal(profileObject(x, n_cells = 30L, n_feat = 10L)$size_from, "user")

    # Giotto-style: cell and feature tables at the same depth
    y <- list(
        expression = list(raw = matrix(rep(0:2, 105L), 15, 21)),
        cell_metadata = list(meta = data.frame(c = seq_len(21))),
        feat_metadata = list(meta = data.frame(f = seq_len(15)))
    )
    py <- profileObject(y)
    expect_equal(py$n_cells, 21L)
    expect_equal(py$n_feat, 15L)
})

test_that("profileObject respects depth and child limits", {
    deep <- list(a = list(b = list(c = list(d = 1))))
    p <- profileObject(deep, max_depth = 2L)
    expect_false("a/b/c" %in% p$nodes$path)
    expect_true(node(p, "a/b")$truncated)

    wide <- stats::setNames(as.list(1:10), letters[1:10])
    p <- profileObject(list(w = wide), max_children = 3L)
    expect_true(node(p, "w")$truncated)
    expect_equal(sum(startsWith(p$nodes$path, "w/")), 3L)
})

test_that("profileObject describes a Seurat v5 object", {
    skip_if_not_installed("SeuratObject")
    p <- profileObject(seurat_fixture())
    expect_equal(p$source$class, "Seurat")
    expect_equal(p$source$package, "SeuratObject")
    expect_equal(p$n_cells, 200L)
    expect_equal(p$n_feat, 50L)
    expect_equal(p$size_from, "assays/Xenium/layers/counts")

    counts <- node(p, "assays/Xenium/layers/counts")
    expect_equal(counts$kind, "sparse_matrix")
    expect_equal(counts$shape, "50x200")
    expect_match(counts$role_hints, "integer_valued")
    expect_match(counts$role_hints, "cols_eq_ncells")
    # Assay5 layers carry no dimnames; names live in LogMap slots
    expect_false(grepl("has_dimnames", counts$role_hints))
    expect_equal(node(p, "assays/Xenium/cells")$class, "LogMap")

    expect_match(node(p, "graphs/Xenium_nn")$role_hints, "binary_valued")
    expect_false(grepl("binary_valued", node(p, "graphs/Xenium_snn")$role_hints))
    expect_equal(node(p, "images/fov/boundaries/segmentation")$class, "Segmentation")
    expect_true(any(grepl("requires package sp", p$notes)))
    expect_output(print(p), "200 cells x 50 features")
})

test_that("profiling a Seurat object does not load SeuratObject", {
    skip_if_not_installed("SeuratObject")
    skip_on_cran()
    rds <- tempfile(fileext = ".rds")
    saveRDS(seurat_fixture(), rds)
    root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
    from_source <- file.exists(file.path(root, "R", "contract_profile.R"))
    load_pkg <- if (from_source && requireNamespace("pkgload", quietly = TRUE)) {
        sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE))", deparse(root))
    } else {
        "suppressMessages(library(GiottoConnect))"
    }
    out <- run_fresh_r(c(
        load_pkg,
        "before <- 'SeuratObject' %in% loadedNamespaces()",
        sprintf("s <- readRDS(%s)", deparse(rds)),
        "p <- GiottoConnect::profileObject(s)",
        "cat('RESULT', before, 'SeuratObject' %in% loadedNamespaces(), p$n_cells, '\\n')"
    ))
    res <- grep("^RESULT", out, value = TRUE)
    expect_length(res, 1L)
    fields <- strsplit(res, " ")[[1]]
    if (fields[2] == "TRUE") skip("SeuratObject is loaded by GiottoConnect's dependencies")
    expect_equal(fields[3], "FALSE")
    expect_equal(fields[4], "200")
})

test_that("fingerprint ignores sizes and templates parameters", {
    skip_if_not_installed("SeuratObject")
    env <- new.env()
    sys.source(contract_script("build_fixtures.R"), envir = env)
    a <- profileObject(seurat_fixture())
    b <- profileObject(env$build_seurat_fixture(n_cells = 120L, n_feat = 30L, seed = 7L))
    c <- profileObject(env$build_seurat_fixture(assay = "Spatial"))
    paths <- c("assays/{assay}/layers/counts", "graphs/{assay}_nn",
        "images/{fov}/boundaries/centroids")

    expect_equal(fingerprint(a), fingerprint(b))
    expect_false(fingerprint(a) == fingerprint(c))
    expect_equal(
        fingerprint(a, paths, params = list(assay = "Xenium", fov = "fov")),
        fingerprint(c, paths, params = list(assay = "Spatial", fov = "fov"))
    )
    expect_equal(
        fingerprint(a, scope = "full", params = list(assay = "Xenium")),
        fingerprint(c, scope = "full", params = list(assay = "Spatial"))
    )
})

test_that("fingerprint over paths ignores unrelated changes but sees missing nodes", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    params <- list(assay = "Xenium", fov = "fov")
    paths <- c("assays/{assay}/layers/counts", "reductions/pca/cell.embeddings")
    fp <- fingerprint(profileObject(s), paths, params)

    s2 <- s
    s2@misc$note <- "unrelated"
    s2@reductions$umap <- NULL
    expect_equal(fingerprint(profileObject(s2), paths, params), fp)
    expect_false(fingerprint(profileObject(s2)) == fingerprint(profileObject(s)))

    s3 <- s
    s3@reductions$pca <- NULL
    expect_false(fingerprint(profileObject(s3), paths, params) == fp)
})

test_that("diffProfiles classifies structural changes", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    old <- profileObject(s)

    s2 <- s
    lyr <- s2@assays$Xenium@layers
    names(lyr)[names(lyr) == "counts"] <- "raw_counts"
    s2@assays$Xenium@layers <- lyr
    s2@misc$old_umap <- s2@reductions$umap
    s2@reductions$umap <- NULL
    s2@reductions$tsne <- s2@reductions$pca
    s2@meta.data$nFeature_Xenium <- as.double(s2@meta.data$nFeature_Xenium)
    d <- diffProfiles(old, profileObject(s2))

    expect_s3_class(d, "giottoProfileDiff")
    row <- function(ch) d[d$change == ch, , drop = FALSE]
    expect_equal(row("renamed_key")$path, "assays/Xenium/layers/counts")
    expect_equal(row("renamed_key")$new_path, "assays/Xenium/layers/raw_counts")
    expect_equal(row("moved")$path, "reductions/umap")
    expect_equal(row("moved")$new_path, "misc/old_umap")
    expect_equal(row("added")$new_path, "reductions/tsne")
    expect_equal(row("dtype_changed")$path, "meta.data/nFeature_Xenium")
    expect_output(print(d), "change")

    expect_equal(nrow(diffProfiles(old, profileObject(s))), 0L)
})

test_that("diffProfiles uses params to absorb dataset-specific names", {
    skip_if_not_installed("SeuratObject")
    env <- new.env()
    sys.source(contract_script("build_fixtures.R"), envir = env)
    a <- profileObject(seurat_fixture())
    b <- profileObject(env$build_seurat_fixture(assay = "Spatial"))
    expect_true(nrow(diffProfiles(a, b)) > 0L)
    d <- diffProfiles(a, b,
        params_old = list(assay = "Xenium"), params_new = list(assay = "Spatial")
    )
    expect_equal(nrow(d), 0L)
})

test_that("fingerprints do not depend on the collation locale", {
    p <- profileObject(list(B = 1, a = list(`{x}` = 2, Z.y = 3)))
    old <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
    Sys.setlocale("LC_COLLATE", "C")
    fp_c <- fingerprint(p)
    if (nzchar(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))) {
        expect_equal(fingerprint(p), fp_c)
    } else {
        skip("en_US.UTF-8 collation not available")
    }
})
