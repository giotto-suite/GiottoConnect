# BayesSpace connection. The offline tests cover the param and the spot grid;
# the full run is opt-in (GIOTTOCONNECT_TOOL_TESTS).

test_that("bayesspaceParam checks and stores its options", {
    p <- bayesspaceParam(q = 4, platform = "ST", nrep = 1000, burn_in = 100)
    expect_s4_class(p, "bayesspaceParam")
    expect_identical(p$q, 4L)
    expect_identical(p$platform, "ST")
    expect_null(p$gamma)
    expect_error(bayesspaceParam())
    expect_error(bayesspaceParam(q = 3, nrep = 100, burn_in = 200))
    expect_error(bayesspaceParam(q = 3, model = "x"))
})

# Visium-like hexagonal lattice: rows 0..r, columns of the row's parity
hex_locs <- function(nr = 6, nc = 10, pitch = 138, swap = FALSE) {
    g <- expand.grid(array_row = 0:(nr - 1), array_col = 0:(2 * nc - 1))
    g <- g[(g$array_row + g$array_col) %% 2 == 0, ]
    x <- 500 + g$array_col * pitch / 2 + stats::rnorm(nrow(g), sd = 1)
    y <- -200 - g$array_row * pitch * sqrt(3) / 2 + stats::rnorm(nrow(g), sd = 1)
    d <- data.frame(cell_ID = paste0("s", seq_len(nrow(g))),
        sdimx = if (swap) y else x, sdimy = if (swap) x else y)
    list(locs = d, truth = g)
}

# neighbour pairs of a grid under BayesSpace's offsets
grid_pairs <- function(row, col, hex = TRUE) {
    off <- if (hex) cbind(c(-2, 2, -1, 1, -1, 1), c(0, 0, -1, -1, 1, 1)) else
        cbind(c(0, 1, 0, -1), c(-1, 0, 1, 0))
    key <- paste(col, row)
    out <- character()
    for (k in seq_len(nrow(off))) {
        j <- match(paste(col + off[k, 1], row + off[k, 2]), key)
        ok <- !is.na(j)
        out <- c(out, paste(which(ok), j[ok]))
    }
    sort(out)
}

test_that("the Visium grid is derived from coordinates in either orientation", {
    set.seed(1)
    for (swap in c(FALSE, TRUE)) {
        h <- hex_locs(swap = swap)
        gr <- GiottoConnect:::.bayesspace_grid(h$locs, "Visium")
        expect_identical(gr$cell_ID, h$locs$cell_ID)
        expect_identical(grid_pairs(gr$array_row, gr$array_col),
            grid_pairs(h$truth$array_row, h$truth$array_col))
    }
})

test_that("given grid positions are kept; square grids are supported", {
    gr <- GiottoConnect:::.bayesspace_grid(data.frame(cell_ID = c("a", "b"),
        array_row = c(1, 2), array_col = c(3, 4)))
    expect_identical(gr$array_col, c(3L, 4L))
    sq <- expand.grid(r = 0:4, c = 0:4)
    locs <- data.frame(cell_ID = paste0("c", seq_len(nrow(sq))),
        sdimx = sq$c * 50, sdimy = sq$r * 50)
    gr <- GiottoConnect:::.bayesspace_grid(locs, "ST")
    expect_identical(grid_pairs(gr$array_row, gr$array_col, hex = FALSE),
        grid_pairs(sq$r, sq$c, hex = FALSE))
    expect_error(GiottoConnect:::.bayesspace_grid(data.frame(cell_ID = "a",
        x = 1)), "sdimx")
})

test_that("runBayesSpace writes clusters and BayesSpace's PCA", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("BayesSpace")
    skip_if_not_installed("GiottoData")
    options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
        giotto.no_python_warn = TRUE)
    g <- GiottoData::loadGiottoMini("visium")
    g <- runBayesSpace(g, param = bayesspaceParam(q = 5, nrep = 1000,
        burn_in = 100))
    md <- GiottoClass::pDataDT(g)
    expect_true(all(c("bayesspace", "bayesspace_init") %in% names(md)))
    expect_false(anyNA(md$bayesspace))
    expect_true(all(as.integer(as.character(md$bayesspace)) %in% 1:5))
    pc <- GiottoClass::getDimReduction(g, reduction_method = "pca",
        name = "bayesspace_pca", output = "matrix")
    expect_identical(dim(pc), c(nrow(md), 15L))
    h <- g@parameters[[length(g@parameters)]]
    expect_s4_class(attr(h, "attachments")$param, "bayesspaceParam")

    g2 <- runBayesSpace(g, param = bayesspaceParam(q = 5, nrep = 1000,
        burn_in = 100), dim_reduction_name = "pca", name = "bs_giotto_pca")
    expect_true("bs_giotto_pca" %in% names(GiottoClass::pDataDT(g2)))
})
