test_that("checkGiottoMethod passes a conforming run and explains failures", {
    skip_if_not_installed("GiottoData")
    skip_if_not_installed("Giotto")
    options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
        giotto.no_python_warn = TRUE)
    g <- GiottoData::loadGiottoMini("visium")
    skip_if_not("runNeighborhoodPCA" %in% getNamespaceExports("Giotto"))
    g2 <- Giotto::runNeighborhoodPCA(g, feats_to_use = "hvf", npcs = 5,
        verbose = FALSE)
    # citations are tested below, independently of Giotto's version
    ok <- checkGiottoMethod(g, g2, slots = "dimension_reduction",
        param_class = "neighborhoodPCAParam", fn = "runNeighborhoodPCA",
        package = "Giotto", cite = FALSE)
    expect_true(all(ok$pass))
    expect_giotto_method(g, g2, slots = "dimension_reduction", cite = FALSE)

    # a param without methodCitation() fails unless cite = FALSE
    setClass("uncitedTestParam", contains = "analyzeParam",
        where = environment())
    g3 <- GiottoClass::update_giotto_params(g2, description = "_uncited",
        attachments = list(param = new("uncitedTestParam")), toplevel = 1)
    res <- checkGiottoMethod(g, g3, slots = "dimension_reduction")
    expect_false(res$pass[res$check == "param declares a citation"])
    expect_match(res$message[res$check == "param declares a citation"],
        "methodCitation")
    res <- checkGiottoMethod(g, g3, slots = "dimension_reduction", cite = FALSE)
    expect_false("param declares a citation" %in% res$check)

    bad <- GiottoClass::addCellMetadata(g, new_metadata = data.frame(
        cell_ID = GiottoClass::spatIDs(g), foo = NA), by_column = TRUE,
        column_cell_ID = "cell_ID")
    res <- checkGiottoMethod(g, bad, slots = "dimension_reduction")
    expect_false(res$pass[res$check == "only declared slots changed"])
    expect_false(res$pass[res$check == "param in history"])
    expect_match(res$message[res$check == "param in history"],
        "update_giotto_params")
})
