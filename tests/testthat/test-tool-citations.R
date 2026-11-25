test_that("every external tool declares its paper", {
    params <- list(
        banksyParam(), bayesspaceParam(q = 3), cellchatParam(),
        nnsvgParam(), ontracParam()
    )
    for (p in params) {
        cit <- methodCitation(p)
        expect_s3_class(cit, "bibentry")
        expect_true(all(nzchar(unlist(cit$doi))), label = class(p))
    }
})

test_that("a recorded tool param is listed by objCitations()", {
    old <- options(giotto.cite = FALSE)
    on.exit(options(old))
    g <- GiottoClass::update_giotto_params(giotto(), description = "_nnSVG",
        attachments = list(param = nnsvgParam()), toplevel = 1)
    expect_match(objCitations(g), "nnSVG")
})
