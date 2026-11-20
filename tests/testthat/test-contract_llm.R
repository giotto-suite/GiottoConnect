llm_cases <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            env <- new.env()
            sys.source(contract_script("build_fixtures.R"), envir = env)
            cache <<- suppressWarnings(env$build_llm_cases())
        }
        cache
    }
})

# a backend whose replies are fixed: patches as list(item = value, ...)
mock_backend <- function(choices, confidence = "high") {
    reply <- lapply(names(choices), function(it) {
        list(item = it, value = choices[[it]], confidence = confidence, rationale = "mock")
    })
    structure(list(provider = "ollama", model = "mock", host = "", max_tokens = 1L,
        reply = reply), class = "llmBackend")
}
with_mock_llm <- function(code) {
    local_mocked_bindings(
        .llm_ollama = function(request, backend) list(patches = backend$reply, input_tokens = 1L,
            output_tokens = 1L),
        .package = "GiottoConnect"
    )
    code
}

test_that("requests offer closed choices per item", {
    skip_if_not_installed("SeuratObject")
    s <- llm_cases()$L1$object
    m <- readMapping("seurat-assay5")
    p <- profileObject(s)
    r <- repairMapping(m, s)
    items <- GiottoConnect:::.llm_items(m, r$unresolved, list(), r$candidates, p,
        GiottoConnect:::.bind_params(m, p))
    req <- GiottoConnect:::.llm_request(m, items, p, r$diff)
    any_of <- req$schema$properties$patches$items$anyOf
    expect_length(any_of, 1L)
    expect_equal(any_of[[1]]$properties$item$const, "expr_raw.from")
    expect_true("misc/stash" %in% unlist(any_of[[1]]$properties$value$enum))
    expect_false(any_of[[1]]$additionalProperties)
    expect_match(req$prompt, "misc/stash")
    # no data values in the prompt
    expect_false(grepl("[0-9]+\\.[0-9]{3,}", req$prompt))
})

test_that("replies outside the offered options or with low confidence are dropped", {
    req <- list(items = list(list(item = "a.from", options = c("x/y", "x/z"))), prompt = "")
    raw <- list(
        list(item = "a.from", value = "x/q", confidence = "high", rationale = "off-list"),
        list(item = "b.from", value = "x/y", confidence = "high", rationale = "unknown item"),
        list(item = "a.from", value = "x/z", confidence = "low", rationale = "unsure"),
        list(item = "a.from", value = "x/y", confidence = "medium", rationale = strrep("r", 400))
    )
    kept <- GiottoConnect:::.llm_check(raw, req)
    expect_equal(nrow(kept), 1L)
    expect_equal(kept$value, "x/y")
    expect_equal(nchar(kept$rationale), 300L)
})

test_that("an LLM choice is applied, validated and read back", {
    skip_if_not_installed("SeuratObject")
    s <- llm_cases()$L1$object
    r <- with_mock_llm(repairMapping("seurat-assay5", s,
        llm = mock_backend(c(expr_raw.from = "misc/stash"))))
    expect_equal(r$status, "repaired")
    expect_equal(r$mapping$mapping$status, "patched-unreviewed")
    expect_equal(r$mapping$mapping$history[[1]]$by, "llm")
    expect_equal(r$llm_log$kept, 1L)
    expect_true(attr(r$readback, "pass"))
    # the patch reads by path and validation flags it for review
    expect_equal(r$mapping$source_nodes$expr_raw_patched$get$fn, "path")
    v <- validateMapping(r$mapping, profileObject(s))
    expect_true(any(v$level == "WARN" & grepl("read by path", v$message)))
})

test_that("param choices switch a part to its declared alternative", {
    skip_if_not_installed("SeuratObject")
    s <- llm_cases()$L3$object
    r <- with_mock_llm(repairMapping("seurat-assay5", s, llm = mock_backend(c(
        locs.locs_from_meta.x_col = "CenterX_global_px",
        locs.locs_from_meta.y_col = "CenterY_global_px"))))
    expect_equal(r$status, "repaired")
    g <- seuratToGiotto(s, mapping = r$mapping)
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    sl <- sl[match(colnames(s), sl$cell_ID), ]
    expect_equal(sl$sdimx, s$CenterX_global_px, ignore_attr = TRUE)
})

test_that("a wrong choice is caught by the semantic gate and handed off", {
    skip_if_not_installed("SeuratObject")
    s <- llm_cases()$L5$object
    f <- tempfile(fileext = ".md")
    r <- with_mock_llm(repairMapping("seurat-assay5", s,
        llm = mock_backend(c(expr_raw.from = "graphs/Xenium_nn")), request_file = f))
    expect_equal(r$status, "failed")
    expect_true(any(grepl("feat_ids|features", r$validation$message)) ||
        !isTRUE(attr(r$readback, "pass")))
    expect_true(file.exists(f))
    md <- readLines(f)
    expect_true(any(grepl("^# Mapping repair request", md)))
    expect_true(any(grepl("Ranked candidates", md)))
})

test_that("the baseline only picks a clearly leading candidate", {
    req <- list(prompt = "", items = list(
        list(item = "a.from", kind = "from", options = c("p", "q"), scores = c(p = 4, q = 3)),
        list(item = "b.from", kind = "from", options = c("p", "q"), scores = c(p = 4, q = 3.8))))
    out <- llmPatch(req, llmBackend("baseline"))
    expect_equal(out$patches$item, "a.from")
    expect_equal(out$log$provider, "baseline")
})

test_that("backends read settings from the environment and never need a key to exist", {
    b <- llmBackend("anthropic", model = "")
    expect_equal(b$model, "claude-haiku-4-5")
    expect_error(llmBackend("ollama", model = ""), "needs a model")
    withr_env <- Sys.getenv(c("ANTHROPIC_API_KEY", "ANTHROPIC_AUTH_TOKEN"))
    on.exit(do.call(Sys.setenv, as.list(withr_env)), add = TRUE)
    Sys.unsetenv(c("ANTHROPIC_API_KEY", "ANTHROPIC_AUTH_TOKEN"))
    out <- llmPatch(list(items = list(), prompt = "x", schema = list()), b)
    expect_match(out$log$error, "no Anthropic credentials")
    expect_equal(out$log$kept, 0L)
})

test_that("conversion never reaches the LLM backends", {
    skip_if_not_installed("SeuratObject")
    local_mocked_bindings(
        .llm_ollama = function(...) stop("LLM called"),
        .llm_anthropic = function(...) stop("LLM called"),
        .package = "GiottoConnect"
    )
    expect_s4_class(seuratToGiotto(seurat_fixture()), "giotto")
    engine <- c(deparse(convertWithMapping), deparse(GiottoConnect:::.prepare_parts))
    expect_false(any(grepl("llm", engine, ignore.case = TRUE)))
})

test_that("live LLM repair (opt-in)", {
    skip_if(Sys.getenv("GIOTTO_LLM_PROVIDER") == "")
    skip_if_not_installed("SeuratObject")
    r <- repairMapping("seurat-assay5", llm_cases()$L1$object, llm = llmBackend())
    expect_true(r$status %in% c("repaired", "failed"))
    expect_equal(nrow(r$llm_log), 1L)
})
