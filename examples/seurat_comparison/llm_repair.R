# LLM repair benchmark: cases L1-L5 (drift the rules cannot resolve alone)
# against the no-LLM baseline and the configured LLM backends.
# Backends: always "baseline"; Ollama when OLLAMA_MODEL is set (default
# qwen3.8:27b if the server answers); Anthropic when ANTHROPIC_API_KEY is set.
# Run from the GiottoConnect package root.
options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE, width = 160)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
env <- new.env()
sys.source(file.path("inst", "scripts", "build_fixtures.R"), envir = env)
cases <- suppressWarnings(env$build_llm_cases())

backends <- list(baseline = llmBackend("baseline"))
ollama_model <- Sys.getenv("OLLAMA_MODEL", "qwen3.8:27b")
if (!identical(Sys.getenv("SKIP_OLLAMA"), "1") && !inherits(try(httr2::req_perform(
    httr2::req_timeout(httr2::request("http://localhost:11434/api/tags"), 3)), silent = TRUE), "try-error")) {
    backends$ollama <- llmBackend("ollama", model = ollama_model)
}
if (nzchar(Sys.getenv("ANTHROPIC_API_KEY")) || nzchar(Sys.getenv("ANTHROPIC_AUTH_TOKEN"))) {
    backends$anthropic <- llmBackend("anthropic", model = Sys.getenv("ANTHROPIC_MODEL", "claude-haiku-4-5"))
}

score <- function(r, truth) {
    chosen <- r$changes[r$changes$rule == "llm", , drop = FALSE]
    got <- if (nrow(chosen)) stats::setNames(sub(" \\(.*", "", sub("^[^=]* = ", "", chosen$detail)),
        sub(" = .*", "", chosen$detail)) else character()
    if (all(is.na(truth))) {
        return(if (r$status == "failed") "correct (refused)" else "WRONG ACCEPTED")
    }
    if (startsWith(truth[[1]], "rules")) {
        alt <- any(r$changes$rule == "alternative")
        return(if (r$status == "repaired" && alt) "correct (rules)" else "not repaired")
    }
    if (r$status != "repaired") return("not repaired")
    if (all(names(truth) %in% names(got)) && all(got[names(truth)] == truth)) "correct" else "WRONG ACCEPTED"
}

rows <- list()
for (cs in names(cases)) {
    x <- cases[[cs]]$object
    truth <- cases[[cs]]$truth
    for (bn in names(backends)) {
        t0 <- Sys.time()
        r <- suppressWarnings(repairMapping("seurat-assay5", x, llm = backends[[bn]],
            request_file = file.path(tempdir(), sprintf("request_%s_%s.md", cs, bn))))
        lg <- r$llm_log
        chosen <- r$changes$detail[r$changes$rule == "llm"]
        rows[[length(rows) + 1L]] <- data.frame(case = cs, backend = bn,
            model = backends[[bn]]$model, status = r$status, outcome = score(r, truth),
            choices = if (length(chosen)) paste(sub(" \\(.*", "", chosen), collapse = "; ") else "-",
            input_tokens = lg$input_tokens %null% NA, output_tokens = lg$output_tokens %null% NA,
            llm_s = lg$latency_s %null% NA, total_s = round(as.numeric(Sys.time() - t0, units = "secs"), 1),
            returned = lg$returned %null% NA, kept = lg$kept %null% NA,
            error = substr(lg$error %null% "", 1, 60), stringsAsFactors = FALSE)
        cat(sprintf("%s %-9s %-9s %-18s %s\n", cs, bn, r$status, rows[[length(rows)]]$outcome,
            rows[[length(rows)]]$choices))
    }
}
tab <- do.call(rbind, rows)
dir.create(file.path("examples", "seurat_comparison", "results"), showWarnings = FALSE)
write.csv(tab, file.path("examples", "seurat_comparison", "results", "phase5_llm.csv"), row.names = FALSE)
cat("\n")
print(tab, row.names = FALSE, right = FALSE)
