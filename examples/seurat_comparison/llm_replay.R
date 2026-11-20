# Replay recorded LLM answers (e.g. from a Claude Code Haiku subagent given the
# exact repair prompt) through the same check and acceptance gate, and append
# the outcome to results/phase5_llm.csv. Answers: <dir>/<case>_answer.json.
# Usage: Rscript examples/seurat_comparison/llm_replay.R <answer dir> <label>
args <- commandArgs(trailingOnly = TRUE)
dir <- args[[1]]
label <- if (length(args) > 1) args[[2]] else "replay"
options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE, width = 160)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
env <- new.env()
sys.source(file.path("inst", "scripts", "build_fixtures.R"), envir = env)
cases <- suppressWarnings(env$build_llm_cases())
source(file.path("examples", "seurat_comparison", "llm_score.R"))

rows <- list()
for (cs in sub("_answer\\.json$", "", list.files(dir, pattern = "_answer\\.json$"))) {
    ans <- jsonlite::read_json(file.path(dir, paste0(cs, "_answer.json")), simplifyVector = FALSE)
    backend <- structure(list(provider = "ollama", model = label, host = "", max_tokens = 1L,
        reply = ans$patches), class = "llmBackend")
    r <- testthat::with_mocked_bindings(
        suppressWarnings(repairMapping("seurat-assay5", cases[[cs]]$object, llm = backend)),
        .llm_ollama = function(request, backend) list(patches = backend$reply),
        .package = "GiottoConnect")
    chosen <- r$changes$detail[r$changes$rule == "llm"]
    rows[[cs]] <- data.frame(case = cs, backend = "replay", model = label, status = r$status,
        outcome = llm_score(r, cases[[cs]]$truth),
        choices = if (length(chosen)) paste(sub(" \\(.*", "", chosen), collapse = "; ") else
            paste("(rejected)", paste(vapply(ans$patches, function(p) paste(p$item, "=", p$value,
                sprintf("[%s]", p$confidence)), ""), collapse = "; ")),
        input_tokens = NA, output_tokens = NA, llm_s = NA, total_s = NA,
        returned = length(ans$patches), kept = r$llm_log$kept, error = "", stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows)
out <- file.path("examples", "seurat_comparison", "results", "phase5_llm.csv")
old <- read.csv(out, stringsAsFactors = FALSE)
old <- old[!(old$backend == "replay" & old$model == label), , drop = FALSE]
write.csv(rbind(old, tab), out, row.names = FALSE)
print(tab[, c("case", "model", "status", "outcome", "choices", "kept")], row.names = FALSE, right = FALSE)
