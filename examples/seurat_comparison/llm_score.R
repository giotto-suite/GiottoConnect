# outcome of a repair against a case's truth (shared by llm_repair.R and llm_replay.R)
llm_score <- function(r, truth) {
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

