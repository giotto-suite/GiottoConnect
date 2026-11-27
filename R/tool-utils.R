# Shared helpers for the external tool connections (R/tool-<name>.R).
# Tools are never imported: R tools are checked with package_check() when
# their function runs, Python tools run in a worker R session with their own
# environment through GiottoUtils::py_run() (adr/0008 in GiottoClass).

# A journal article for methodCitation(). Keep it in sync with @references.
.cite_article <- function(author, title, journal, year, volume, pages, doi) {
    utils::bibentry("Article",
        author = author, title = title, journal = journal, year = year,
        volume = volume, pages = pages, doi = doi
    )
}

# Features a tool tests, from `feats_to_use`: NULL (all), a feature metadata
# column (features with "yes" or TRUE) or feature IDs.
.tool_feats <- function(feats_to_use, fm, all_feats, fn) {
    if (is.null(feats_to_use)) return(all_feats)
    if (length(feats_to_use) == 1L && feats_to_use %in% names(fm)) {
        v <- fm[[feats_to_use]]
        sel <- fm$feat_ID[v %in% c("yes", TRUE)]
    } else {
        sel <- feats_to_use
    }
    sel <- intersect(all_feats, sel)
    if (!length(sel)) {
        stop("[", fn, "] `feats_to_use` selects no features", call. = FALSE)
    }
    sel
}
