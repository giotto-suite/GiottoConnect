# LLM-assisted repair. The model only selects: for each unresolved item it
# picks one option from a closed list built by the repair rules (a candidate
# source path, or a value for an unbound param). The schema sent with each
# request enumerates exactly those options; replies are checked again here.
# Nothing in the conversion path calls this module.

.LLM_DEFAULT_MODEL <- c(anthropic = "claude-haiku-4-5")

# what each Giotto part expects, for the prompt
.PART_EXPECTS <- c(
    expression = "an expression matrix, features x cells; raw counts are integer-valued and non-negative",
    cell_metadata = "a table with one row per cell",
    feat_metadata = "a table with one row per feature",
    spatial_locs = "x and y coordinates of cell centroids",
    spatial_info = "cell boundary polygons (vertex coordinates per cell)",
    dim_reduction = "a cells x dimensions embedding",
    nn_network = "a cell-cell nearest-neighbour graph"
)

#' @title LLM backend for mapping repair
#' @name llmBackend
#' @description Configure the model [repairMapping()] may ask to choose among
#' repair options. Settings come from the arguments or the environment:
#' `GIOTTO_LLM_PROVIDER` (`"anthropic"`, `"ollama"` or `"baseline"`),
#' `GIOTTO_LLM_MODEL`, and for Anthropic `ANTHROPIC_API_KEY` (or
#' `ANTHROPIC_AUTH_TOKEN`), for Ollama `OLLAMA_HOST` (default
#' `http://localhost:11434`). Keys are read at call time and never stored or
#' logged.
#'
#' `"baseline"` uses no model: it takes the top-ranked candidate path when its
#' score leads the next one by at least `min_lead`, and makes no param choices.
#' It is the reference an LLM backend has to beat.
#' @param provider `"anthropic"`, `"ollama"` or `"baseline"`
#' @param model model id. Anthropic defaults to `"claude-haiku-4-5"`; Ollama
#' needs one (e.g. `"qwen3.8:27b"`).
#' @param host Ollama host
#' @param max_tokens maximum response tokens
#' @param min_lead baseline only: required score lead of the top candidate
#' @returns an `llmBackend` list
#' @export
llmBackend <- function(provider = Sys.getenv("GIOTTO_LLM_PROVIDER"),
    model = Sys.getenv("GIOTTO_LLM_MODEL"),
    host = Sys.getenv("OLLAMA_HOST", "http://localhost:11434"),
    max_tokens = 2048L,
    min_lead = 0.5) {
    provider <- match.arg(provider, c("anthropic", "ollama", "baseline"))
    if (!nzchar(model)) model <- unname(.LLM_DEFAULT_MODEL[provider]) %na% ""
    if (provider == "ollama" && !nzchar(model)) {
        stop(wrap_txt("Ollama needs a model: set GIOTTO_LLM_MODEL or `model`"), call. = FALSE)
    }
    if (provider != "baseline") package_check("httr2")
    structure(list(provider = provider, model = if (provider == "baseline") "top-candidate" else model,
        host = host, max_tokens = as.integer(max_tokens), min_lead = min_lead),
        class = "llmBackend")
}

#' @export
print.llmBackend <- function(x, ...) {
    cat(sprintf("<llmBackend> %s / %s\n", x$provider, x$model))
    invisible(x)
}

#' @title Ask an LLM backend for a repair patch
#' @name llmPatch
#' @description Send one repair request and return the patches the backend
#' chose, after checking each against the request's options. Patches naming an
#' option that was not offered, or with `"low"` confidence, are dropped.
#' @param request a repair request built by [repairMapping()] (list with
#' `items`, `prompt` and `schema`)
#' @param backend an [llmBackend()]
#' @returns list with `patches` (data.frame: `item`, `value`, `confidence`,
#' `rationale`) and `log` (one-row data.frame: provider, model, prompt
#' characters, input/output tokens, latency, patches returned and kept, error)
#' @export
llmPatch <- function(request, backend) {
    checkmate::assert_class(backend, "llmBackend")
    t0 <- Sys.time()
    res <- tryCatch(
        switch(backend$provider,
            anthropic = .llm_anthropic(request, backend),
            ollama = .llm_ollama(request, backend),
            baseline = .llm_baseline(request, backend)
        ),
        error = function(e) list(error = conditionMessage(e), patches = list())
    )
    raw <- res$patches %null% list()
    kept <- .llm_check(raw, request)
    log <- data.frame(provider = backend$provider, model = backend$model,
        prompt_chars = nchar(request$prompt), input_tokens = res$input_tokens %null% NA_integer_,
        output_tokens = res$output_tokens %null% NA_integer_,
        latency_s = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2),
        returned = length(raw), kept = nrow(kept), error = res$error %null% "",
        stringsAsFactors = FALSE)
    list(patches = kept, log = log)
}


# request ####

# One request covering every unresolved item. An item is either the source
# of a part (options: candidate paths) or an unbound param (options: names
# from the param's listing path).
.llm_request <- function(m, items, new, diff) {
    options <- lapply(items, `[[`, "options")
    names(options) <- vapply(items, `[[`, character(1), "item")
    any_of <- lapply(names(options), function(id) {
        list(type = "object",
            properties = list(
                item = list(type = "string", const = id),
                value = list(type = "string", enum = as.list(options[[id]])),
                confidence = list(type = "string", enum = list("high", "medium", "low")),
                rationale = list(type = "string")
            ),
            required = list("item", "value", "confidence", "rationale"),
            additionalProperties = FALSE)
    })
    schema <- list(type = "object",
        properties = list(patches = list(type = "array", items = list(anyOf = any_of))),
        required = list("patches"), additionalProperties = FALSE)
    list(items = items, schema = schema, prompt = .llm_prompt(m, items, new, diff))
}

.llm_prompt <- function(m, items, new, diff) {
    nd <- new$nodes
    describe <- function(path) {
        r <- nd[nd$path == path, , drop = FALSE]
        if (!nrow(r)) return(sprintf("- `%s`", path))
        sprintf("- `%s`: %s %s%s; roles: %s%s", path, r$kind, r$class,
            if (is.na(r$shape)) "" else paste0(" ", r$shape),
            r$role_hints %na% "none",
            if (r$kind %in% c("vector", "factor")) "" else
                sprintf("; names: %s", r$sample %na% "-"))
    }
    blocks <- vapply(items, function(it) {
        opts <- vapply(it$option_paths, describe, character(1))
        paste0(
            sprintf("### item `%s`\n", it$item),
            sprintf("%s\n", it$question),
            sprintf("Giotto expects %s.\n", .PART_EXPECTS[[it$part]]),
            "Options:\n", paste(opts, collapse = "\n"), "\n"
        )
    }, character(1))
    changes <- if (!is.null(diff) && nrow(diff)) {
        d <- as.data.frame(unclass(diff))
        paste(sprintf("- %s: %s%s", d$change, ifelse(is.na(d$path), "", d$path), ifelse(is.na(d$new_path), "",
            paste0(" -> ", d$new_path))), collapse = "\n")
    } else "- (no stored profile to compare with)"
    paste0(
        "A Seurat object is being converted to Giotto with a mapping, but its structure ",
        "differs from the one the mapping was verified on. Rules could not resolve the ",
        "items below. For each item, choose the option that holds the requested data. ",
        "Only choose when the option's path, name, kind and roles fit; give confidence ",
        "\"low\" when unsure, and leave an item out if no option fits. Each choice is ",
        "validated and the converted object is read back before anything is accepted.\n\n",
        sprintf("Mapping: %s v%s\n\nStructural changes since the mapping was verified:\n%s\n\n",
            m$mapping$id, m$mapping$version, changes),
        paste(blocks, collapse = "\n"),
        "\nAnswer with JSON matching the schema: one entry per item you resolve."
    )
}

# drop patches for options that were not offered, low confidence, duplicates
.llm_check <- function(raw, request) {
    empty <- data.frame(item = character(), value = character(), confidence = character(),
        rationale = character(), stringsAsFactors = FALSE)
    if (!length(raw)) return(empty)
    opts <- stats::setNames(lapply(request$items, `[[`, "options"),
        vapply(request$items, `[[`, character(1), "item"))
    rows <- lapply(raw, function(p) {
        if (!is.list(p) || is.null(p$item) || is.null(p$value)) return(NULL)
        if (!p$item %in% names(opts) || !p$value %in% opts[[p$item]]) return(NULL)
        if (!identical(p$confidence, "high") && !identical(p$confidence, "medium")) return(NULL)
        data.frame(item = p$item, value = p$value, confidence = p$confidence,
            rationale = substr(p$rationale %null% "", 1L, 300L), stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    if (is.null(out)) return(empty)
    out[!duplicated(out$item), , drop = FALSE]
}


# providers ####

.llm_anthropic <- function(request, backend) {
    key <- Sys.getenv("ANTHROPIC_API_KEY")
    token <- Sys.getenv("ANTHROPIC_AUTH_TOKEN")
    if (!nzchar(key) && !nzchar(token)) {
        stop("no Anthropic credentials: set ANTHROPIC_API_KEY (or ANTHROPIC_AUTH_TOKEN)")
    }
    body <- list(
        model = backend$model,
        max_tokens = backend$max_tokens,
        messages = list(list(role = "user", content = request$prompt)),
        output_config = list(format = list(type = "json_schema", schema = request$schema))
    )
    # sampling parameters are rejected by newer models
    if (grepl("haiku-4-5|-4-6|-4-5", backend$model)) body$temperature <- 0
    req <- httr2::request("https://api.anthropic.com/v1/messages")
    req <- httr2::req_headers(req, `anthropic-version` = "2023-06-01",
        `content-type` = "application/json")
    req <- if (nzchar(key)) {
        httr2::req_headers(req, `x-api-key` = key, .redact = "x-api-key")
    } else {
        httr2::req_headers(req, Authorization = paste("Bearer", token),
            `anthropic-beta` = "oauth-2025-04-20", .redact = "Authorization")
    }
    req <- httr2::req_body_json(req, body, auto_unbox = TRUE, null = "null")
    req <- httr2::req_error(req, is_error = function(resp) FALSE)
    req <- httr2::req_timeout(req, 120)
    resp <- httr2::req_perform(req)
    out <- httr2::resp_body_json(resp, simplifyVector = FALSE)
    if (httr2::resp_status(resp) >= 400) {
        stop(sprintf("Anthropic API %d: %s", httr2::resp_status(resp),
            out$error$message %null% "error"))
    }
    if (identical(out$stop_reason, "refusal") || identical(out$stop_reason, "max_tokens")) {
        stop(sprintf("Anthropic API stopped: %s", out$stop_reason))
    }
    txt <- paste(vapply(Filter(function(b) identical(b$type, "text"), out$content),
        `[[`, character(1), "text"), collapse = "")
    list(patches = jsonlite::fromJSON(txt, simplifyVector = FALSE)$patches,
        input_tokens = out$usage$input_tokens, output_tokens = out$usage$output_tokens)
}

.llm_ollama <- function(request, backend) {
    body <- list(
        model = backend$model,
        messages = list(list(role = "user", content = request$prompt)),
        format = request$schema,
        stream = FALSE,
        think = FALSE,
        options = list(temperature = 0, num_predict = backend$max_tokens)
    )
    req <- httr2::request(paste0(sub("/$", "", backend$host), "/api/chat"))
    req <- httr2::req_body_json(req, body, auto_unbox = TRUE, null = "null")
    req <- httr2::req_error(req, is_error = function(resp) FALSE)
    req <- httr2::req_timeout(req, 600)
    resp <- httr2::req_perform(req)
    out <- httr2::resp_body_json(resp, simplifyVector = FALSE)
    if (httr2::resp_status(resp) >= 400) {
        stop(sprintf("Ollama %d: %s", httr2::resp_status(resp), out$error %null% "error"))
    }
    list(patches = jsonlite::fromJSON(out$message$content, simplifyVector = FALSE)$patches,
        input_tokens = out$prompt_eval_count, output_tokens = out$eval_count)
}

# no model: the top-ranked candidate path when it clearly leads
.llm_baseline <- function(request, backend) {
    patches <- list()
    for (it in request$items) {
        if (it$kind != "from" || !length(it$scores)) next
        s <- sort(it$scores, decreasing = TRUE)
        lead <- if (length(s) > 1L) s[[1]] - s[[2]] else Inf
        if (lead >= backend$min_lead) {
            patches[[length(patches) + 1L]] <- list(item = it$item, value = names(s)[1],
                confidence = "high", rationale = sprintf("top candidate, score lead %.2f", lead))
        }
    }
    list(patches = patches)
}
