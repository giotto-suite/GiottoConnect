# Template: connection to an external tool (GiottoConnect/R/tool-<name>.R).
# Worked example: R/tool-ontrac.R (Python CLI). Replace <name>/<Tool>.

#' @title <Tool> parameters
#' @name <name>Param
#' @description
#' Parameters for [GiottoClass::analyzeData()] that run <Tool> (<Authors> <year>),
#' which <one line on what it does>. <Tool> is <an R package from ... /
#' a Python package installed on first use>.
#' **If you use it, please cite the original paper** (see References).
#' @param ... <tool options>
#' @returns A `<name>Param` object.
#' @references <Authors>. <Title>. *<Journal>* <vol>, <pages> (<year>).
#'   \doi{<doi>}
#' @seealso [run<Tool>()]
#' @md
#' @export
<name>Param <- function(<options with the tool's defaults>, ...) {
    # checkmate::assert_*() each option
    p <- methods::new("<name>Param", param = list(...))
    # p$<option> <- <option>
    p
}

#' @rdname <name>Param
#' @exportClass <name>Param
setClass("<name>Param", contains = "analyzeParam")

#' @rdname <name>Param
#' @param x <input data class and its required columns>
#' @param param a `<name>Param`
#' @export
setMethod("analyzeData", signature(x = "<data.frame|allMatrix|...>", param = "<name>Param"),
    function(x, param, ...) {
        # R tool:      package_check("<pkg>", repository = "<CRAN|Bioc|github:owner/repo>")
        #              fn <- getExportedValue("<pkg>", "<fun>")   # for GitHub-only tools
        # Python tool: exe <- .tool_python_cli("<pkg>==<ver>", "<module>", "<cli>")
        #              .tool_run(exe, args, log_file)
        # return a data.table with cell_ID (+ attributes for extra outputs)
    }
)

#' @title Run <Tool>
#' @name run<Tool>
#' @description Runs <Tool> on <inputs> and stores <outputs> in <slots>.
#' **If you use it, please cite the original paper** (see [<name>Param()]).
#' @param gobject giotto object
#' @param spat_unit,feat_type spatial unit and feature type
#' @param param a [<name>Param()]; its defaults when `NULL`
#' @param name prefix of the results
#' @returns The giotto object with the results (or, when the results fit no
#'   slot, the result with `attr(, "param")`).
#' @inherit <name>Param references
#' @md
#' @export
run<Tool> <- function(gobject, spat_unit = NULL, feat_type = NULL, param = NULL,
    name = "<name>", verbose = NULL) {
    param <- param %null% <name>Param()
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject, spat_unit = spat_unit,
        feat_type = feat_type)
    # x <- <inputs from accessors>
    res <- analyzeData(x, param)
    # write with the setters of the result -> slot table
    update_giotto_params(gobject, description = "_<name>",
        attachments = list(param = param))
}
