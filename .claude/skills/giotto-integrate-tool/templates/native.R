# Template: native reimplementation (Giotto/R/<method>.R).
# Worked example: Giotto R/neighborhood-pca.R (BANKSY embedding as reduceData).
# Give it a descriptive name, not the original tool's.

#' @include pca-param.R
NULL

#' @title <Method> parameters
#' @name <method>Param
#' @description
#' Parameters for [reduceData()] that compute <what> as in <Tool>
#' (<Authors> <year>). This is a reimplementation of <Tool>'s <part>,
#' <differences, e.g. without X>. **If you use it, please cite the original
#' paper** (see References). [run<Tool>()] runs the original package.
#' @references <full reference with \doi{}>
#' @md
#' @export
<method>Param <- function(<options with the original's defaults>) {
    p <- new("<method>Param", param = list())
    p
}

#' @rdname <method>Param
#' @exportClass <method>Param
setClass("<method>Param", contains = "reduceParam")

#' @rdname <method>Param
#' @export
setMethod("reduceData", signature(x = "allMatrix", param = "<method>Param"),
    function(x, param, ...) {
        # compose existing blocks (createNetwork, reduceData(pcaParam), ...);
        # private helpers only for what does not exist yet
    }
)

#' @title <Method>
#' @name run<Method>
#' @description <what>; reimplements <Tool> (please cite it; see References).
#' @references <full reference with \doi{}>
#' @md
#' @export
run<Method> <- function(gobject, spat_unit = NULL, feat_type = NULL, ...,
    name = "<method>", verbose = NULL) {
    # defaults, inputs, res <- reduceData(x, param, ...), createDimObj + setDimReduction
    update_giotto_params(gobject, description = "_<method>",
        attachments = list(param = param))
}

# Parity (examples/<method>_comparison/): run the original on the same input;
# compare intermediate matrices (max abs diff), embeddings (|cor| of PCs) and
# labels (aricode::ARI). Report differences you find in the original.
