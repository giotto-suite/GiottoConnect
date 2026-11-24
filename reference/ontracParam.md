# ONTraC parameters

Parameters for
[`GiottoClass::analyzeData()`](https://giotto-suite.github.io/GiottoClass/reference/analyzeData.html)
that run ONTraC (Wang et al. 2025), which builds niche trajectories with
a graph neural network on cell-type composition of each cell's
neighbourhood. ONTraC is a Python package; it is installed in a cached
environment the first time it runs (through
[`reticulate::py_require()`](https://rstudio.github.io/reticulate/reference/py_require.html),
needs reticulate \>= 1.41). **If you use it, please cite the original
paper** (see References).

## Usage

``` r
ontracParam(
  k = 6L,
  n_neighbors = 50L,
  n_local = 20L,
  epochs = 1000L,
  patience = 100L,
  min_epochs = 50L,
  lr = 0.03,
  hidden_feats = 4L,
  n_gcn_layers = 2L,
  modularity_loss_weight = 1,
  purity_loss_weight = 30,
  regularization_loss_weight = 0.1,
  beta = 0.3,
  seed = 42L,
  device = c("cpu", "cuda"),
  n_cpu = 4L,
  version = "2.1.2",
  out_dir = NULL,
  ...
)

# S4 method for class 'data.frame,ontracParam'
analyzeData(x, param, ...)
```

## Arguments

- k:

  integer, number of niche clusters (at most 8)

- n_neighbors:

  integer, neighbours for the niche network

- n_local:

  integer, index of the local neighbour used to normalize distances

- epochs, patience, min_epochs:

  integer, training length

- lr:

  numeric, learning rate

- hidden_feats, n_gcn_layers:

  integer, model size

- modularity_loss_weight, purity_loss_weight,
  regularization_loss_weight, beta:

  numeric, loss weights (ONTraC's defaults)

- seed:

  integer, seed

- device:

  `"cpu"` or `"cuda"`

- n_cpu:

  integer, CPUs for building the niche network

- version:

  ONTraC version installed on first use

- out_dir:

  directory for ONTraC's files; `NULL` uses a temporary one

- ...:

  further ONTraC command line options as `name = value`, e.g.
  `"trajectory-construct" = "BF"`

- x:

  `data.frame` with columns `cell_ID`, `sdimx`, `sdimy`, `cell_type` and
  optionally `sample`

- param:

  an `ontracParam`

## Value

An `ontracParam` object.

## Details

`analyzeData(x, param)` returns a `data.table` with `cell_ID`,
`NTScore`, `niche_cluster` (the most likely cluster) and one probability
column per niche cluster (`NicheCluster_0`, ...). Its attributes
`niche_clusters` (each cluster's NT score), `connectivity` (cluster x
cluster) and `out_dir` hold the niche-level results.

## References

Wang W, Zheng S, Shin SC, Chávez-Fuentes JC, Yuan GC. ONTraC
characterizes spatially continuous variations of tissue microenvironment
through niche trajectory analysis. *Genome Biology* 26, 117 (2025).
[doi:10.1186/s13059-025-03588-5](https://doi.org/10.1186/s13059-025-03588-5)

## See also

[`runONTraC()`](https://giotto-suite.github.io/GiottoConnect/reference/runONTraC.md)

## Examples

``` r
p <- ontracParam(k = 6, epochs = 200)
p$k
```
