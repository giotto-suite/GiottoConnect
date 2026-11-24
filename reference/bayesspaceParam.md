# BayesSpace parameters

Parameters for
[`GiottoClass::clusterData()`](https://giotto-suite.github.io/GiottoClass/reference/clusterData.html)
that run BayesSpace (Zhao et al. 2021), which clusters spots of
array-based spatial transcriptomics (Visium, ST, Visium HD bins) with a
Bayesian model whose Potts prior favours the same cluster for
neighbouring spots. BayesSpace is a Bioconductor package
(`BiocManager::install("BayesSpace")`). **If you use it, please cite the
original paper** (see References).

## Usage

``` r
bayesspaceParam(
  q,
  platform = c("Visium", "ST", "VisiumHD"),
  d = 15L,
  n_pcs = 15L,
  n_hvgs = 2000L,
  log_normalize = TRUE,
  init_method = c("mclust", "kmeans"),
  model = c("t", "normal"),
  precision = c("equal", "variable"),
  nrep = 50000L,
  burn_in = 1000L,
  thin = 100L,
  gamma = NULL,
  mu0 = NULL,
  lambda0 = NULL,
  alpha = 1,
  beta = 0.01,
  seed = 1234L
)

# S4 method for class 'allMatrix,bayesspaceParam'
clusterData(x, param, spatial_locs, pcs = NULL, ...)
```

## Arguments

- q:

  integer, number of clusters (no default in BayesSpace; its `qTune()`
  helps choose one)

- platform:

  `"Visium"` (hexagonal grid, 6 neighbours), `"ST"` or `"VisiumHD"`
  (square grid, 4 neighbours)

- d:

  integer, number of principal components used

- n_pcs, n_hvgs:

  integer, components and highly variable genes of BayesSpace's
  preprocessing (`spatialPreprocess()`); not used when a stored
  dimension reduction is given

- log_normalize:

  logical, log-normalize counts before the PCA

- init_method:

  `"mclust"` or `"kmeans"`, initial clustering

- model:

  `"t"` or `"normal"`, error model

- precision:

  `"equal"` or `"variable"`, cluster covariance

- nrep, burn_in, thin:

  integer, MCMC iterations, burn-in and thinning

- gamma:

  numeric, smoothing strength; `NULL` uses BayesSpace's default (3 for
  Visium, 2 otherwise)

- mu0, lambda0:

  prior mean and precision of the cluster means; `NULL` uses
  BayesSpace's defaults

- alpha, beta:

  numeric, priors of the t-model scaling

- seed:

  integer seed set before preprocessing and clustering; `NULL` leaves
  the random state alone

- x:

  count matrix, features x cells (dense or sparse)

- spatial_locs:

  `data.frame` with `cell_ID` and either the array grid positions
  `array_row`, `array_col` or the coordinates `sdimx`, `sdimy` (the grid
  is then derived from the spot spacing)

- pcs:

  optional matrix, cells x dimensions (rownames are cell IDs) to cluster
  on instead of BayesSpace's own preprocessing

## Value

A `bayesspaceParam` object.

## Details

`clusterData(x, param, spatial_locs)` returns a `data.table` with
`cell_ID`, `cluster` (BayesSpace's `spatial.cluster`) and `cluster_init`
(the initial non-spatial clustering). Its attributes `pcs` (the
components clustered) and `hvf` (BayesSpace's highly variable features,
when it preprocessed) hold the intermediate results.

## References

Zhao E, Stone MR, Ren X, Guenthoer J, Smythe KS, Pulliam T, Williams SR,
Uytingco CR, Taylor SEB, Nghiem P, Bielas JH, Gottardo R. Spatial
transcriptomics at subspot resolution with BayesSpace. *Nature
Biotechnology* 39, 1375-1384 (2021).
[doi:10.1038/s41587-021-00935-2](https://doi.org/10.1038/s41587-021-00935-2)

## See also

[`runBayesSpace()`](https://giotto-suite.github.io/GiottoConnect/reference/runBayesSpace.md)

## Examples

``` r
p <- bayesspaceParam(q = 7, nrep = 10000)
p$q
```
