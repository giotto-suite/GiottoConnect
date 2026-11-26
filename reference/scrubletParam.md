# Scrublet parameters

Parameters for
[`GiottoClass::analyzeData()`](https://giotto-suite.github.io/GiottoClass/reference/analyzeData.html)
that run Scrublet (Wolock, Lopez and Klein 2019), which detects doublets
by simulating them from pairs of observed cells and scoring each cell by
the share of simulated doublets among its nearest neighbors. It is meant
for single-cell data, not for spots that hold several cells. Scrublet is
a Python package: the first run installs it with reticulate (reticulate
\>= 1.41, no conda). **If you use it, please cite the original paper**
(see References).

## Usage

``` r
scrubletParam(
  expected_doublet_rate = 0.06,
  min_counts = 1L,
  min_cells = 1L,
  min_gene_variability_pctl = 85,
  n_prin_comps = 30L,
  threshold = NULL,
  seed = 1234L
)

# S4 method for class 'allMatrix,scrubletParam'
analyzeData(x, param, verbose = FALSE, ...)
```

## Arguments

- expected_doublet_rate:

  expected fraction of doublets (0.06 is from the 10x Chromium guide)

- min_counts, min_cells:

  Scrublet's feature filter before its PCA: a feature needs at least
  `min_counts` counts in at least `min_cells` cells

- min_gene_variability_pctl:

  percentile of feature variability above which features are used

- n_prin_comps:

  number of principal components

- threshold:

  `NULL` (Scrublet sets the doublet score threshold itself) or a score
  above which cells are called doublets. Scrublet cannot always set one,
  e.g. on spot data, and then calls no doublets.

- seed:

  integer random state of Scrublet

- x:

  count matrix, features x cells

- param:

  a `scrubletParam`

- verbose:

  show Scrublet's progress messages

## Value

A `scrubletParam` object.

## Details

`analyzeData(x, param)` returns a `data.table` with one row per cell:
`cell_ID`, `doublet_score` and `doublet` (Scrublet's call, `NA` when
there is no `threshold` and Scrublet could not set one).

## References

Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of
cell doublets in single-cell transcriptomic data. *Cell Systems* 8,
281-291.e9 (2019).
[doi:10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005)

## See also

[`runScrublet()`](https://giotto-suite.github.io/GiottoConnect/reference/runScrublet.md)

## Examples

``` r
p <- scrubletParam(expected_doublet_rate = 0.08)
p$expected_doublet_rate
```
