# SPARK-X parameters

Parameters for
[`GiottoClass::analyzeData()`](https://giotto-suite.github.io/GiottoClass/reference/analyzeData.html)
that run SPARK-X (Zhu, Sun and Zhou 2021), which finds spatially
variable genes with a non-parametric covariance test between each gene's
counts and spatial kernels built from the locations. It scales to
hundreds of thousands of cells. SPARK-X is part of the SPARK package on
GitHub, install it with `remotes::install_github("xzhoulab/SPARK")`.
**If you use it, please cite the original paper** (see References).

## Usage

``` r
sparkxParam(option = c("mixture", "single"), n_cores = 1L)

# S4 method for class 'allMatrix,sparkxParam'
analyzeData(x, param, spatial_locs, X = NULL, verbose = FALSE, ...)
```

## Arguments

- option:

  `"mixture"` (SPARK-X's default: the projection kernel, five Gaussian
  and five cosine kernels, combined with the Cauchy combination rule) or
  `"single"` (projection kernel only)

- n_cores:

  integer, cores SPARK-X uses

- x:

  count matrix, features x cells

- param:

  a `sparkxParam`

- spatial_locs:

  two-column matrix of coordinates, rows in the column order of `x`

- X:

  `NULL` or a covariate model matrix (cells x covariates, with an
  intercept column), e.g. `model.matrix(~ cell_type)`

- verbose:

  show SPARK-X's progress messages

## Value

A `sparkxParam` object.

## Details

`analyzeData(x, param, spatial_locs, X)` returns a `data.table` with one
row per tested feature (`feat_ID`): `pval`, the combined p-value over
the kernels, and `padj`, adjusted with Benjamini-Yekutieli, as SPARK-X
reports them. SPARK-X leaves out cells and features without counts.

## References

Zhu J, Sun S, Zhou X. SPARK-X: non-parametric modeling enables scalable
and robust detection of spatial expression patterns for large spatial
transcriptomic studies. *Genome Biology* 22, 184 (2021).
[doi:10.1186/s13059-021-02404-0](https://doi.org/10.1186/s13059-021-02404-0)

## See also

[`runSparkX()`](https://giotto-suite.github.io/GiottoConnect/reference/runSparkX.md)

## Examples

``` r
p <- sparkxParam(option = "single")
p$option
```
