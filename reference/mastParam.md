# MAST parameters

Parameters for
[`GiottoClass::analyzeData()`](https://giotto-suite.github.io/GiottoClass/reference/analyzeData.html)
that run MAST (Finak et al. 2015), which tests each feature for
differential expression between two groups of cells with a hurdle model:
a logistic regression for whether the feature is detected and a linear
regression for its level where it is. The likelihood ratio test combines
both parts. MAST is a Bioconductor package, install it with
`BiocManager::install("MAST")`. **If you use it, please cite the
original paper** (see References).

## Usage

``` r
mastParam(
  method = c("bayesglm", "glm", "ridge", "blmer", "glmer"),
  ebayes = TRUE
)

# S4 method for class 'allMatrix,mastParam'
analyzeData(
  x,
  param,
  groups,
  test,
  ref,
  covariates = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- method:

  fitting method passed to `MAST::zlm()`: `"bayesglm"` (MAST's default),
  `"glm"`, `"ridge"`, `"blmer"` or `"glmer"`

- ebayes:

  logical, shrink the variances with empirical Bayes (MAST's default
  `TRUE`)

- x:

  expression matrix, features x cells (MAST expects log-scale values; it
  is made dense)

- param:

  a `mastParam`

- groups:

  group label of each cell. Named by cell ID (matched to `colnames(x)`)
  or in the column order of `x`. Cells labelled `NA` are left out.

- test, ref:

  the labels in `groups` to compare: `test` against `ref`

- covariates:

  `NULL` or a data.frame with one row per column of `x` and the
  covariates to adjust for (e.g. the detection rate)

- verbose:

  show MAST's progress messages

## Value

A `mastParam` object.

## Details

`analyzeData(x, param, groups, test, ref, covariates)` returns a
`data.table` with one row per feature, sorted by `fdr`: `feats`,
`Pr(>Chisq)` (the hurdle likelihood ratio test), `coef`, `ci.hi` and
`ci.lo` (the log fold change of `test` over `ref` and its confidence
interval), `fdr` (Benjamini-Hochberg) and `cluster`
(`"<test>_vs_<ref>"`).

## References

Finak G, McDavid A, Yajima M, Deng J, Gersuk V, Shalek AK, Slichter CK,
Miller HW, McElrath MJ, Prlic M, Linsley PS, Gottardo R. MAST: a
flexible statistical framework for assessing transcriptional changes and
characterizing heterogeneity in single-cell RNA sequencing data. *Genome
Biology* 16, 278 (2015).
[doi:10.1186/s13059-015-0844-5](https://doi.org/10.1186/s13059-015-0844-5)

## See also

[`runMAST()`](https://giotto-suite.github.io/GiottoConnect/reference/runMAST.md)

## Examples

``` r
p <- mastParam(method = "glm")
p$method
```
