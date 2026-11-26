# Run MAST differential expression

Runs MAST (Finak et al. 2015) to find marker features between clusters
of a giotto object and returns its per-feature tests (see
[`mastParam()`](https://giotto-suite.github.io/GiottoConnect/reference/mastParam.md)
for the columns). The giotto object is not changed.

- With `group_1` and `group_2`, one comparison: the cells of the
  `group_1` clusters against those of the `group_2` clusters. All
  features are returned.

- Without them, each cluster against all other clusters (those in
  `subset_clusters`, or all). Per cluster, the `min_feats` best features
  are kept and any other with `fdr < pval` and `coef > logFC`; the
  `ranking` column gives the rank within the cluster.

MAST fits a model per feature on a dense matrix, so it is slow on large
data. **If you use it, please cite the original paper** (see
[`mastParam()`](https://giotto-suite.github.io/GiottoConnect/reference/mastParam.md)).

## Usage

``` r
runMAST(
  gobject,
  cluster_column,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = "normalized",
  group_1 = NULL,
  group_2 = NULL,
  group_1_name = NULL,
  group_2_name = NULL,
  subset_clusters = NULL,
  adjust_columns = NULL,
  pval = 0.001,
  logFC = 1,
  min_feats = 10L,
  param = NULL,
  verbose = NULL
)
```

## Arguments

- gobject:

  giotto object

- cluster_column:

  cell metadata column with the clusters

- spat_unit, feat_type:

  spatial unit and feature type

- expression_values:

  expression values to use (MAST expects log-scale values)

- group_1, group_2:

  clusters of the two sides of a single comparison

- group_1_name, group_2_name:

  names of the two sides (default: the cluster IDs joined by `_`)

- subset_clusters:

  clusters compared one against the others (default: all)

- adjust_columns:

  cell metadata columns to adjust for (e.g. the detection rate)

- pval, logFC, min_feats:

  filters of the one-against-the-others results

- param:

  a
  [`mastParam()`](https://giotto-suite.github.io/GiottoConnect/reference/mastParam.md);
  its defaults when `NULL`

- verbose:

  be verbose

## Value

A `data.table` of MAST's results.

## References

Finak G, McDavid A, Yajima M, Deng J, Gersuk V, Shalek AK, Slichter CK,
Miller HW, McElrath MJ, Prlic M, Linsley PS, Gottardo R. MAST: a
flexible statistical framework for assessing transcriptional changes and
characterizing heterogeneity in single-cell RNA sequencing data. *Genome
Biology* 16, 278 (2015).
[doi:10.1186/s13059-015-0844-5](https://doi.org/10.1186/s13059-015-0844-5)

## See also

[`mastParam()`](https://giotto-suite.github.io/GiottoConnect/reference/mastParam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
g <- GiottoData::loadGiottoMini("visium")
# cluster 1 against clusters 2 and 3
runMAST(g, cluster_column = "leiden_clus", group_1 = 1, group_2 = c(2, 3))
# every cluster against the others
runMAST(g, cluster_column = "leiden_clus")
} # }
```
