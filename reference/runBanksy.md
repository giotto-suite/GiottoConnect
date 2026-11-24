# Spatial domains with BANKSY

Runs BANKSY (Singhal et al. 2024) through the Bioconductor package
Banksy (see
[`banksyParam()`](https://giotto-suite.github.io/GiottoConnect/reference/banksyParam.md)),
stores its PCs as a dimension reduction (method `"banksy"`) and clusters
them into domains (cell metadata `name`): by default with Giotto's sNN
network and Leiden, as
`Giotto::findSpatialDomains(method = "neighborhood")` does with the
native embedding, so the two can be compared directly.

The input is expression by default, which is what Banksy is designed
for. With `input = "dim_reduction"` (e.g. PCs) Banksy discards the sign
of the neighbour means (see
[`banksyParam()`](https://giotto-suite.github.io/GiottoConnect/reference/banksyParam.md))
and a warning is given;
[`Giotto::runNeighborhoodPCA()`](https://giottosuite.com/reference/runNeighborhoodPCA.html)
keeps it.

**If you use it, please cite the original paper** (see References).

## Usage

``` r
runBanksy(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  input = c("expression", "dim_reduction"),
  expression_values = "normalized",
  feats_to_use = "hvf",
  dim_reduction_to_use = "pca",
  dim_reduction_name = NULL,
  dimensions_to_use = NULL,
  spat_loc_name = "raw",
  lambda = 0.8,
  k_geom = 15L,
  npcs = 20L,
  compute_agf = FALSE,
  seed = 1234L,
  cluster_method = c("giotto", "banksy", "none"),
  k_neighbors = 30L,
  resolution = 0.8,
  name = "banksy",
  verbose = NULL
)
```

## Arguments

- gobject:

  giotto object

- spat_unit, feat_type:

  spatial unit and feature type

- input:

  `"expression"` (default) or `"dim_reduction"`

- expression_values:

  expression values with `input = "expression"`

- feats_to_use:

  features with `input = "expression"`: a character vector or the name
  of a logical feature metadata column (default `"hvf"`); `NULL` uses
  all

- dim_reduction_to_use, dim_reduction_name, dimensions_to_use:

  the embedding with `input = "dim_reduction"`

- spat_loc_name:

  name of the spatial locations

- lambda, k_geom, npcs, compute_agf, seed:

  see
  [`banksyParam()`](https://giotto-suite.github.io/GiottoConnect/reference/banksyParam.md)

- cluster_method:

  `"giotto"` (sNN network + Leiden), `"banksy"`
  (`Banksy::clusterBanksy()`) or `"none"`

- k_neighbors:

  neighbours of the sNN network

- resolution:

  Leiden resolution

- name:

  name of the PCs, the network and the domain column

- verbose:

  be verbose

## Value

The giotto object with the BANKSY PCs and the domains.

## References

Singhal V, Chou N, Lee J, et al. BANKSY unifies cell typing and tissue
domain segmentation for scalable spatial omics data analysis. *Nature
Genetics* 56, 431–441 (2024).
[doi:10.1038/s41588-024-01664-3](https://doi.org/10.1038/s41588-024-01664-3)

## See also

[`banksyParam()`](https://giotto-suite.github.io/GiottoConnect/reference/banksyParam.md),
[`Giotto::runNeighborhoodPCA()`](https://giottosuite.com/reference/runNeighborhoodPCA.html),
[`Giotto::findSpatialDomains()`](https://giottosuite.com/reference/findSpatialDomains.html)

## Examples

``` r
if (FALSE) { # \dontrun{
g <- GiottoData::loadGiottoMini("visium")
g <- runBanksy(g, lambda = 0.8)
} # }
```
