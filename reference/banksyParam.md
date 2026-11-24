# BANKSY parameters

Parameters for
[`GiottoClass::reduceData()`](https://giotto-suite.github.io/GiottoClass/reference/reduceData.html)
that run BANKSY (Singhal et al. 2024) through the Bioconductor package
Banksy: each cell's values are joined by the mean (and, with
`compute_agf = TRUE`, the azimuthal gradient) of its spatial neighbours,
weighted by `lambda`, and reduced by PCA. Banksy is suggested and
checked when this runs.

Banksy takes the absolute value of the neighbour means, which is
harmless for expression (never negative) but discards the sign of signed
inputs such as PCs. Use expression as the input, or the native
[`Giotto::runNeighborhoodPCA()`](https://giottosuite.com/reference/runNeighborhoodPCA.html),
which keeps the sign.

**If you use it, please cite the original paper** (see References).

## Usage

``` r
banksyParam(
  lambda = 0.8,
  k_geom = 15L,
  npcs = 20L,
  compute_agf = FALSE,
  seed = 1234L
)

# S4 method for class 'ANY,banksyParam'
reduceData(x, param, ..., spatial_locs)
```

## Arguments

- lambda:

  numeric in `[0, 1]`, weight of the neighbourhood (about 0.8 for
  spatial domains, about 0.2 for cell typing)

- k_geom:

  integer (or one per harmonic with `compute_agf = TRUE`), spatial
  neighbours

- npcs:

  integer, principal components

- compute_agf:

  logical, also use the azimuthal Gabor filter features

- seed:

  integer, seed of the PCA

- x:

  features x cells values with cell IDs as column names

- param:

  a `banksyParam`

- ...:

  unused

- spatial_locs:

  cells x 2 coordinates with cell IDs as row names

## Value

A `banksyParam` object.

## Details

`reduceData(x, param, spatial_locs = )` returns a list with `u` (cells x
`npcs` BANKSY PCs) and `se`, the `SpatialExperiment` Banksy worked on
(for `Banksy::clusterBanksy()`).

## References

Singhal V, Chou N, Lee J, et al. BANKSY unifies cell typing and tissue
domain segmentation for scalable spatial omics data analysis. *Nature
Genetics* 56, 431–441 (2024).
[doi:10.1038/s41588-024-01664-3](https://doi.org/10.1038/s41588-024-01664-3)

## See also

[`runBanksy()`](https://giotto-suite.github.io/GiottoConnect/reference/runBanksy.md)

## Examples

``` r
p <- banksyParam(lambda = 0.8)
```
