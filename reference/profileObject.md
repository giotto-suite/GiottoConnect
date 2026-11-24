# Profile the structure of an object

Describe the structure of any R object (S4, lists, environments,
matrices, tables) as a flat node table, without loading the package that
defines it. Profiles are used to choose a mapping, to compute structural
fingerprints and to detect drift between object versions.

## Usage

``` r
profileObject(
  x,
  max_depth = 8L,
  max_children = 50L,
  sample_n = 10000L,
  n_cells = NULL,
  n_feat = NULL
)
```

## Arguments

- x:

  object to profile

- max_depth:

  maximum depth walked

- max_children:

  maximum number of named children walked per node; unnamed lists are
  represented by their first element only

- sample_n:

  maximum number of values inspected per matrix or vector for role hints

- n_cells, n_feat:

  optional known numbers of cells and features. When `NULL` they are
  inferred from the largest count-like matrix.

## Value

a `giottoProfile`: a list with `source`, `profiled_at`,
`profiler_version`, `n_cells`, `n_feat`, `size_from`, `notes` and
`nodes` (a data.frame with columns `path`, `kind`, `class`, `dtype`,
`shape`, `n_children`, `truncated`, `role_hints`, `sample`, `value`;
`value` holds scalar vectors and is not part of the fingerprint)

## Examples

``` r
p <- profileObject(list(counts = matrix(1:6, 2), meta = data.frame(a = 1:3)))
p
```
