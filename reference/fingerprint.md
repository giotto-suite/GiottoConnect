# Structural fingerprint of a profile

Hash the structure of a profiled object: the sorted set of (templated
path, kind, dtype, ndim). Sizes and values are excluded, so two datasets
with the same structure share a fingerprint, and integer and double
storage count as the same type. Parameter values in paths are replaced
by their `{param}` placeholder first.

## Usage

``` r
fingerprint(
  profile,
  paths = NULL,
  params = NULL,
  scope = if (is.null(paths)) "full" else "paths"
)
```

## Arguments

- profile:

  a `giottoProfile`

- paths:

  templated paths to restrict the fingerprint to (plus their ancestors).
  Missing paths are part of the hash. Ignored when `scope = "full"`.

- params:

  named list of bound parameter values, e.g. `list(assay = "Xenium")`.
  For `scope = "paths"` they fill the `{param}` placeholders of `paths`.
  For `scope = "full"` their values are replaced by placeholders in the
  profile, as whole path segments or as `_`/`.`-delimited tokens
  (`Xenium_nn`, `FindNeighbors.Xenium.pca`); pass only object names
  there, since a layer name like `data` would also match `meta.data`.

- scope:

  `"paths"` (default when `paths` is given) or `"full"`

## Value

character, an xxhash64 digest

## Examples

``` r
p <- profileObject(list(a = list(counts = matrix(1:4, 2))))
fingerprint(p, scope = "full")
fingerprint(p, paths = "a/counts")
```
