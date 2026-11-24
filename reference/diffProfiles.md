# Compare two profiles

Classify structural changes between two profiles. Removed and added
subtrees are reported once, at their root. A removed subtree that
reappears elsewhere with the same kind, dtype, dimensionality and size
roles is reported as `renamed_key` (same parent) or `moved` (different
parent).

## Usage

``` r
diffProfiles(old, new, params_old = NULL, params_new = NULL)
```

## Arguments

- old, new:

  `giottoProfile` objects

- params_old, params_new:

  optional bound object names (e.g. assay names), templated as in
  [`fingerprint()`](https://giotto-suite.github.io/GiottoConnect/reference/fingerprint.md)
  before comparing

## Value

a `giottoProfileDiff` data.frame with columns `change` (`added`,
`removed`, `moved`, `renamed_key`, `kind_changed`, `dtype_changed`),
`path`, `new_path`, `kind`, `class`, `detail`

## Examples

``` r
a <- profileObject(list(x = list(counts = matrix(1:4, 2))))
b <- profileObject(list(x = list(raw = matrix(1:4, 2))))
diffProfiles(a, b)
```
