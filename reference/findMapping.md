# Find the mapping for an object

Choose the mapping that applies to a profiled object. Mappings are
filtered by their `applies_to` declaration (class, package version and
structural predicates after parameter binding); among those that apply,
a fingerprint match is preferred, then the more specific mapping (more
`applies_to` predicates satisfied), then the one whose source paths
(optional parts included) are best covered; without an exact match the
status is `"repair_needed"`. Fingerprints cover the mapping's required
parts: an optional part that is absent from the object does not change
it.

## Usage

``` r
findMapping(profile, mappings = NULL, ecosystem = NULL, target = NULL)
```

## Arguments

- profile:

  `giottoProfile` of the object

- mappings:

  list of mappings (`giottoMapping`, files or ids). Defaults to all
  bundled mappings.

- ecosystem:

  optional ecosystem filter, e.g. `"seurat"`

- target:

  optional target filter (`mapping.target.contract`), e.g. `"seurat"`
  for mappings from giotto objects to Seurat

## Value

list with `status` (`"match"`, `"repair_needed"` or `"none"`),
`mapping`, `params`, `fingerprint`, `candidates` (one row per mapping
considered) and `reasons`
