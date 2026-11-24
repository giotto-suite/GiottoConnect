# Convert an object with a mapping

Convert an external object into a giotto object using only a mapping
file and the engine: parameters are bound on the object's profile, each
source node is read through its declared accessor, ops or a hook are
applied and the result is built with the GiottoClass constructor for the
target part. A part whose source is absent switches to the first of its
declared alternatives that is present; otherwise an optional part is
skipped. With `source`, the giotto object is backed by a GiottoDisk
project and disk-backable parts (expression, polygons, networks) are
written there.

## Usage

``` r
convertWithMapping(
  x,
  mapping,
  profile = NULL,
  source = NULL,
  validate = TRUE,
  verbose = FALSE,
  params = NULL,
  output = NULL
)
```

## Arguments

- x:

  object to convert

- mapping:

  a `giottoMapping`, mapping file or bundled mapping id

- profile:

  optional `giottoProfile` of `x`

- source:

  `NULL` for an in-memory giotto object, or a directory for a GiottoDisk
  project (requires GiottoDisk)

- validate:

  logical. Run
  [`validateMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/validateMapping.md)
  before (structural) and after (semantic) conversion; FAIL stops, WARN
  warns.

- verbose:

  logical

- params:

  optional named list of param values that replace the mapping's rules,
  e.g. `list(spat_unit = "nucleus", feat_type = "protein")` to export a
  chosen spatial unit and feature type; recorded in the provenance and
  reused by
  [`readbackCheck()`](https://giotto-suite.github.io/GiottoConnect/reference/readbackCheck.md)

- output:

  file to write, for file targets (e.g. `.h5ad`)

## Value

giotto object. The conversion (mapping id, version, fingerprint, params,
losses, alternatives used and absent optional parts) is recorded in the
object history.
