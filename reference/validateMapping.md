# Validate a mapping against an object

Check a mapping against a profiled object (structural checks) and, when
`result` is given, check the converted giotto object (semantic checks).
No conversion is run here.

Structural: schema, parameter binding, source paths, target keys under
the Giotto contract, hooks, and data in the source that is neither
mapped nor declared lossy. Semantic: raw expression is integer-valued
and non-negative, cell IDs agree across parts, dimensions match the
profile, coordinate ranges are plausible for microns, reductions and
networks index the cells.

## Usage

``` r
validateMapping(
  mapping,
  profile,
  params = NULL,
  result = NULL,
  checks = c("all", "structural", "semantic")
)
```

## Arguments

- mapping:

  a `giottoMapping`, mapping file or bundled mapping id

- profile:

  `giottoProfile` of the source object

- params:

  optional bound params (as returned internally by the engine)

- result:

  optional giotto object produced with this mapping

- checks:

  `"all"`, `"structural"` or `"semantic"` (the latter requires `result`)

## Value

a `giottoValidation` data.frame with columns `check`, `part_id`, `level`
(OK, WARN, FAIL) and `message`
