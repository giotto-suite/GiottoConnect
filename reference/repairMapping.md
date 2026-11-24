# Repair a mapping for a drifted object

Adapt a mapping to an object whose structure differs from the one the
mapping was verified on. Parameters are re-bound first; then rules
rewrite moved or renamed source paths, add a `cast` when only the
storage of a node changed, and switch a part to one of its declared
alternatives when its source is gone. Whatever is still broken gets up
to five ranked candidate paths. A repaired mapping is validated (and,
when `x` is the object, converted and read back) before it is accepted.

## Usage

``` r
repairMapping(
  mapping,
  x,
  llm = NULL,
  old_profile = NULL,
  file = NULL,
  request_file = NULL
)
```

## Arguments

- mapping:

  a `giottoMapping`, mapping file or bundled mapping id

- x:

  the drifted object, or its `giottoProfile`

- llm:

  optional
  [`llmBackend()`](https://giotto-suite.github.io/GiottoConnect/reference/llmBackend.md).
  Parts the rules cannot repair, and declared alternatives whose params
  do not bind by rule, are offered to it as closed choices (see
  [`llmPatch()`](https://giotto-suite.github.io/GiottoConnect/reference/llmPatch.md)).
  Its choices pass the same validation, conversion and readback as rule
  repairs.

- old_profile:

  profile the mapping was verified on. Defaults to the
  `<id>.profile.json` stored next to the mapping file.

- file:

  optional path to write the repaired mapping to (YAML)

- request_file:

  optional path for a
  [`writeMappingRequest()`](https://giotto-suite.github.io/GiottoConnect/reference/writeMappingRequest.md)
  handoff, written when the repair fails

## Value

a `giottoRepair`: list with `status` (`"not_needed"`, `"repaired"` or
`"failed"`), `mapping` (the possibly patched mapping), `changes`,
`unresolved`, `candidates`, `validation` and `readback`
