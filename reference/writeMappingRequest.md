# Write a mapping repair request for a human or a stronger model

Write a markdown handoff for a repair that failed: what the mapping
expected, what changed in the object, which parts are unresolved with
their ranked candidates, the validation and readback results, and the
contract lines for the affected parts.

## Usage

``` r
writeMappingRequest(repair, file)
```

## Arguments

- repair:

  a `giottoRepair` from
  [`repairMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/repairMapping.md)

- file:

  output path (`.md`)

## Value

`file`, invisibly
