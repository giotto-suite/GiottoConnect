# Read a mapping file

Read a YAML mapping file (or a mapping already parsed into a list) and
validate it against the mapping schema. All schema problems are reported
together, each with the YAML path it concerns.

## Usage

``` r
readMapping(x)
```

## Arguments

- x:

  path to a YAML mapping file, the id of a mapping bundled with
  GiottoConnect (e.g. `"seurat-assay5"`), or a parsed mapping list

## Value

a `giottoMapping` (the parsed mapping as a list)

## Examples

``` r
f <- system.file("mappings", package = "GiottoConnect")
list.files(f, pattern = "\\.yaml$")
```
