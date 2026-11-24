# The Giotto contract

Machine-readable description of the Giotto parts a mapping can target:
for each part the subobject class, constructor (and its arguments),
setter, nesting keys, matrix orientation and whether the part is written
to disk when the giotto object has a GiottoDisk source.

## Usage

``` r
giottoContract(source = c("auto", "bundled", "live"))
```

## Arguments

- source:

  `"auto"` (bundled contract for the installed GiottoClass version,
  otherwise generated live), `"bundled"` or `"live"`

## Value

a `giottoContract` list

## Examples

``` r
ct <- giottoContract("bundled")
names(ct$parts)
```
