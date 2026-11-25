# Check a method integration

Compares a giotto object before and after running a method's `run*()`
wrapper and checks the rules of GiottoClass ADR 0008:

- only the declared slots changed (via
  [`GiottoClass::manifestDiff()`](https://giotto-suite.github.io/GiottoClass/reference/manifestDiff.html));

- something was stored;

- stored dimension reductions and spatial enrichments cover every cell
  of their spatial unit, and new cell metadata columns are not all `NA`;

- the last history entry carries the method's param
  (`update_giotto_params(attachments = list(param = param))`);

- the param declares its paper
  ([`GiottoClass::methodCitation()`](https://giotto-suite.github.io/GiottoClass/reference/methodCitation.html)),
  which is then shown once per session and listed by
  [`GiottoClass::objCitations()`](https://giotto-suite.github.io/GiottoClass/reference/objCitations.html);

- the help page of `fn` has references and asks users to cite them.

Methods whose results fit no slot (e.g. an interaction table) return
them instead: pass `slots = "none"` and the returned object as `result`.
Then the object must be unchanged and `attr(result, "param")` must hold
the param.

## Usage

``` r
checkGiottoMethod(
  before,
  after,
  slots,
  param_class = NULL,
  fn = NULL,
  package = NULL,
  result = NULL,
  cite = TRUE
)

expect_giotto_method(...)
```

## Arguments

- before, after:

  giotto objects before and after the run

- slots:

  slots the method may change, e.g.
  `c("cell_metadata", "spatial_enrichment")`, or `"none"` for methods
  that return their results

- param_class:

  class of the param expected in the history; `NULL` accepts any param

- fn, package:

  name of the `run*()` function and its package, for the documentation
  check; skipped when `NULL`

- result:

  the returned result, with \`slots = "none"

- cite:

  logical. Whether the param must declare a paper. `FALSE` only for a
  Suite-native method that wraps or reimplements no published method

- ...:

  arguments of `checkGiottoMethod()`

## Value

A `giottoMethodCheck`: a `data.table` with columns `check`, `pass` and
`message`, printed as a report.

## Details

`expect_giotto_method()` is the testthat expectation: it fails with the
report when any check fails.

## Examples

``` r
if (FALSE) { # \dontrun{
g <- GiottoData::loadGiottoMini("visium")
g2 <- Giotto::runNeighborhoodPCA(g, feats_to_use = "hvf")
checkGiottoMethod(g, g2, slots = "dimension_reduction",
    param_class = "neighborhoodPCAParam", fn = "runNeighborhoodPCA",
    package = "Giotto")
} # }
```
