# Run Scrublet doublet detection

Runs Scrublet (Wolock, Lopez and Klein 2019) on the counts of a giotto
object and stores each cell's doublet score and call as cell metadata,
in columns `<name>_score` and `<name>_doublet`. When Scrublet cannot set
a doublet score threshold and the param gives none, only the scores are
stored, with a warning. Scrublet is meant for single-cell data, not for
spots that hold several cells. **If you use it, please cite the original
paper** (see
[`scrubletParam()`](https://giotto-suite.github.io/GiottoConnect/reference/scrubletParam.md)).

## Usage

``` r
runScrublet(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = "raw",
  param = NULL,
  name = "scrublet",
  verbose = NULL
)
```

## Arguments

- gobject:

  giotto object

- spat_unit, feat_type:

  spatial unit and feature type

- expression_values:

  expression values to use; Scrublet expects counts

- param:

  a
  [`scrubletParam()`](https://giotto-suite.github.io/GiottoConnect/reference/scrubletParam.md);
  its defaults when `NULL`

- name:

  prefix of the cell metadata columns

- verbose:

  be verbose

## Value

The giotto object with Scrublet's results in its cell metadata.

## References

Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of
cell doublets in single-cell transcriptomic data. *Cell Systems* 8,
281-291.e9 (2019).
[doi:10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005)

## See also

[`scrubletParam()`](https://giotto-suite.github.io/GiottoConnect/reference/scrubletParam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
g <- GiottoData::loadGiottoMini("visium")
g <- runScrublet(g)
pDataDT(g)[, c("cell_ID", "scrublet_score", "scrublet_doublet")]
} # }
```
