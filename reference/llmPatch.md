# Ask an LLM backend for a repair patch

Send one repair request and return the patches the backend chose, after
checking each against the request's options. Patches naming an option
that was not offered, or with `"low"` confidence, are dropped.

## Usage

``` r
llmPatch(request, backend)
```

## Arguments

- request:

  a repair request built by
  [`repairMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/repairMapping.md)
  (list with `items`, `prompt` and `schema`)

- backend:

  an
  [`llmBackend()`](https://giotto-suite.github.io/GiottoConnect/reference/llmBackend.md)

## Value

list with `patches` (data.frame: `item`, `value`, `confidence`,
`rationale`) and `log` (one-row data.frame: provider, model, prompt
characters, input/output tokens, latency, patches returned and kept,
error)
