# LLM backend for mapping repair

Configure the model
[`repairMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/repairMapping.md)
may ask to choose among repair options. Settings come from the arguments
or the environment: `GIOTTO_LLM_PROVIDER` (`"anthropic"`, `"ollama"` or
`"baseline"`), `GIOTTO_LLM_MODEL`, and for Anthropic `ANTHROPIC_API_KEY`
(or `ANTHROPIC_AUTH_TOKEN`), for Ollama `OLLAMA_HOST` (default
`http://localhost:11434`). Keys are read at call time and never stored
or logged.

`"baseline"` uses no model: it takes the top-ranked candidate path when
its score leads the next one by at least `min_lead`, and makes no param
choices. It is the reference an LLM backend has to beat.

## Usage

``` r
llmBackend(
  provider = Sys.getenv("GIOTTO_LLM_PROVIDER"),
  model = Sys.getenv("GIOTTO_LLM_MODEL"),
  host = Sys.getenv("OLLAMA_HOST", "http://localhost:11434"),
  max_tokens = 2048L,
  min_lead = 0.5
)
```

## Arguments

- provider:

  `"anthropic"`, `"ollama"` or `"baseline"`

- model:

  model id. Anthropic defaults to `"claude-haiku-4-5"`; Ollama needs one
  (e.g. `"qwen3.8:27b"`).

- host:

  Ollama host

- max_tokens:

  maximum response tokens

- min_lead:

  baseline only: required score lead of the top candidate

## Value

an `llmBackend` list
