# Convert an iconic_data object to a MultiAssayExperiment

Exports the standardised analysis data held in an `iconic_data` object
as a
[MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html),
the Bioconductor container for multi-omic data on shared samples. Use it
to store, share, subset, or inspect the exact data an ICONIC analysis
was run on with standard Bioconductor tools. The result can be passed
straight back to
[`as_iconic_data()`](https://seantbresnahan.com/iconic/reference/as_iconic_data.md)
(the role mapping is recorded in the metadata).

## Usage

``` r
as_MultiAssayExperiment(data)
```

## Arguments

- data:

  An `iconic_data` object (from
  [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md)
  or
  [`as_iconic_data()`](https://seantbresnahan.com/iconic/reference/as_iconic_data.md)).

## Value

A `MultiAssayExperiment` with

- experiments:

  `SummarizedExperiment`s (assay `value`, features x samples) named `Y`
  (outcome panel; continuous outcomes only), `M` (mediators), `Gm`
  (mediator instruments), `W`, `W1`, `W2` (negative-control panels),
  each present only when supplied.

- colData:

  `X` (exposure), `G` (exposure instrument), `surv_time` / `surv_event`
  (survival), `Y_bin` (binary outcome), and the encoded covariates
  (prefixed `cov_` if a name collides).

- metadata:

  `iconic`: the outcome type, scaling parameters, feature / mediator
  names, availability flags, the role mapping used by
  [`as_iconic_data()`](https://seantbresnahan.com/iconic/reference/as_iconic_data.md),
  and the attached GAN (if any).

Values are those stored in `data`, i.e. centred and scaled when
`iconic_data(scale = TRUE)` (the default); original-scale parameters are
kept in `metadata(mae)$iconic$scaling`.

## Examples

``` r
d <- iconic_data(X = rnorm(60), Y = matrix(rnorm(60 * 5), 5, 60),
    M = rnorm(60), G = rnorm(60), Gm = rnorm(60),
    W = matrix(rnorm(60 * 5), 5, 60))
mae <- as_MultiAssayExperiment(d)
mae
#> A MultiAssayExperiment object of 4 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 4:
#>  [1] Y: SummarizedExperiment with 5 rows and 60 columns
#>  [2] M: SummarizedExperiment with 1 rows and 60 columns
#>  [3] Gm: SummarizedExperiment with 1 rows and 60 columns
#>  [4] W: SummarizedExperiment with 5 rows and 60 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
d2 <- as_iconic_data(mae) # round trip
all.equal(d$Y, d2$Y, check.attributes = FALSE)
#> [1] TRUE
```
