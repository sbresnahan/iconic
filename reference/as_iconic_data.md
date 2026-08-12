# Convert external data containers to iconic_data

S3 generic bridging other data interfaces to the estimation interface.
Methods:

## Usage

``` r
as_iconic_data(input, ...)

# Default S3 method
as_iconic_data(input, ...)

# S3 method for class 'SummarizedExperiment'
as_iconic_data(
  input,
  assay = 1,
  mediator_assay = NULL,
  exposure,
  instrument = NULL,
  mediator_instrument = NULL,
  negative_controls = NULL,
  covariates = NULL,
  surv_time = NULL,
  surv_event = NULL,
  outcome_type = c("continuous", "survival"),
  ...
)
```

## Arguments

- input:

  An object to convert: a list returned by
  [`load_real_input_data()`](https://seantbresnahan.com/iconic/reference/load_real_input_data.md),
  an `iconic_data` object, an exposure vector (named-argument form), or
  a SummarizedExperiment.

- ...:

  Named arguments passed to
  [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md)
  when using the named-argument form, or to the method.

- assay:

  Name or index of the assay of `input` holding the outcome panel
  (features x samples). Default `1`. Set to `NULL` for survival outcomes
  (when `surv_time`/`surv_event` are colData columns and no continuous
  outcome panel is used).

- mediator_assay:

  Optional name or index of an assay holding the mediator panel
  (mediators x samples).

- exposure:

  Character: name of the `colData` column holding the exposure X.

- instrument:

  Optional character: name of the `colData` column holding the exposure
  instrument G (e.g. a polygenic score).

- mediator_instrument:

  Optional character vector of `colData` column names holding the
  mediator instrument(s) Gm (one column per mediator).

- negative_controls:

  Optional character vector of `colData` column names forming the
  negative-control panel W.

- covariates:

  Optional character vector of `colData` column names to carry through
  as covariates.

- surv_time, surv_event:

  Optional character: names of `colData` columns holding follow-up time
  and the 0/1 event indicator; set together with
  `outcome_type = "survival"`.

- outcome_type:

  `\"continuous\"` (default) or `\"survival\"`.

## Value

An `iconic_data` S3 object.

## Details

- `default`: a list returned by
  [`load_real_input_data()`](https://seantbresnahan.com/iconic/reference/load_real_input_data.md),
  an existing `iconic_data` object (returned as-is), or the exposure
  vector `X` with named arguments matching
  [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md)
  (Y, M, G, Gm, W, W1, W2, covariates, feature_names, mediator_names),
  which delegates to
  [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md).

- `SummarizedExperiment`: extracts the outcome panel from an assay and
  sample-level fields (exposure, instruments, negative controls,
  covariates) from `colData`. Requires the SummarizedExperiment package
  (listed under `Suggests`).

## Examples

``` r
# From a load_real_input_data() result
input <- load_real_input_data(example = TRUE)
data <- as_iconic_data(input)
#> Warning: W has 30 features but Y has 1. Using row recycling.
print(data)
#> <iconic_data> 200 samples, 1 outcome features
#>  Available: W (negative controls), W1/W2 (path-specific NCs) 
#>  Covariates: sex, GA, mother_ethnicity_indian, mother_ethnicity_malay 
#>  Mode: total effect 

# From named components (delegates to iconic_data())
data <- as_iconic_data(rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
G = rnorm(100), Gm = rnorm(100),
W = matrix(rnorm(100*10), 10, 100))
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = matrix(rnorm(20 * 60), 20, 60,
                                dimnames = list(paste0("gene", 1:20),
                                                paste0("S", 1:60)))),
    colData = S4Vectors::DataFrame(
      bmi = rnorm(60), prs = rnorm(60),
      nc1 = rnorm(60), nc2 = rnorm(60), age = rnorm(60))
  )
  data <- as_iconic_data(se, assay = "expr", exposure = "bmi",
                         instrument = "prs",
                         negative_controls = c("nc1", "nc2"),
                         covariates = "age")
  print(data)
}
#> Warning: replacing previous import ‘S4Arrays::makeNindexFromArrayViewport’ by ‘DelayedArray::makeNindexFromArrayViewport’ when loading ‘SummarizedExperiment’
#> Warning: W has 2 features but Y has 20. Using row recycling.
#> <iconic_data> 60 samples, 20 outcome features
#>  Available: G (exposure instrument), W (negative controls), W1/W2 (path-specific NCs) 
#>  Covariates: age 
#>  Mode: total effect 
```
