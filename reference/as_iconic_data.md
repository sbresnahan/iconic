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
  bin_outcome = NULL,
  outcome_type = c("continuous", "survival", "binary"),
  ...
)

# S3 method for class 'MultiAssayExperiment'
as_iconic_data(
  input,
  exposure,
  outcome = NULL,
  mediator = NULL,
  instrument = NULL,
  mediator_instrument = NULL,
  negative_controls = NULL,
  W1 = NULL,
  W2 = NULL,
  covariates = NULL,
  surv_time = NULL,
  surv_event = NULL,
  bin_outcome = NULL,
  outcome_type = c("continuous", "survival", "binary"),
  assay = 1L,
  verbose = TRUE,
  ...
)
```

## Arguments

- input:

  An object to convert: a MultiAssayExperiment, a SummarizedExperiment,
  a list returned by
  [`load_real_input_data()`](https://seantbresnahan.com/iconic/reference/load_real_input_data.md),
  an `iconic_data` object, or an exposure vector (named-argument form).

- ...:

  Named arguments passed to
  [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md)
  when using the named-argument form, or to the method.

- assay:

  SummarizedExperiment method: name or index of the assay of `input`
  holding the outcome panel (features x samples). Default `1`. Set to
  `NULL` for survival outcomes (when `surv_time`/`surv_event` are
  colData columns and no continuous outcome panel is used).
  MultiAssayExperiment method: name or index of the assay used within
  each SummarizedExperiment experiment. Default `1`.

- mediator_assay:

  Optional name or index of an assay holding the mediator panel
  (mediators x samples).

- exposure:

  Character: name of the `colData` column holding the exposure X
  (MultiAssayExperiment method: or of a single experiment). For the
  MultiAssayExperiment method, when `exposure` is missing and `input`
  was produced by
  [`as_MultiAssayExperiment()`](https://seantbresnahan.com/iconic/reference/as_MultiAssayExperiment.md),
  the roles recorded in `metadata(input)$iconic` are used (round trip).

- instrument:

  Optional character: name of the `colData` column holding the exposure
  instrument G (e.g. a polygenic score).

- mediator_instrument:

  Optional character vector of `colData` column names holding the
  mediator instrument(s) Gm (one column per mediator); for the
  MultiAssayExperiment method, alternatively the name of an experiment
  (instruments x samples).

- negative_controls:

  Optional character vector of `colData` column names forming the
  negative-control panel W; for the MultiAssayExperiment method,
  alternatively the name of an experiment (e.g. principal components
  from
  [`build_w_pcs()`](https://seantbresnahan.com/iconic/reference/build_w_pcs.md),
  features x samples).

- covariates:

  Optional character vector of `colData` column names to carry through
  as covariates.

- surv_time, surv_event:

  Optional character: names of `colData` columns holding follow-up time
  and the 0/1 event indicator; set together with
  `outcome_type = "survival"`.

- bin_outcome:

  Optional character: name of a `colData` column holding the 0/1
  outcome; set together with `outcome_type = "binary"`.

- outcome_type:

  `\"continuous\"` (default), `\"survival\"`, or `\"binary\"`.

- outcome:

  For the MultiAssayExperiment method: name of the experiment holding
  the outcome panel Y (features x samples). Required when
  `outcome_type = "continuous"`.

- mediator:

  Optional: name of the experiment (or character vector of `colData`
  columns) holding the mediator panel M.

- W1, W2:

  Optional: names of experiments (or character vectors of `colData`
  columns) holding the path-specific negative-control panels for the
  X-\>M (`W1`) and M-\>Y (`W2`) paths.

- verbose:

  Logical: report samples dropped during alignment. Default `TRUE`.

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

- `MultiAssayExperiment` (recommended for multi-omic data): omics panels
  (outcome, mediator, negative controls, mediator instruments) are
  experiments, sample-level fields (exposure, instrument, covariates,
  survival / binary outcome) are `colData` columns, and samples are
  matched through the `sampleMap`. See Details and
  [`as_MultiAssayExperiment()`](https://seantbresnahan.com/iconic/reference/as_MultiAssayExperiment.md)
  for the reverse conversion.

- `SummarizedExperiment`: extracts the outcome panel from an assay and
  sample-level fields (exposure, instruments, negative controls,
  covariates) from `colData`.

**MultiAssayExperiment method.** Each role (`outcome`, `mediator`,
`negative_controls`, `W1`, `W2`, `mediator_instrument`, and optionally
`exposure` / `instrument`) names either one experiment of `input` or one
or more `colData` columns; experiment names take precedence. For
SummarizedExperiment experiments the assay selected by `assay` is used;
matrix experiments are used directly. Samples are matched through
`sampleMap(input)` to primary IDs (the rownames of `colData`). Only
primaries present in every experiment used, and with non-missing values
in the sample-level `colData` columns used, are kept. Experiments with
more than one column per primary (e.g. technical replicates) must be
collapsed first, for example with
[`MultiAssayExperiment::mergeReplicates()`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment-helpers.html).
The aligned components are passed to
[`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md),
so validation and scaling are identical to calling it directly.

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
data <- as_iconic_data(rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
    G = rnorm(100), Gm = rnorm(100),
    W = matrix(rnorm(100 * 10), 10, 100))
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
#> Warning: W has 2 features but Y has 20. Using row recycling.
print(data)
#> <iconic_data> 60 samples, 20 outcome features
#>  Available: G (exposure instrument), W (negative controls), W1/W2 (path-specific NCs) 
#>  Covariates: age 
#>  Mode: total effect 
# From a MultiAssayExperiment: omics panels are experiments,
# sample-level variables are colData columns
set.seed(1)
ids <- paste0("S", 1:80)
mk <- function(p, pre) matrix(rnorm(p * 80), p, 80,
    dimnames = list(paste0(pre, seq_len(p)), ids))
mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(
        rna = SummarizedExperiment::SummarizedExperiment(
            assays = list(expr = mk(10, "gene"))),
        methyl = SummarizedExperiment::SummarizedExperiment(
            assays = list(m = mk(1, "cg"))),
        nc = SummarizedExperiment::SummarizedExperiment(
            assays = list(pcs = mk(10, "PC")))),
    colData = S4Vectors::DataFrame(
        bmi = rnorm(80), prs = rnorm(80), cis_score = rnorm(80),
        age = rnorm(80), row.names = ids))
data <- as_iconic_data(mae, exposure = "bmi", outcome = "rna",
    mediator = "methyl", instrument = "prs",
    mediator_instrument = "cis_score", negative_controls = "nc",
    covariates = "age")
print(data)
#> <iconic_data> 80 samples, 10 outcome features, 1 mediator(s)
#>  Available: G (exposure instrument), Gm (mediator instrument), W (negative controls), W1/W2 (path-specific NCs) 
#>  Covariates: age 
#>  Mode: mediation 
```
