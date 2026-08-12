# Subset an iconic_data object to a random panel for confounding inference

Draws a random subset of mediators and (when more than one) outcome
features, capping each at `max_infer_tasks`, and returns a thinned
`iconic_data` object suitable for
[`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md).
Used to make the `estimate = NULL` auto-run in
[`infer_confounding()`](https://seantbresnahan.com/iconic/reference/infer_confounding.md)
fast on large panels: the confounding-strength gaps are averages across
the mediator x feature grid, so a random subset is an unbiased Monte
Carlo estimate of the full-panel value.

## Usage

``` r
.subset_data_for_inference(data, max_infer_tasks = 50)
```

## Arguments

- data:

  An `iconic_data` object.

- max_infer_tasks:

  Cap on mediators and on features.

## Value

A thinned `iconic_data` object with an `inference_subset` attribute
recording the subset sizes.
