# Convert a linear predictor to survival data

Given an n-vector (or n x p matrix) of linear predictors `eta` and a
target event fraction, returns a list with `surv_time` and `surv_event`
vectors. When `eta` is a matrix, the first column is used (survival
outcomes are scalar — n_features = 1).

## Usage

``` r
.linpred_to_surv(eta, h0 = 0.1, event_frac = 0.6, censor_rate = NULL)
```

## Arguments

- eta:

  numeric vector or matrix of linear predictors.

- h0:

  baseline hazard (default 0.1).

- event_frac:

  target event fraction; `censor_rate` is tuned so that approximately
  this fraction of subjects are observed events.

- censor_rate:

  optional explicit censoring rate; if NULL it is solved from
  `event_frac`.

## Value

list(surv_time, surv_event, true_h0, true_censor_rate)
