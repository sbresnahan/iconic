# iconic: Causal Model Selection with Genetic Instruments and Negative Controls

Provides a model selection workflow for causal inference with genetic
instruments and negative controls in observational omics data. The
package fits eight estimators of the natural direct and indirect effects
(NDE/NIE), diagnoses which are valid for the user's data via
instrument-strength checks and negative-control validity screens (A1,
A2, A2'), stress-tests them against confounding and pleiotropy
violations, and recommends the estimator most likely to be unbiased. A
hybrid generative texture model (a torch GAN for sample-level structure
and a Gaussian copula for the mediator panel) lets the sensitivity
analysis mirror the marginal and joint structure of the user's cohort.
The modular API comprises iconic_data() to standardize data,
iconic_diagnose() to assess estimator eligibility, iconic_estimate() to
fit all eligible estimators, iconic_sensitivity() to sweep instrument
exogeneity, iconic_recommend() to rank estimators by identification
strength and robustness, and iconic_prospect() to plan future studies.
infer_confounding() estimates held-fixed confounding parameters from the
data. Multi-omic cohorts are supplied as MultiAssayExperiment objects
(omics panels as experiments, sample-level exposure, instruments and
covariates as colData), and analysis data can be exported back to
MultiAssayExperiment. Helper functions construct exposure instruments
(polygenic scores from GWAS summary statistics or PGS Catalog scoring
files), mediator instruments (GenomicRanges-based cis-eQTL scans and
elastic-net predicted-expression composites from
RangedSummarizedExperiment or matrix inputs), and negative-control
panels (principal components of residualized omics matrices).
Time-to-event outcomes are supported on the Cox log-hazard-ratio and
restricted-mean-survival-time scales. Binary outcomes are supported on
the log-odds-ratio (logistic two-stage predictor-substitution outcome
stage) and risk-difference (linear probability model) scales. The torch
package is required for the generative texture model (sensitivity and
prospective analysis); all other functionality works without it. The
methods are described in Bresnahan et al. (2026)
[doi:10.64898/2026.08.26.26361466](https://doi.org/10.64898/2026.08.26.26361466)
.

## References

Bresnahan ST, Xiong C, Head ST, Chang Y-H, Bhattacharya A, Huang JY.
ICONIC: An R Package for Integrating Instrumental Variable- and
Negative-Control-Informed Causal Discovery and Diagnostics in Multiomic
Studies. medRxiv (2026).
[doi:10.64898/2026.08.26.26361466](https://doi.org/10.64898/2026.08.26.26361466)

## See also

Useful links:

- <https://seantbresnahan.com/iconic/>

- <https://github.com/sbresnahan/iconic>

- Report bugs at <https://github.com/sbresnahan/iconic/issues>

## Author

**Maintainer**: Sean T. Bresnahan <seantbresnahan3@gmail.com>
([ORCID](https://orcid.org/0000-0001-6685-1930))

Authors:

- Sean T. Bresnahan <seantbresnahan3@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-6685-1930))

- Charis Xiong <charis.xiong@gmail.com>
