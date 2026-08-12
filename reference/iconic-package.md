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
data. Helper functions construct exposure instruments (polygenic scores
from GWAS summary statistics or PGS Catalog scoring files), mediator
instruments (cis-eQTL scans and elastic-net predicted-expression
composites), and negative-control panels (principal components of
residualized omics matrices), and data can be imported directly from
SummarizedExperiment objects. Time-to-event outcomes are supported on
the Cox log-hazard-ratio and restricted-mean-survival-time scales. The
torch package is required for the generative texture model (sensitivity
and prospective analysis); all other functionality works without it.

## See also

Useful links:

- <https://github.com/sbresnahan/iconic>

- Report bugs at <https://github.com/sbresnahan/iconic/issues>

## Author

**Maintainer**: Sean T. Bresnahan <seantbresnahan3@gmail.com>

Authors:

- Sean T. Bresnahan <seantbresnahan3@gmail.com>

- Charis Xiong <charis.xiong@gmail.com>
