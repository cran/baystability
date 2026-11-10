

# baystability 1.0.0

## Changes visible to users

- `bayes_ammi()` keeps its arguments (`.data`, `.y`, `.gen`, `.env`,
  `.rep`, `.nIter`) and the components of its result (`mu1`, `tau1`,
  `tao1`, `delta1`, `lambdas1`, `alphas1`, `gammas1`, `Mean_Alphas`,
  `Mean_Gammas`), and gains:
  - `.prior_data`: the earlier trial used as prior (Perez-Elizalde et
    al. 2012, Section 2.2), and `.n0`, its weight (default its
    replicates); without it the priors are vague (`bys_priors()`);
  - `.model = "SREG"` (sites regression, GGE) besides `"AMMI"`;
  - `.k` (number of terms; default all, as before), `.burnin`, `.thin`,
    `.chains`, `.seed`, `.level`, `.align`, `.waic`;
  - the components `summary` (posterior mean, SD, Monte Carlo error, HPD
    interval, effective sample size, potential scale reduction factor),
    `scores`, `draws` (`coda::mcmc.list`), `diagnostics`, `ls`, `prior`,
    `fit` (DIC, WAIC) and `settings`;
  - `print()`, `summary()` and `plot()` methods.
- New names that follow the package **stability** and its Bayesian
  version **bayesstability**: the prefix `bayes_` for Bayesian analyses
  and `bys_` for tools. The least-squares tables `ge_mean()`,
  `ge_var()`, `g_eff()`, `e_eff()`, `ge_eff()`, `ge_model()` and
  `ge_ammi()` are `bys_ls_means()`, `bys_ls_var()`,
  `bys_ls_gen_effects()`, `bys_ls_env_effects()`, `bys_ls_ge_effects()`,
  `bys_ls_model()` and `bys_ls_ammi()`; `matrix_k()` and `orthnorm()`
  are `bys_contrasts()` and `bys_orthonormal()`. The old names still
  work, with a deprecation warning, and will be removed in a later
  version (`?"bys-deprecated"`). Every analysis is an S3 generic with a
  `.default` method.
- New `bayes_gge()` (the SREG model, `bayes_ammi(.model = "SREG")`),
  `bayes_ammi_biplot()` and `bayes_gge_biplot()` (HPD regions),
  `bys_priors()`, `bys_compare()`, `bys_interpret()`, `bys_cite()`;
  every printed result ends with the sources of its method.
- New `bys_app()`, a built-in ‘shiny’ application: data and earlier
  trial, model, biplot, model choice and references, each page with an
  interpretation, tables, plots and the R code that repeats the fit.
- The sampler is compiled (Rcpp, RcppArmadillo); a plain R version gives
  the same draws from the same seed.
- Six vignettes in HTML and PDF, among them a reproduction of the
  published analysis.

## Corrections

- `bayes_ammi()` read the earlier trial from objects in the user’s
  workspace (`n0`, `mu0`, `tau0`, `tao0`, `delta0`, `lambdas0`,
  `alphas0`, `gammas0`, `ge_means0`) and stopped with “object ‘n0’ not
  found” without them.
- The singular values were computed once before the loop and never
  sampled.
- `tau1` held a precision drawn from a residual matrix computed once
  before the loop and never updated; it now holds the draws of the error
  variance sigma^2.
- The scores were drawn without their constraints and orthogonalised
  afterwards; they are now drawn from their von Mises-Fisher
  conditionals on the constrained spheres.
- With all terms fitted (11 terms for 12 genotypes) a column-wise
  sampler cannot move the genotype scores: each column is fixed by the
  others up to its sign. Exact moves of every pair of columns (von Mises
  draws of a rotation angle) now let them mix.
- The columns of `delta1` were named after the genotypes.
- The help page of `cultivo2009` described another data set.

## Reproduction of Perez-Elizalde et al. (2012)

- The 38 linear effects of Table B.1 agree to within 0.03 posterior
  standard deviations, their standard deviations to within 2%, and sigma
  of Table 2 to within 0.2 standard deviations. The singular values of
  Table 2 are not reproduced, and neither scores that cannot move nor a
  zero prior mean explains them: the published values disagree with eq.
  2.15 of the paper (see the vignette on the published analysis).

## Internal

- The tidyverse, `rlang`, `rstiefel` and `lme4` (now suggested, only for
  `bys_ls_model()`) are no longer needed.

# baystability 0.2.0

- Updated the DESCRIPTION file to meet the requirements of CRAN.

# baystability 0.1.0

- First version of baystability.
