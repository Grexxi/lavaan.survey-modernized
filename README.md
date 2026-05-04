# lavaan.survey

This repository is an unofficial modernization fork of the original
`lavaan.survey` package by Daniel Oberski. It is not an official CRAN release
and has not yet been reviewed or endorsed by the original author/maintainer.
The original copyright and GPL (>= 2) license are preserved.

`lavaan.survey` fits structural equation models with complex survey designs by
combining `lavaan` model fitting with design information from the `survey`
package.

This 1.2.0 modernization snapshot updates the original package for current
`lavaan` versions and adds initial support for ordinal survey SEM through
`lavaan.survey.ordinal()`.

## Installation

Install the local source package with:

```r
install.packages("lavaan.survey_1.2.0.tar.gz", repos = NULL, type = "source")
```

Then load it as usual:

```r
library(lavaan.survey)
```

## Continuous SEM with complex survey designs

The original workflow is unchanged: fit the model in `lavaan`, define the survey
design with `survey::svydesign()`, then call `lavaan.survey()`.

```r
library(lavaan)
library(survey)
library(lavaan.survey)

fit <- cfa(model, data = dat, estimator = "MLM")

des <- svydesign(
  ids = ~psu,
  strata = ~stratum,
  weights = ~weight,
  data = dat,
  nest = TRUE
)

fit_svy <- lavaan.survey(fit, des)
summary(fit_svy, fit.measures = TRUE, standardized = TRUE)
```

## Continuous SEM with multiple imputation

For continuous models, `lavaan.survey()` also supports survey designs whose
`data` argument is a `mitools::imputationList`. Fit the naive `lavaan` model on
one completed dataset, then pass the imputed survey design to `lavaan.survey()`.
The sample means, covariances, and their design-based covariance matrix are
pooled across imputations using Rubin-style pooling.

```r
library(mitools)

imputed_data <- list(dat_imp1, dat_imp2, dat_imp3)
imp <- imputationList(imputed_data)

des_imp <- svydesign(
  ids = ~psu,
  strata = ~stratum,
  weights = ~weight,
  data = imp,
  nest = TRUE
)

fit_imp <- cfa(model, data = imputed_data[[1]], estimator = "MLM")
fit_svy_imp <- lavaan.survey(fit_imp, des_imp)
summary(fit_svy_imp, fit.measures = TRUE, standardized = TRUE)
```

As with any survey-weighted MI analysis, the imputation model should be chosen
with care and should reflect important design information such as weights,
clusters, and strata.

## Ordinal SEM with complex survey designs

For ordinal indicators, fit the naive ordinal model first with `lavaan`, then
pass the fitted object and the survey design to `lavaan.survey.ordinal()`.

```r
items <- paste0("v", 1:10)

model <- "
  f1 =~ v1 + v2 + v3 + v4 + v5
  f2 =~ v6 + v7 + v8 + v9 + v10
"

fit_naive <- cfa(
  model = model,
  data = dat,
  ordered = items,
  estimator = "WLSMV"
)

des <- svydesign(
  ids = ~cluster,
  strata = ~stratum,
  weights = ~weight,
  data = dat,
  nest = TRUE
)

fit_svy <- lavaan.survey.ordinal(
  lavaan.fit = fit_naive,
  survey.design = des,
  estimator = "WLSMV"
)

summary(fit_svy, fit.measures = TRUE, standardized = TRUE)
```

Internally, traditional `svydesign` objects are converted to replicate-weight
designs with `survey::as.svrepdesign()`. Thresholds and polychoric correlations
are estimated with sampling weights, and their design-based covariance matrix is
estimated from the replicate weights.

## Multiple-group ordinal models

Multiple-group ordinal models are supported when all observed model variables
are ordered.

```r
fit_group <- cfa(
  model = model,
  data = dat,
  ordered = items,
  group = "gender",
  estimator = "WLSMV"
)

fit_group_svy <- lavaan.survey.ordinal(fit_group, des)
summary(fit_group_svy, fit.measures = TRUE, standardized = TRUE)
```

## Ordinal measurement invariance

Equality constraints from the original `lavaan` fit are preserved. For example,
loading and threshold invariance can be fitted with:

```r
fit_scalar <- cfa(
  model = model,
  data = dat,
  ordered = items,
  group = "gender",
  group.equal = c("loadings", "thresholds"),
  estimator = "WLSMV"
)

fit_scalar_svy <- lavaan.survey.ordinal(fit_scalar, des)
summary(fit_scalar_svy, fit.measures = TRUE, standardized = TRUE)
```

## Current limitations

`lavaan.survey.ordinal()` is an initial implementation. It currently supports
single-group and multiple-group models where all observed model variables are
ordered. Multiple imputation and mixed continuous/ordinal observed-variable sets
are not yet implemented for the ordinal workflow.

## Package checks

The package can be checked with:

```sh
R CMD check lavaan.survey_1.2.0.tar.gz
```

## Cross-software validation

The repository includes a small Mplus Demo validation workflow under
`validation/mplus-demo/`. It creates a six-indicator ordinal survey CFA that fits
within the Mplus Demo limits and compares the resulting Mplus output with
`lavaan.survey.ordinal()` results. For WLSMV, the Mplus fit block is compared
primarily with lavaan's scaled fit measures; lavaan's additional robust CFI, TLI,
and RMSEA are treated as sensitivity checks rather than direct Mplus targets.
