# Mplus Demo Validation Results

This note records cross-software validation checks for the modernized
`lavaan.survey` package against Mplus Demo 9, plus the continuous package
regression tests that anchor the original workflow.

## Note on Scaled and Robust Fit Measures

Mplus WLSMV output reports one adjusted model-fit block for the WLSMV estimator.
It does not print a second lavaan-style column of `robust` CFI, TLI, and RMSEA.
For this validation, the closest direct comparison is therefore Mplus WLSMV
against lavaan's `*.scaled` fit measures. lavaan's `*.robust` measures are shown
as an additional sensitivity check, not as the primary Mplus target.

Mplus marks WLSMV chi-square values with a note that they cannot be used for
ordinary chi-square difference testing and points WLSMV users to `DIFFTEST` for
nested model comparisons. The Mplus authors also describe the WLSMV/ULSMV/MLMV
second-order chi-square correction as a scaled and shifted approximation to a
chi-square distribution.

## Validation Ladder

The validation results are organized from the original continuous workflow to the
new ordinal workflow:

| Step | Target | Result status |
| --- | --- | --- |
| 1 | Continuous survey SEM | Covered by package regression tests and the bundled Roosma example; no separate Mplus Demo result section yet. |
| 2 | Continuous survey SEM with multiple imputation | Cross-checked against Mplus Demo. |
| 3 | Continuous multiple-group / invariance models | Planned next validation target. |
| 4 | Ordinal survey SEM | Cross-checked against Mplus Demo on simulated data and ESS4 GB. |
| 5 | Ordinal survey SEM with multiple imputation | Cross-checked against Mplus Demo. |
| 6 | Ordinal multiple-group / invariance models | Cross-checked against Mplus Demo. |
| 7 | Ordinal multiple-group / invariance models with multiple imputation | Cross-checked against Mplus Demo. |
| 8 | Mixed ordinal/continuous multiple-group models with multiple imputation | Diagnostic workflow runs in Mplus Demo and compares pooled sample statistics with experimental Rubin parameter pooling. |


## Validation 1: Continuous Survey SEM

The basic continuous survey SEM path is covered by package regression tests
rather than by a separate Mplus Demo workflow in this folder.

- `tests/testthat/test_roosma.R` checks the bundled ESS4/Roosma example against
  legacy numerical targets for `lavaan.survey()`.
- `tests/testthat/test_continuous_robustness.R` checks modern continuous-path
  robustness cases, including Yuan-Bentler without a mean structure,
  multiple-group refits, multiple imputation sample sizes, and `pval.pFsum()`.

A dedicated continuous Mplus Demo comparison without multiple imputation remains
a useful future validation target.


## Validation 2: Continuous Multiple-Imputation Survey CFA

### Scope

This validation checks the modernized continuous multiple-imputation path in the
original `lavaan.survey()` function. The script simulates a two-factor
continuous CFA model with six indicators, 600 observations, 60 clusters, 6
strata, and sampling weights. Artificial missingness is introduced in three
indicators:

```text
y2  MAR, depending on y1 and stratum
y5  MAR, depending on y4 and weight
y6  MCAR, about 15 percent
```

The missing data are imputed ten times with `mice`. The same ten completed
datasets are then analyzed by both `lavaan.survey()` and Mplus Demo, so the
comparison focuses on the SEM/survey/MI analysis layer rather than on different
imputation engines.

The model is:

```text
f1 =~ y1 + y2 + y3
f2 =~ y4 + y5 + y6
f1 ~~ f2
```

The corresponding Mplus input uses:

```text
DATA: TYPE = IMPUTATION;
VARIABLE: WEIGHT IS wgt; CLUSTER IS clu; STRATIFICATION IS str;
ANALYSIS: TYPE = COMPLEX; ESTIMATOR = MLR;
```

### Commands

Prepare the imputed datasets, Mplus input, and `lavaan.survey()` results:

```r
source("validation/mplus-demo/prepare_continuous_mi_validation_files.R")
```

Run Mplus Demo from `validation/mplus-demo`:

```sh
/Applications/MplusDemo/mpdemo continuous_mi_complex.inp
```

Parse and compare the output:

```r
source("validation/mplus-demo/compare_mplus_continuous_mi_output.R")
```

### Fit Measures

| Measure | lavaan scaled | Mplus MLR imputation mean |
| --- | ---: | ---: |
| Chi-square | 8.875 | 15.327 |
| df | 8 | 8 |
| p-value | 0.353 | -- |
| CFI | 0.999 | 0.995 |
| TLI | 0.998 | 0.990 |
| RMSEA | 0.014 | 0.036 |
| SRMR | 0.0205 | 0.022 |

Mplus prints the mean and standard deviation over the ten imputed-data fit
statistics; it does not print the same single pooled lavaan-style scaled test
statistic used by `lavaan.survey()`.

### Parameter Agreement

The comparison matched 21 unstandardized parameters between the two outputs.

| Quantity | Value |
| --- | ---: |
| Matched parameters | 21 |
| Maximum absolute estimate difference | 0.0012 |
| Maximum absolute standard error difference | 0.0015 |

Largest estimate differences:

| Parameter | lavaan.survey | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `y2 ~~ y2` | 0.3542 | 0.3530 | -0.0012 |
| `f2 ~~ f2` | 0.7452 | 0.7460 | 0.0008 |
| `f2 =~ y6` | 0.8383 | 0.8390 | 0.0007 |
| `f1 ~~ f2` | 0.2865 | 0.2860 | -0.0005 |
| `f2 =~ y5` | 0.8645 | 0.8650 | 0.0005 |
| `y4 ~1` | -0.1065 | -0.1060 | 0.0005 |

Largest standard-error differences:

| Parameter | lavaan.survey | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `f2 =~ y5` | 0.0645 | 0.0630 | -0.0015 |
| `f1 =~ y3` | 0.0524 | 0.0510 | -0.0014 |
| `y2 ~~ y2` | 0.0376 | 0.0390 | 0.0014 |
| `y5 ~~ y5` | 0.0472 | 0.0460 | -0.0012 |
| `f1 =~ y2` | 0.0670 | 0.0680 | 0.0010 |

### Interpretation

This is a strong sanity check for the continuous MI modernization. When both
programs receive the same ten imputed datasets and the same complex survey
variables, unstandardized parameter estimates and standard errors agree to
roughly the third decimal place.

The fit measures are less direct. `lavaan.survey()` pools the sample statistics
and their design-based covariance matrix, then fits one lavaan model. Mplus
with `TYPE = IMPUTATION` runs the model across completed datasets and reports
means for several fit quantities. The parameter agreement is therefore the more
important validation target here.

## Validation 3: Continuous Multiple-Group / Invariance Survey CFA

This validation cell is intentionally listed as planned. The repository already
checks continuous multiple-group refitting in `test_continuous_robustness.R`,
but it does not yet contain a dedicated Mplus Demo comparison for continuous
measurement invariance.

A future workflow should mirror the ordinal invariance validation with
continuous indicators, survey weights, clusters, strata, and equality
constraints such as loadings and intercepts. This would complete the continuous
side of the validation ladder before extending it further.


## Validation 4a: Simulated Ordinal Survey CFA

### Scope

The Mplus Demo version is limited to six dependent variables. The validation
therefore uses a compact ordinal survey CFA model with six four-category
indicators and two factors:

```text
f1 =~ y1 + y2 + y3
f2 =~ y4 + y5 + y6
f1 ~~ f2
```

The generated dataset has 600 observations, survey weights, 60 clusters, and 6
strata. Both Mplus and `lavaan.survey.ordinal()` use WLSMV/DWLS-style estimation
for ordered categorical indicators with a complex survey correction.

### Commands

Prepare the validation files and `lavaan.survey.ordinal()` results:

```r
source("validation/mplus-demo/prepare_validation_files.R")
```

Run Mplus Demo from `validation/mplus-demo`:

```sh
/Applications/MplusDemo/mpdemo ordinal_survey_complex.inp
```

Parse and compare the output:

```r
source("validation/mplus-demo/compare_mplus_output.R")
```

### Fit Measures

| Measure | lavaan scaled | lavaan robust | Mplus WLSMV |
| --- | ---: | ---: | ---: |
| Adjusted chi-square | 8.481 | -- | 7.925 |
| df | 8 | -- | 8 |
| p-value | 0.388 | -- | 0.441 |
| CFI | 1.000 | 1.000 | 1.000 |
| TLI | 1.000 | 1.004 | 1.000 |
| RMSEA | 0.010 | 0.000 | 0.000 |
| SRMR | 0.0206 | 0.0206 | 0.015 |

### Parameter Agreement

The comparison matched 27 parameters between the two outputs.

| Quantity | Value |
| --- | ---: |
| Matched parameters | 27 |
| Maximum absolute estimate difference | 0.0028 |
| Maximum absolute standard error difference | 0.0028 |

Largest estimate differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `f1 ~~ f1` | 0.6572 | 0.6600 | 0.0028 |
| `f1 =~ y2` | 1.0614 | 1.0590 | -0.0024 |
| `f2 =~ y6` | 0.8429 | 0.8450 | 0.0021 |
| `f2 ~~ f2` | 0.7309 | 0.7290 | -0.0019 |
| `f2 =~ y5` | 0.8892 | 0.8910 | 0.0018 |
| `f1 =~ y3` | 0.9017 | 0.9000 | -0.0017 |
| `f1 ~~ f2` | 0.2519 | 0.2510 | -0.0009 |
| `y3 | t3` | 0.7745 | 0.7750 | 0.0005 |

Threshold estimates agreed especially closely; most differences are at the
third or fourth decimal place.

### Interpretation

This first validation run is encouraging. The core parameter estimates,
thresholds, and standard errors from `lavaan.survey.ordinal()` closely reproduce
the Mplus Demo results for a demo-compatible ordinal survey CFA. The fit
statistics are also in the same range.

Exact equality should not be expected. Mplus and `lavaan` differ in internal
rounding, weight scaling, robust test statistic corrections, and output
precision. The goal of this validation is therefore numerical agreement within a
small tolerance, not bit-for-bit identity.

## Validation 4b: ESS4 GB Ordinal Survey CFA

### Data Source

The `ess4.gb` dataset is bundled with `lavaan.survey`. Its package
documentation describes it as European Social Survey (ESS) round 4 data from the
2008 United Kingdom sample, downloaded from the ESS data portal and converted to
an R dataset. The variables used here measure attitudes about government
responsibility for welfare-state tasks and were used in the Roosma, Gelissen,
and van Oorschot (2013) welfare-state-attitudes application.

The original six selected ESS variables use 0-10 response scales. Directly
using all 11 categories produced sparse-category numerical instability in the
replicate-weight polychoric correlations. For this Mplus Demo validation, the
six original variables are therefore collapsed to four ordered categories:

```text
0-4   -> 1
5-6   -> 2
7-8   -> 3
9-10  -> 4
```

The variables and one-factor structure are unchanged; only the response scale is
collapsed for numerical stability.

### Scope

The ESS4 GB validation uses 2,194 complete cases, ESS design weights, PSUs, and
strata:

```text
range =~ gvjbevn + gvhlthc + gvslvol + gvslvue + gvcldcr + gvpdlwk
```

The corresponding Mplus input uses:

```text
WEIGHT IS dweight;
CLUSTER IS psu;
STRATIFICATION IS strat;
TYPE = COMPLEX;
ESTIMATOR = WLSMV;
PARAMETERIZATION = DELTA;
```

### Commands

Prepare the ESS4 validation files and `lavaan.survey.ordinal()` results:

```r
source("validation/mplus-demo/prepare_ess4_validation_files.R")
```

Run Mplus Demo from `validation/mplus-demo`:

```sh
/Applications/MplusDemo/mpdemo ess4_range_complex.inp
```

Parse and compare the output:

```r
source("validation/mplus-demo/compare_mplus_ess4_output.R")
```

### Fit Measures

| Measure | lavaan scaled | lavaan robust | Mplus WLSMV |
| --- | ---: | ---: | ---: |
| Adjusted chi-square | 223.100 | -- | 290.691 |
| df | 9 | -- | 9 |
| p-value | 0.000 | -- | 0.000 |
| CFI | 0.920 | 0.854 | 0.923 |
| TLI | 0.866 | 0.756 | 0.872 |
| RMSEA | 0.104 | 0.184 | 0.119 |
| SRMR | 0.0657 | 0.0657 | 0.050 |

### Parameter Agreement

The comparison matched 25 parameters between the two outputs.

| Quantity | Value |
| --- | ---: |
| Matched parameters | 25 |
| Maximum absolute estimate difference | 0.0615 |
| Maximum absolute standard error difference | 0.0089 |

Largest estimate differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `range =~ gvhlthc` | 1.2015 | 1.2630 | 0.0615 |
| `range =~ gvslvol` | 1.2687 | 1.3050 | 0.0363 |
| `range =~ gvcldcr` | 1.1754 | 1.1910 | 0.0156 |
| `range =~ gvslvue` | 1.0249 | 1.0390 | 0.0141 |
| `range ~~ range` | 0.3470 | 0.3440 | -0.0030 |
| `range =~ gvpdlwk` | 1.1533 | 1.1560 | 0.0027 |
| `gvhlthc | t1` | -2.0825 | -2.0820 | 0.0005 |
| `gvpdlwk | t3` | 0.6504 | 0.6500 | -0.0004 |

Threshold estimates again agree very closely. The larger differences are mainly
in a few factor loadings for highly skewed ESS items.

### Interpretation

The ESS4 GB comparison is a harder real-data check than the simulated
validation. It uses genuine survey weights, PSUs, and strata from the package
example data, with strongly skewed welfare-attitude items. Parameter estimates
are still in the same range, and thresholds are almost identical.

The global fit measures are more nuanced. The Mplus WLSMV CFI and TLI align very
closely with lavaan's `scaled` CFI and TLI, while they differ from lavaan's
additional `robust` CFI and TLI. RMSEA is also much closer to lavaan's scaled
RMSEA than to lavaan's robust RMSEA, although the Mplus value remains somewhat
higher. The adjusted chi-square itself is more sensitive to implementation
details such as weight scaling, the baseline-model correction, and the
second-order WLSMV test statistic correction.

This suggests that the ordinal implementation is reproducing the core
threshold/polychoric part of the Mplus workflow well. The earlier-looking
fit-index discrepancy is mostly an apples-to-oranges comparison between Mplus
WLSMV and lavaan's extra robust fit-index column. More real-data comparisons are
still useful before treating the ordinal extension as production-grade.

## Validation 5: Ordinal Multiple-Imputation Survey CFA

### Scope

This validation checks the initial multiple-imputation path in
`lavaan.survey.ordinal()`. The script simulates the same two-factor,
six-indicator CFA structure used in the first ordinal validation, but stores the
observed indicators as four-category ordered variables. The dataset has 600
observations, 60 clusters, 6 strata, and sampling weights.

Artificial missingness is introduced in three indicators:

```text
y2  MAR, depending on y1 and stratum
y5  MAR, depending on y4 and weight
y6  MCAR, about 15 percent
```

The missing ordinal responses are imputed ten times with `mice` using
proportional-odds imputation (`polr`). The same ten completed datasets are then
analyzed by `lavaan.survey.ordinal()` and Mplus Demo.

The model is:

```text
f1 =~ y1 + y2 + y3
f2 =~ y4 + y5 + y6
f1 ~~ f2
```

The corresponding Mplus input uses:

```text
DATA: TYPE = IMPUTATION;
VARIABLE: CATEGORICAL ARE y1-y6;
VARIABLE: WEIGHT IS wgt; CLUSTER IS clu; STRATIFICATION IS str;
ANALYSIS: TYPE = COMPLEX; ESTIMATOR = WLSMV; PARAMETERIZATION = DELTA;
```

### Commands

Prepare the imputed datasets, Mplus input, and `lavaan.survey.ordinal()`
results:

```r
source("validation/mplus-demo/prepare_ordinal_mi_validation_files.R")
```

Run Mplus Demo from `validation/mplus-demo`:

```sh
/Applications/MplusDemo/mpdemo ordinal_mi_complex.inp
```

Parse and compare the output:

```r
source("validation/mplus-demo/compare_mplus_ordinal_mi_output.R")
```

### Fit Measures

| Measure | lavaan scaled | lavaan robust | Mplus WLSMV imputation mean |
| --- | ---: | ---: | ---: |
| Chi-square | 7.045 | -- | 10.300 |
| df | 8 | -- | 8 |
| p-value | 0.532 | -- | -- |
| CFI | 1.000 | 1.000 | 0.999 |
| TLI | 1.001 | 1.015 | 0.998 |
| RMSEA | 0.000 | 0.000 | 0.020 |
| SRMR | 0.0220 | 0.0220 | 0.017 |

As in the continuous MI validation, Mplus reports means over imputed-data fit
statistics, while `lavaan.survey.ordinal()` pools the ordinal sample statistics
and their design-based covariance matrix before refitting one lavaan model.

### Parameter Agreement

The comparison matched 27 unstandardized parameters between the two outputs.

| Quantity | Value |
| --- | ---: |
| Matched parameters | 27 |
| Maximum absolute estimate difference | 0.0073 |
| Maximum absolute standard error difference | 0.0032 |

Largest estimate differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `f2 =~ y5` | 0.9117 | 0.9190 | 0.0073 |
| `f2 =~ y6` | 0.7129 | 0.7180 | 0.0051 |
| `f1 =~ y3` | 0.8677 | 0.8640 | -0.0037 |
| `f2 ~~ f2` | 0.7135 | 0.7110 | -0.0025 |
| `f1 ~~ f2` | 0.2473 | 0.2450 | -0.0023 |
| `f1 ~~ f1` | 0.6962 | 0.6980 | 0.0018 |

Largest standard-error differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `f1 =~ y3` | 0.0562 | 0.0530 | -0.0032 |
| `f2 =~ y6` | 0.0701 | 0.0680 | -0.0021 |
| `f1 =~ y2` | 0.0592 | 0.0610 | 0.0018 |
| `f1 ~~ f1` | 0.0463 | 0.0470 | 0.0007 |
| `y5 | t2` | 0.0555 | 0.0550 | -0.0005 |

Threshold estimates are especially close, with most differences far below
0.001. The largest differences are in a few factor loadings, which is also
where the ESS4 real-data validation was most sensitive.

### Interpretation

This validation is encouraging for the new ordinal MI implementation. It checks
the full path from imputed ordered indicators through complex survey correction
and WLSMV-style ordinal SEM in both programs. Parameter estimates and standard
errors agree closely, especially for thresholds.

The remaining fit-measure differences are expected to be less diagnostic than
the parameter comparison because Mplus and `lavaan.survey.ordinal()` summarize
multiple-imputation fit differently.

## Validation 6: Ordinal Multiple-Group / Invariance Survey CFA

This validation checks `lavaan.survey.ordinal()` for multiple-group ordinal CFA
with survey weights, clusters, and strata. It uses a demo-compatible one-factor
model with four four-category ordered indicators and two groups.

The model is:

```text
f =~ y1 + y2 + y3 + y4
```

Mplus applies the model statement across groups. To match that parameterization
in lavaan, the validation fits loading, threshold, and intercept invariance:

```r
group.equal = c("loadings", "thresholds", "intercepts")
```

This gives the same model degrees of freedom in Mplus and lavaan.

The corresponding Mplus input uses:

```text
VARIABLE: CATEGORICAL ARE y1-y4;
VARIABLE: GROUPING IS grp (1 = boys 2 = girls);
VARIABLE: WEIGHT IS wgt; CLUSTER IS clu; STRATIFICATION IS str;
ANALYSIS: TYPE = COMPLEX; ESTIMATOR = WLSMV; PARAMETERIZATION = DELTA;
MODEL: f BY y1 y2 y3 y4;
```

### Commands

Prepare the dataset, Mplus input, and `lavaan.survey.ordinal()` results:

```r
source("validation/mplus-demo/prepare_ordinal_group_validation_files.R")
```

Run Mplus Demo from `validation/mplus-demo`:

```sh
/Applications/MplusDemo/mpdemo ordinal_group_complex.inp
```

Parse and compare the output:

```r
source("validation/mplus-demo/compare_mplus_ordinal_group_output.R")
```

### Fit Measures

| Measure | lavaan scaled | lavaan robust | Mplus WLSMV |
| --- | ---: | ---: | ---: |
| Chi-square | 10.948 | -- | 11.323 |
| df | 14 | -- | 14 |
| p-value | 0.690 | -- | 0.661 |
| CFI | 1.000 | 1.000 | 1.000 |
| TLI | 1.002 | 1.038 | 1.000 |
| RMSEA | 0.000 | 0.000 | 0.000 |
| SRMR | 0.0172 | 0.0172 | 0.0170 |

### Parameter Agreement

The comparison matched 34 unstandardized parameters between the two outputs.

| Quantity | Value |
| --- | ---: |
| Matched parameters | 34 |
| Maximum absolute estimate difference | 0.0050 |
| Maximum absolute standard error difference | 0.0058 |

Largest estimate differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `girls: f ~~ f` | 0.7340 | 0.7390 | 0.0050 |
| `boys: f =~ y2` | 0.8691 | 0.8650 | -0.0041 |
| `girls: f =~ y2` | 0.8691 | 0.8650 | -0.0041 |
| `boys: f =~ y3` | 0.8538 | 0.8510 | -0.0028 |
| `girls: f =~ y3` | 0.8538 | 0.8510 | -0.0028 |

Largest standard-error differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `boys: y4 | t3` | 0.0788 | 0.0730 | -0.0058 |
| `girls: y4 | t3` | 0.0788 | 0.0730 | -0.0058 |
| `boys: f =~ y4` | 0.0590 | 0.0540 | -0.0050 |
| `girls: f =~ y4` | 0.0590 | 0.0540 | -0.0050 |
| `boys: f =~ y3` | 0.0601 | 0.0560 | -0.0041 |

### Interpretation

This is a strong cross-software check for the multiple-group ordinal path. With
the same invariance parameterization, the global fit measures are very close and
all compared loading, threshold, and factor-variance parameters match closely.

## Validation 7: Ordinal Multiple-Group / Invariance Multiple-Imputation Survey CFA

This validation combines the previous two extensions: grouped ordinal indicators
and multiple imputation. It uses the same one-factor, two-group ordinal CFA
structure as the multiple-group validation and introduces artificial missingness
in two ordered indicators.

Artificial missingness is introduced in:

```text
y2  MAR, depending on y1 and group
y4  MAR, depending on y3 and weight
```

The observed missing proportions are:

| Variable | Missing proportion |
| --- | ---: |
| y1 | 0.000 |
| y2 | 0.255 |
| y3 | 0.000 |
| y4 | 0.233 |

The missing ordered responses are imputed ten times with `mice` using
proportional-odds imputation (`polr`). The same ten completed datasets are then
analyzed by `lavaan.survey.ordinal()` and Mplus Demo.

The lavaan model uses the same invariance parameterization as the Mplus
multiple-group model:

```r
group.equal = c("loadings", "thresholds", "intercepts")
```

The corresponding Mplus input uses:

```text
DATA: TYPE = IMPUTATION;
VARIABLE: CATEGORICAL ARE y1-y4;
VARIABLE: GROUPING IS grp (1 = boys 2 = girls);
VARIABLE: WEIGHT IS wgt; CLUSTER IS clu; STRATIFICATION IS str;
ANALYSIS: TYPE = COMPLEX; ESTIMATOR = WLSMV; PARAMETERIZATION = DELTA;
MODEL: f BY y1 y2 y3 y4;
```

### Commands

Prepare the imputed datasets, Mplus input, and `lavaan.survey.ordinal()`
results:

```r
source("validation/mplus-demo/prepare_ordinal_group_mi_validation_files.R")
```

Run Mplus Demo from `validation/mplus-demo`:

```sh
/Applications/MplusDemo/mpdemo ordinal_group_mi_complex.inp
```

Parse and compare the output:

```r
source("validation/mplus-demo/compare_mplus_ordinal_group_mi_output.R")
```

### Fit Measures

| Measure | lavaan scaled | lavaan robust | Mplus WLSMV imputation mean |
| --- | ---: | ---: | ---: |
| Chi-square | 8.732 | -- | 12.479 |
| df | 14 | -- | 14 |
| p-value | 0.848 | -- | -- |
| CFI | 1.000 | 1.000 | 1.000 |
| TLI | 1.005 | 1.033 | 1.001 |
| RMSEA | 0.000 | 0.076 | 0.006 |
| SRMR | 0.0217 | 0.0217 | 0.019 |

As in the other MI validations, fit measures should be interpreted with care:
Mplus reports summaries across imputed-data analyses, while
`lavaan.survey.ordinal()` pools thresholds, polychoric correlations, and their
design-based covariance matrix before fitting one lavaan model.

### Parameter Agreement

The comparison matched 34 unstandardized parameters between the two outputs.

| Quantity | Value |
| --- | ---: |
| Matched parameters | 34 |
| Maximum absolute estimate difference | 0.0123 |
| Maximum absolute standard error difference | 0.0043 |

Largest estimate differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `girls: f ~~ f` | 0.8923 | 0.8800 | -0.0123 |
| `boys: f =~ y4` | 0.7348 | 0.7440 | 0.0092 |
| `girls: f =~ y4` | 0.7348 | 0.7440 | 0.0092 |
| `boys: y1 | t3` | 1.1638 | 1.1570 | -0.0068 |
| `girls: y1 | t3` | 1.1638 | 1.1570 | -0.0068 |

Largest standard-error differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `boys: y4 | t3` | 0.0743 | 0.0700 | -0.0043 |
| `girls: y4 | t3` | 0.0743 | 0.0700 | -0.0043 |
| `boys: f =~ y4` | 0.0633 | 0.0600 | -0.0033 |
| `girls: f =~ y4` | 0.0633 | 0.0600 | -0.0033 |
| `girls: f ~~ f` | 0.1550 | 0.1580 | 0.0030 |

### Interpretation

This is the broadest Mplus Demo validation currently in the repository. It
checks ordered indicators, grouping, invariance constraints, complex survey
correction, and multiple imputation together. The parameter agreement is strong
given that the two programs summarize MI fit differently.

## Validation 8: Mixed Ordinal/Continuous Multiple-Group MI Survey CFA

This diagnostic extends the mixed-indicator proof of concept to the hardest
case currently supported by the package: two ordered indicators, two continuous
indicators, two groups, loading/threshold/intercept invariance, artificial
missingness, ten imputations, and a complex survey design.

Artificial missingness is introduced in:

```text
y2  MAR, depending on y1 and group
x2  MAR, depending on x1 and weight
```

The observed missing proportions are:

| Variable | Missing proportion |
| --- | ---: |
| y1 | 0.000 |
| y2 | 0.255 |
| x1 | 0.000 |
| x2 | 0.194 |

The ordered response is imputed with proportional-odds imputation (`polr`), and
the continuous response is imputed with predictive mean matching (`pmm`). The
same ten completed datasets are then analyzed by `lavaan.survey.ordinal()` and
Mplus Demo.

### Commands

Prepare the imputed datasets, Mplus input, and `lavaan.survey.ordinal()`
results:

```r
source("validation/mplus-demo/prepare_mixed_group_mi_validation_files.R")
```

Run Mplus Demo from `validation/mplus-demo`:

```sh
/Applications/MplusDemo/mpdemo mixed_group_mi_complex.inp
```

Parse and compare the output:

```r
source("validation/mplus-demo/compare_mplus_mixed_group_mi_output.R")
```

### Mplus Status

Mplus Demo successfully fits the model with:

| Quantity | Value |
| --- | ---: |
| Requested imputations | 10 |
| Completed imputations | 10 |
| Dependent variables | 4 |
| Groups | 2 |
| Free parameters | 18 |

This means the Mplus Demo limit is not the blocker for the mixed
multiple-group MI validation.

### Fit Measures

The mixed MI workflow now writes both lavaan algorithms:

- `parameter_pooling`: the Mplus-nearer mixed-MI default under the `auto`
  settings, fitting each imputation with lavaan sampling weights and pooling
  parameters with Rubin's rules.
- `sample_statistics`: the original sensitivity path, pooling WLS sample
  statistics and their design-based covariance matrix before one refit.

| Measure | lavaan sample-stat scaled | lavaan parameter-pooling mean | Mplus WLSMV imputation mean |
| --- | ---: | ---: | ---: |
| Chi-square | 88.606 | 13.273 | 12.848 |
| df | 10 | 10 | 10 |
| p-value | 0.000 | -- | -- |
| CFI | 0.909 | 0.9958 | 0.997 |
| TLI | 0.891 | -- | 0.996 |
| RMSEA | 0.165 | 0.0317 | 0.030 |
| SRMR | 0.0271 | 0.0275 | 0.026 |

### Parameter Agreement

The comparison matched all 22 selected unstandardized loading, threshold,
intercept, and factor-variance parameters. Agreement is not adequate for this
to count as a passed external validation for the default sample-statistic
algorithm.

| Quantity | Value |
| --- | ---: |
| Matched parameters | 22 |
| Maximum absolute estimate difference | 0.4790 |
| Maximum absolute standard error difference | 0.0350 |

Largest estimate differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `boys: x1 ~ 1` | -0.3390 | 0.1400 | 0.4790 |
| `girls: x1 ~ 1` | -0.3390 | 0.1400 | 0.4790 |
| `boys: x2 ~ 1` | 0.1192 | -0.2680 | -0.3872 |
| `girls: x2 ~ 1` | 0.1192 | -0.2680 | -0.3872 |
| `girls: y1 | t2` | 0.0004 | 0.2790 | 0.2786 |

Largest standard-error differences:

| Parameter | lavaan.survey.ordinal | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `girls: f ~~ f` | 0.0960 | 0.1310 | 0.0350 |
| `boys: f =~ y2` | 0.0677 | 0.0590 | -0.0087 |
| `girls: f =~ y2` | 0.0677 | 0.0590 | -0.0087 |
| `boys: f =~ x1` | 0.0706 | 0.0620 | -0.0086 |
| `girls: f =~ x1` | 0.0706 | 0.0620 | -0.0086 |

The experimental parameter-pooling path is much closer to Mplus:

| Quantity | Value |
| --- | ---: |
| Matched free parameters | 20 |
| Maximum absolute estimate difference | 0.0221 |
| Maximum absolute standard error difference | 0.0109 |

Largest estimate differences for parameter pooling:

| Parameter | lavaan parameter pooling | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `boys: f =~ x1` | 0.9341 | 0.9120 | -0.0221 |
| `girls: f =~ x1` | 0.9341 | 0.9120 | -0.0221 |
| `boys: y1 | t2` | 0.3008 | 0.2790 | -0.0218 |
| `girls: y1 | t2` | 0.3008 | 0.2790 | -0.0218 |
| `boys: y1 | t1` | -0.6439 | -0.6280 | 0.0159 |

Largest standard-error differences for parameter pooling:

| Parameter | lavaan parameter pooling | Mplus Demo | Mplus - lavaan |
| --- | ---: | ---: | ---: |
| `girls: f ~~ f` | 0.1419 | 0.1310 | -0.0109 |
| `boys: y2 | t2` | 0.0778 | 0.0670 | -0.0108 |
| `girls: y2 | t2` | 0.0778 | 0.0670 | -0.0108 |
| `boys: f ~~ f` | 0.0601 | 0.0510 | -0.0091 |
| `boys: y2 | t1` | 0.0778 | 0.0860 | 0.0082 |

### Interpretation

This is now a useful fork in the evidence rather than a simple failure. Mplus
fits the intended mixed, grouped, imputed, complex-survey model and reports
pooled imputation-results. The default `lavaan.survey.ordinal()` algorithm
pools the WLS sample statistics and their design-based covariance matrix first,
then fits one lavaan model. For all-ordinal MI examples this approximation
matched Mplus closely, but in this mixed multiple-group case the distinction is
substantive.

The experimental parameter-pooling path is much closer to Mplus for both
parameter estimates and fit-measure means. It should still be treated as
experimental because it currently mirrors Mplus's imputation layer more closely
than its full complex-survey correction layer, but it gives us a practical,
testable path for mixed ordinal/continuous MI models.

## Files Produced Locally

The scripts produce these generated files, which are intentionally ignored by
Git:

- `ordinal_survey_complex.dat`
- `ordinal_survey_complex.inp`
- `ordinal_survey_complex.out`
- `lavaan_survey_complex_parameters.csv`
- `lavaan_survey_complex_fit.csv`
- `mplus_lavaan_fit_comparison.csv`
- `mplus_complex_parameters_raw.csv`
- `mplus_complex_fit_summary.csv`
- `mplus_lavaan_parameter_comparison.csv`
- `ess4_range_complex.dat`
- `ess4_range_complex.inp`
- `ess4_range_complex.out`
- `ess4_collapsed_category_counts.csv`
- `ess4_lavaan_survey_parameters.csv`
- `ess4_lavaan_survey_fit.csv`
- `ess4_mplus_lavaan_fit_comparison.csv`
- `ess4_mplus_parameters_raw.csv`
- `ess4_mplus_fit_summary.csv`
- `ess4_mplus_lavaan_parameter_comparison.csv`
- `continuous_mi_imp01.dat` through `continuous_mi_imp10.dat`
- `continuous_mi_implist.dat`
- `continuous_mi_complex.inp`
- `continuous_mi_complex.out`
- `continuous_mi_missing_summary.csv`
- `continuous_mi_lavaan_survey_parameters.csv`
- `continuous_mi_lavaan_survey_fit.csv`
- `continuous_mi_mplus_fit_summary.csv`
- `continuous_mi_mplus_parameters_raw.csv`
- `continuous_mi_mplus_lavaan_fit_comparison.csv`
- `continuous_mi_mplus_lavaan_parameter_comparison.csv`
- `ordinal_mi_imp01.dat` through `ordinal_mi_imp10.dat`
- `ordinal_mi_implist.dat`
- `ordinal_mi_complex.inp`
- `ordinal_mi_complex.out`
- `ordinal_mi_missing_summary.csv`
- `ordinal_mi_lavaan_survey_parameters.csv`
- `ordinal_mi_lavaan_survey_fit.csv`
- `ordinal_mi_mplus_fit_summary.csv`
- `ordinal_mi_mplus_parameters_raw.csv`
- `ordinal_mi_mplus_lavaan_fit_comparison.csv`
- `ordinal_mi_mplus_lavaan_parameter_comparison.csv`
- `ordinal_group_complex.dat`
- `ordinal_group_complex.inp`
- `ordinal_group_complex.out`
- `ordinal_group_lavaan_survey_parameters.csv`
- `ordinal_group_lavaan_survey_fit.csv`
- `ordinal_group_mplus_fit_summary.csv`
- `ordinal_group_mplus_parameters_raw.csv`
- `ordinal_group_mplus_lavaan_fit_comparison.csv`
- `ordinal_group_mplus_lavaan_parameter_comparison.csv`
- `ordinal_group_mi_imp01.dat` through `ordinal_group_mi_imp10.dat`
- `ordinal_group_mi_implist.dat`
- `ordinal_group_mi_complex.inp`
- `ordinal_group_mi_complex.out`
- `ordinal_group_mi_missing_summary.csv`
- `ordinal_group_mi_lavaan_survey_parameters.csv`
- `ordinal_group_mi_lavaan_survey_fit.csv`
- `ordinal_group_mi_mplus_fit_summary.csv`
- `ordinal_group_mi_mplus_parameters_raw.csv`
- `ordinal_group_mi_mplus_lavaan_fit_comparison.csv`
- `ordinal_group_mi_mplus_lavaan_parameter_comparison.csv`
- `mixed_group_mi_imp01.dat` through `mixed_group_mi_imp10.dat`
- `mixed_group_mi_implist.dat`
- `mixed_group_mi_complex.inp`
- `mixed_group_mi_complex.out`
- `mixed_group_mi_missing_summary.csv`
- `mixed_group_mi_lavaan_survey_parameters.csv`
- `mixed_group_mi_lavaan_survey_fit.csv`
- `mixed_group_mi_lavaan_survey_parameters_parameter_pooling.csv`
- `mixed_group_mi_lavaan_survey_fit_parameter_pooling.csv`
- `mixed_group_mi_mplus_fit_summary.csv`
- `mixed_group_mi_mplus_parameters_raw.csv`
- `mixed_group_mi_mplus_lavaan_fit_comparison.csv`
- `mixed_group_mi_mplus_lavaan_parameter_comparison.csv`
- `mixed_group_mi_mplus_lavaan_parameter_comparison_parameter_pooling.csv`
- `mixed_group_mi_mplus_lavaan_parameter_comparison_all_algorithms.csv`

## References

- Mplus User's Guide examples show WLSMV output with a single fit block and the
  WLSMV chi-square difference-testing note:
  https://www.statmodel.com/usersguide/chap3/ex3.14.html
- Asparouhov and Muthen's note on the WLSMV/ULSMV/MLMV second-order
  chi-square correction:
  https://www.statmodel.com/download/WLSMV_new_chi21.pdf
