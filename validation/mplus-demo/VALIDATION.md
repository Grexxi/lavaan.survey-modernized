# Mplus Demo Validation Results

This note records a first cross-software validation run for
`lavaan.survey.ordinal()` against Mplus Demo 9.

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

## Validation 1: Simulated Ordinal Survey CFA

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

## Validation 2: ESS4 GB Ordinal Survey CFA

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

## Validation 3: Continuous Multiple-Imputation Survey CFA

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

## References

- Mplus User's Guide examples show WLSMV output with a single fit block and the
  WLSMV chi-square difference-testing note:
  https://www.statmodel.com/usersguide/chap3/ex3.14.html
- Asparouhov and Muthen's note on the WLSMV/ULSMV/MLMV second-order
  chi-square correction:
  https://www.statmodel.com/download/WLSMV_new_chi21.pdf
