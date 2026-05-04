# Mplus Demo Validation Results

This note records a first cross-software validation run for
`lavaan.survey.ordinal()` against Mplus Demo 9.

## Scope

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

## Commands

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

## Fit Measures

| Measure | lavaan.survey.ordinal | Mplus Demo |
| --- | ---: | ---: |
| Scaled chi-square | 8.481 | 7.925 |
| df | 8 | 8 |
| p-value | 0.388 | 0.441 |
| CFI | 1.000 | 1.000 |
| TLI | 1.004 | 1.000 |
| RMSEA | 0.000 | 0.000 |
| SRMR | 0.0206 | 0.015 |

## Parameter Agreement

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

## Interpretation

This first validation run is encouraging. The core parameter estimates,
thresholds, and standard errors from `lavaan.survey.ordinal()` closely reproduce
the Mplus Demo results for a demo-compatible ordinal survey CFA. The fit
statistics are also in the same range.

Exact equality should not be expected. Mplus and `lavaan` differ in internal
rounding, weight scaling, robust test statistic corrections, and output
precision. The goal of this validation is therefore numerical agreement within a
small tolerance, not bit-for-bit identity.

## Files Produced Locally

The scripts produce these generated files, which are intentionally ignored by
Git:

- `ordinal_survey_complex.dat`
- `ordinal_survey_complex.inp`
- `ordinal_survey_complex.out`
- `lavaan_survey_complex_parameters.csv`
- `lavaan_survey_complex_fit.csv`
- `mplus_complex_parameters_raw.csv`
- `mplus_complex_fit_summary.csv`
- `mplus_lavaan_parameter_comparison.csv`

