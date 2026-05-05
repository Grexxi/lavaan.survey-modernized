# Mplus Demo Validation

This folder contains cross-software validation workflows for the modernized
`lavaan.survey` package. The validation ladder is ordered from the original
continuous workflow to the new ordinal workflow:

| Step | Target | Status | Prepare script | Mplus input | Compare script |
| --- | --- | --- | --- | --- | --- |
| 1 | Continuous survey SEM | Covered by package tests and the Roosma example; no separate Mplus Demo workflow yet | -- | -- | -- |
| 2 | Continuous survey SEM with multiple imputation | Implemented | `prepare_continuous_mi_validation_files.R` | `continuous_mi_complex.inp` | `compare_mplus_continuous_mi_output.R` |
| 3 | Continuous multiple-group / invariance models | Planned next validation target | -- | -- | -- |
| 4a | Ordinal survey SEM, simulated data | Implemented | `prepare_validation_files.R` | `ordinal_survey_complex.inp` | `compare_mplus_output.R` |
| 4b | Ordinal survey SEM, ESS4 GB data | Implemented | `prepare_ess4_validation_files.R` | `ess4_range_complex.inp` | `compare_mplus_ess4_output.R` |
| 5 | Ordinal survey SEM with multiple imputation | Implemented | `prepare_ordinal_mi_validation_files.R` | `ordinal_mi_complex.inp` | `compare_mplus_ordinal_mi_output.R` |
| 6 | Ordinal multiple-group / invariance models | Implemented | `prepare_ordinal_group_validation_files.R` | `ordinal_group_complex.inp` | `compare_mplus_ordinal_group_output.R` |
| 7 | Ordinal multiple-group / invariance models with multiple imputation | Implemented | `prepare_ordinal_group_mi_validation_files.R` | `ordinal_group_mi_complex.inp` | `compare_mplus_ordinal_group_mi_output.R` |
| 8 | Mixed ordinal/continuous multiple-group models with multiple imputation | Diagnostic; compares the default pooled-statistic algorithm with experimental Rubin parameter pooling, which is much closer to Mplus `TYPE = IMPUTATION` | `prepare_mixed_group_mi_validation_files.R` | `mixed_group_mi_complex.inp` | `compare_mplus_mixed_group_mi_output.R` |

The Mplus Demo version is limited to six dependent variables. The simulated
ordinal and continuous validation models therefore use compact CFA structures
that fit inside those limits. The ESS4 GB workflow uses six collapsed ordinal
items from the bundled `ess4.gb` dataset.

The folder also includes `mplus_ex510.inp` and `ex5.10.dat`, downloaded from the
official Mplus User Guide example 5.10. This is a plain Mplus sanity check for
categorical CFA and is not a direct `lavaan.survey` comparison.

## 1. Continuous Survey SEM

The basic continuous path is covered by package tests rather than a dedicated
Mplus Demo workflow:

- `tests/testthat/test_roosma.R` checks the bundled ESS4/Roosma example against
  legacy numerical targets.
- `tests/testthat/test_continuous_robustness.R` checks modern robustness cases
  such as Yuan-Bentler without a mean structure, multiple-group refits, and
  `pval.pFsum()`.

## 2. Continuous Survey SEM With Multiple Imputation

From the package root, prepare the imputed datasets, Mplus input, and
`lavaan.survey()` results:

```r
source("validation/mplus-demo/prepare_continuous_mi_validation_files.R")
```

Run Mplus Demo from this folder:

```sh
/Applications/MplusDemo/mpdemo continuous_mi_complex.inp
```

Compare outputs from the package root:

```r
source("validation/mplus-demo/compare_mplus_continuous_mi_output.R")
```

## 3. Continuous Multiple-Group / Invariance Models

This is the next missing validation cell. A future workflow should mirror the
ordinal invariance checks with continuous indicators, survey weights, clusters,
strata, and equality constraints such as loadings and intercepts.

## 4. Ordinal Survey SEM

### Simulated Ordinal Survey CFA

Prepare files and `lavaan.survey.ordinal()` results:

```r
source("validation/mplus-demo/prepare_validation_files.R")
```

Run Mplus Demo:

```sh
/Applications/MplusDemo/mpdemo ordinal_survey_complex.inp
```

Compare outputs:

```r
source("validation/mplus-demo/compare_mplus_output.R")
```

### ESS4 GB Ordinal Survey CFA

Prepare files and `lavaan.survey.ordinal()` results:

```r
source("validation/mplus-demo/prepare_ess4_validation_files.R")
```

Run Mplus Demo:

```sh
/Applications/MplusDemo/mpdemo ess4_range_complex.inp
```

Compare outputs:

```r
source("validation/mplus-demo/compare_mplus_ess4_output.R")
```

## 5. Ordinal Survey SEM With Multiple Imputation

Prepare imputed datasets, Mplus input, and `lavaan.survey.ordinal()` results:

```r
source("validation/mplus-demo/prepare_ordinal_mi_validation_files.R")
```

Run Mplus Demo:

```sh
/Applications/MplusDemo/mpdemo ordinal_mi_complex.inp
```

Compare outputs:

```r
source("validation/mplus-demo/compare_mplus_ordinal_mi_output.R")
```

## 6. Ordinal Multiple-Group / Invariance Models

Prepare the grouped ordinal dataset, Mplus input, and
`lavaan.survey.ordinal()` results:

```r
source("validation/mplus-demo/prepare_ordinal_group_validation_files.R")
```

Run Mplus Demo:

```sh
/Applications/MplusDemo/mpdemo ordinal_group_complex.inp
```

Compare outputs:

```r
source("validation/mplus-demo/compare_mplus_ordinal_group_output.R")
```

## 7. Ordinal Multiple-Group / Invariance Models With Multiple Imputation

Prepare imputed grouped ordinal datasets, Mplus input, and
`lavaan.survey.ordinal()` results:

```r
source("validation/mplus-demo/prepare_ordinal_group_mi_validation_files.R")
```

Run Mplus Demo:

```sh
/Applications/MplusDemo/mpdemo ordinal_group_mi_complex.inp
```

Compare outputs:

```r
source("validation/mplus-demo/compare_mplus_ordinal_group_mi_output.R")
```

## 8. Mixed Ordinal/Continuous Multiple-Group Models With Multiple Imputation

Prepare imputed grouped mixed-indicator datasets, Mplus input, and
`lavaan.survey.ordinal()` results:

```r
source("validation/mplus-demo/prepare_mixed_group_mi_validation_files.R")
```

Run Mplus Demo:

```sh
/Applications/MplusDemo/mpdemo mixed_group_mi_complex.inp
```

Compare outputs:

```r
source("validation/mplus-demo/compare_mplus_mixed_group_mi_output.R")
```

This workflow is intentionally diagnostic. It confirms that Mplus Demo can fit
the mixed, grouped, imputed, complex-survey CFA model and compares two
`lavaan.survey.ordinal()` algorithms. For all-ordinal MI models, the default
path pools sample statistics and their design-based covariance matrix before
refitting one model. For mixed ordinal/continuous MI models, the default
`auto` settings use the Mplus-nearer parameter-pooling path: each imputed
dataset is fitted with lavaan sampling weights, and model parameters are pooled
with Rubin's rules. The original sample-statistic algorithm remains available
with `point.wls = "design", mi.pooling = "sample.statistics"`.

## Interpreting Comparisons

Do not expect exact equality. The main target is close agreement in loadings,
thresholds, intercepts, factor variances/covariances, and standard errors.

For WLSMV, compare the Mplus fit block primarily with lavaan's scaled fit
measures. lavaan's additional robust CFI, TLI, and RMSEA are useful sensitivity
checks, but Mplus does not print a separate lavaan-style robust fit-index column
for WLSMV.

For multiple imputation, Mplus prints means and standard deviations over the
imputed-data analyses for several fit measures. The mixed-indicator default now
mirrors Mplus's imputation layer more closely through Rubin parameter pooling,
while the sample-statistic algorithm remains a sensitivity check. Parameter
estimates should therefore be the primary validation target for MI workflows,
with fit measures treated as a useful diagnostic summary.
