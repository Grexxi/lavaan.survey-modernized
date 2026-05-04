# Mplus Demo Validation

This folder contains a small cross-software validation workflow for the
experimental `lavaan.survey.ordinal()` function.

The Mplus Demo version is limited to six dependent variables, so the validation
model deliberately uses six ordinal indicators and two factors:

```text
f1 =~ y1 + y2 + y3
f2 =~ y4 + y5 + y6
f1 ~~ f2
```

The workflow creates two sets of files:

1. `mplus_ex510.inp` plus `ex5.10.dat`, downloaded from the official Mplus User
   Guide example 5.10. This is a plain Mplus sanity check for categorical CFA.
2. `ordinal_survey_complex.inp` plus `ordinal_survey_complex.dat`, a simulated
   four-category ordinal survey dataset with weights, clusters, and strata.

The second dataset is the actual comparison target for `lavaan.survey.ordinal()`.

The folder also contains an ESS4 GB validation workflow. The `ess4.gb` dataset
is bundled with `lavaan.survey`; its documentation states that it comes from the
European Social Survey round 4 United Kingdom sample, downloaded from the ESS
data portal and converted to an R dataset.

## Prepare files and lavaan.survey results

From the package root, run:

```r
source("validation/mplus-demo/prepare_validation_files.R")
```

For the ESS4 GB validation, run:

```r
source("validation/mplus-demo/prepare_ess4_validation_files.R")
```

This writes:

- `ex5.10.dat`
- `mplus_ex510.inp`
- `ordinal_survey_complex.dat`
- `ordinal_survey_complex.inp`
- `lavaan_survey_complex_parameters.csv`
- `lavaan_survey_complex_fit.csv`

## Run Mplus Demo

If `mpdemo` is on your `PATH`, run from this folder:

```sh
mpdemo ordinal_survey_complex.inp
```

For the ESS4 GB validation:

```sh
mpdemo ess4_range_complex.inp
```

You can also run the official Mplus example sanity check:

```sh
mpdemo mplus_ex510.inp
```

## Compare outputs

After Mplus has created `ordinal_survey_complex.out`, run from the package root:

```r
source("validation/mplus-demo/compare_mplus_output.R")
```

For the ESS4 GB validation:

```r
source("validation/mplus-demo/compare_mplus_ess4_output.R")
```

The comparison script uses `MplusAutomation` if available:

```r
install.packages("MplusAutomation")
```

It writes:

- `mplus_lavaan_parameter_comparison.csv`
- `mplus_complex_parameters_raw.csv`
- `mplus_complex_fit_summary.csv`

Do not expect exact equality. The main target is close agreement in loadings,
thresholds, factor variances/covariances, and the overall pattern of fit. Robust
standard errors and scaled test statistics can differ because the two
implementations use different internal corrections and weight scaling details.
