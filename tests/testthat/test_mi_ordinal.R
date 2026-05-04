suppressPackageStartupMessages(suppressWarnings(library(lavaan.survey)))

context("ordinal multiple-imputation survey models")

make_ordinal_mi_data <- function() {
  set.seed(2402)
  n <- 240
  cluster <- rep(seq_len(40), each=6)
  weight <- round(stats::runif(n, 0.7, 1.5), 4)
  factor_score <- stats::rnorm(n)

  to_ordered <- function(x) {
    ordered(cut(x, breaks=c(-Inf, -0.5, 0.5, Inf), labels=1:3),
            levels=1:3)
  }

  dat <- data.frame(
    y1=to_ordered(0.85 * factor_score + stats::rnorm(n, sd=0.55)),
    y2=to_ordered(0.75 * factor_score + stats::rnorm(n, sd=0.65)),
    y3=to_ordered(0.80 * factor_score + stats::rnorm(n, sd=0.60)),
    y4=to_ordered(0.70 * factor_score + stats::rnorm(n, sd=0.70)),
    cluster=cluster,
    weight=weight
  )

  missing_y3 <- sample(seq_len(n), 48)
  dat$y3[missing_y3] <- NA

  lapply(seq_len(3), function(i) {
    imp <- dat
    set.seed(9100 + i)
    observed <- imp$y3[!is.na(imp$y3)]
    imp$y3[is.na(imp$y3)] <- sample(observed, sum(is.na(imp$y3)),
                                    replace=TRUE)
    imp$y3 <- ordered(imp$y3, levels=levels(dat$y3))
    imp
  })
}

test_that("ordinal MI survey models pool thresholds and correlations", {
  skip_if_not_installed("mitools")

  items <- paste0("y", 1:4)
  imputed_data <- make_ordinal_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    weights=~weight,
    data=imputation_list
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=8
  )

  fit_naive <- lavaan::cfa(
    model="f =~ y1 + y2 + y3 + y4",
    data=imputed_data[[1]],
    ordered=items
  )

  fit_mi <- suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit_naive,
    survey.design=rep_design,
    estimator="WLSMV"
  ))

  per_imputation <- lapply(
    rep_design$designs,
    lavaan.survey:::get.ordinal.survey.stats,
    ov.names=items,
    ordered=items,
    ngroups=1
  )
  expected <- lavaan.survey:::pool.ordinal.mi.stats(
    per_imputation,
    ngroups=1
  )

  sampstat <- lavaan::lavInspect(fit_mi, "sampstat")
  Gamma <- lavaan::lavTech(fit_mi, "gamma")[[1]]

  expect_true(lavaan::lavInspect(fit_mi, "converged"))
  expect_equal(unclass(as.matrix(sampstat$cov)),
               unclass(as.matrix(expected$point.stats$sample.cov)),
               tolerance=1e-10, check.attributes=FALSE)
  expect_equal(as.numeric(sampstat$th),
               as.numeric(expected$point.stats$sample.th),
               tolerance=1e-10, check.attributes=FALSE)
  expect_equal(Gamma, expected$Gamma,
               tolerance=1e-10, check.attributes=FALSE)
  expect_true(is.finite(lavaan::fitMeasures(fit_mi, "chisq.scaled")))
})
