suppressPackageStartupMessages(suppressWarnings(library(lavaan.survey)))

make_mixed_indicator_data <- function(n=300) {
  set.seed(2601)
  factor_score <- stats::rnorm(n)
  y1star <- 0.90 * factor_score + stats::rnorm(n, sd=0.55)
  y2star <- 0.75 * factor_score + stats::rnorm(n, sd=0.65)

  data.frame(
    y1=ordered(cut(y1star, breaks=c(-Inf, -0.6, 0.2, Inf), labels=1:3)),
    y2=ordered(cut(y2star, breaks=c(-Inf, -0.4, 0.5, Inf), labels=1:3)),
    x1=0.80 * factor_score + stats::rnorm(n, sd=0.60),
    x2=0.70 * factor_score + stats::rnorm(n, sd=0.70),
    cluster=factor(rep(seq_len(30), each=n / 30)),
    stratum=factor(rep(rep(seq_len(3), each=10), length.out=n)),
    weight=round(stats::runif(n, 0.6, 1.8), 4)
  )
}

fit_mixed_indicator_model <- function(data) {
  lavaan::cfa(
    model="f =~ y1 + y2 + x1 + x2",
    data=data,
    ordered=c("y1", "y2"),
    estimator="WLSMV"
  )
}

test_that("mixed ordinal and continuous survey CFA works", {
  data <- make_mixed_indicator_data()
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=data,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=12
  )
  fit <- fit_mixed_indicator_model(data)

  fit_survey <- suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit,
    survey.design=rep_design,
    estimator="WLSMV"
  ))

  pe <- lavaan::parameterEstimates(fit_survey)
  loadings <- pe[pe$op == "=~", ]
  thresholds <- pe[pe$op == "|", ]

  expect_true(lavaan::lavInspect(fit_survey, "converged"))
  expect_true(is.finite(lavaan::fitMeasures(fit_survey, "chisq.scaled")))
  expect_equal(sort(loadings$rhs), c("x1", "x2", "y1", "y2"))
  expect_equal(sort(unique(thresholds$lhs)), c("y1", "y2"))
  expect_true(all(is.finite(pe$se[!is.na(pe$se) & pe$se > 0])))
})

test_that("mixed WLS observations use lavaan's ordering", {
  data <- make_mixed_indicator_data()
  ov.names <- c("y1", "y2", "x1", "x2")
  ordered <- c("y1", "y2")

  full <- lavaan.survey:::get.ordinal.stats(
    data=data,
    ov.names=ov.names,
    ordered=ordered,
    weights=data$weight,
    full=TRUE
  )
  lite <- lavaan.survey:::get.ordinal.stats(
    data=data,
    ov.names=ov.names,
    ordered=ordered,
    weights=data$weight,
    wls.names=names(full$wls.obs)
  )

  expect_equal(names(lite$wls.obs), names(full$wls.obs))
  expect_equal(as.numeric(lite$wls.obs), as.numeric(full$wls.obs),
               tolerance=1e-10)
  expect_true(all(c("y1|t1", "y2|t2", "x1~1", "x2~1", "x1~~x1",
                    "x2~~x2", "y1~~x1", "y2~~x2", "x1~~x2") %in%
                    names(full$wls.obs)))
})

test_that("mixed ordered argument must be a subset of observed variables", {
  data <- make_mixed_indicator_data()
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=data,
    nest=TRUE
  )
  fit <- fit_mixed_indicator_model(data)

  expect_error(
    lavaan.survey.ordinal(
      lavaan.fit=fit,
      survey.design=design,
      ordered=c("y1", "not_in_model"),
      estimator="WLSMV",
      replicates=4
    ),
    "Ordered variables are not observed model variables"
  )
})
