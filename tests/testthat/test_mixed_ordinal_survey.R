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

make_mixed_indicator_group_data <- function(n=400) {
  set.seed(2602)
  group <- factor(rep(c("Girls", "Boys"), each=n / 2),
                  levels=c("Girls", "Boys"))
  factor_score <- stats::rnorm(n) + ifelse(group == "Girls", 0.10, -0.10)
  y1star <- 0.90 * factor_score + stats::rnorm(n, sd=0.55)
  y2star <- 0.78 * factor_score + stats::rnorm(n, sd=0.62)

  data.frame(
    y1=ordered(cut(y1star, breaks=c(-Inf, -0.70, 0.15, Inf), labels=1:3)),
    y2=ordered(cut(y2star, breaks=c(-Inf, -0.45, 0.55, Inf), labels=1:3)),
    x1=0.82 * factor_score + stats::rnorm(n, sd=0.60),
    x2=0.68 * factor_score + stats::rnorm(n, sd=0.72),
    cluster=factor(rep(seq_len(40), each=n / 40)),
    stratum=factor(rep(rep(seq_len(4), each=10), length.out=n)),
    weight=round(stats::runif(n, 0.6, 1.8), 4),
    group=group
  )
}

fit_mixed_indicator_model <- function(data, group=NULL, group.equal=NULL,
                                      meanstructure=NULL) {
  args <- list(
    model="f =~ y1 + y2 + x1 + x2",
    data=data,
    ordered=c("y1", "y2"),
    estimator="WLSMV"
  )
  if(!is.null(group)) args$group <- group
  if(!is.null(group.equal)) args$group.equal <- group.equal
  if(!is.null(meanstructure)) args$meanstructure <- meanstructure

  do.call(lavaan::cfa, args)
}

make_mixed_indicator_design <- function(data) {
  survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=data,
    nest=TRUE
  )
}

make_mixed_indicator_rep_design <- function(data, replicates=12) {
  survey::as.svrepdesign(
    make_mixed_indicator_design(data),
    type="bootstrap",
    replicates=replicates
  )
}

test_that("mixed ordinal and continuous survey CFA works", {
  data <- make_mixed_indicator_data()
  rep_design <- make_mixed_indicator_rep_design(data)
  fit <- fit_mixed_indicator_model(data, meanstructure=FALSE)

  fit_survey <- suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit,
    survey.design=rep_design,
    estimator="WLSMV"
  ))

  pe <- lavaan::parameterEstimates(fit_survey)
  loadings <- pe[pe$op == "=~", ]
  thresholds <- pe[pe$op == "|", ]

  expect_true(lavaan::lavInspect(fit_survey, "converged"))
  expect_true(lavaan::lavInspect(fit_survey, "options")$meanstructure)
  expect_true(is.finite(lavaan::fitMeasures(fit_survey, "chisq.scaled")))
  expect_equal(sort(loadings$rhs), c("x1", "x2", "y1", "y2"))
  expect_equal(sort(unique(thresholds$lhs)), c("y1", "y2"))
  expect_true(all(is.finite(pe$se[!is.na(pe$se) & pe$se > 0])))
})

test_that("mixed ordinal and continuous survey CFA works for multiple groups", {
  data <- make_mixed_indicator_group_data()
  rep_design <- make_mixed_indicator_rep_design(data)
  fit <- fit_mixed_indicator_model(data, group="group")

  fit_survey <- suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit,
    survey.design=rep_design,
    estimator="WLSMV"
  ))

  pe <- lavaan::parameterEstimates(fit_survey)
  loadings <- pe[pe$op == "=~", ]
  thresholds <- pe[pe$op == "|", ]

  expect_true(lavaan::lavInspect(fit_survey, "converged"))
  expect_equal(lavaan::lavInspect(fit_survey, "ngroups"), 2)
  expect_equal(lavaan::lavInspect(fit_survey, "nobs"), c(200L, 200L))
  expect_equal(length(lavaan::lavTech(fit_survey, "gamma")), 2)
  expect_true(is.finite(lavaan::fitMeasures(fit_survey, "chisq.scaled")))
  expect_equal(sort(unique(loadings$group)), c(1L, 2L))
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

test_that("mixed grouped WLS observations use lavaan's ordering by group", {
  data <- make_mixed_indicator_group_data()
  ov.names <- c("y1", "y2", "x1", "x2")
  ordered <- c("y1", "y2")
  group_labels <- c("Girls", "Boys")

  full <- lavaan.survey:::get.ordinal.stats(
    data=data,
    ov.names=ov.names,
    ordered=ordered,
    weights=data$weight,
    group="group",
    group.labels=group_labels,
    full=TRUE
  )
  lite <- lavaan.survey:::get.ordinal.stats(
    data=data,
    ov.names=ov.names,
    ordered=ordered,
    weights=data$weight,
    group="group",
    group.labels=group_labels,
    wls.names=lapply(full$wls.obs, names)
  )

  expect_equal(names(lite$wls.obs), group_labels)
  expect_equal(names(full$wls.obs), group_labels)
  for(g in group_labels) {
    expect_equal(names(lite$wls.obs[[g]]), names(full$wls.obs[[g]]))
    expect_equal(as.numeric(lite$wls.obs[[g]]), as.numeric(full$wls.obs[[g]]),
                 tolerance=1e-10)
    expect_true(all(c("y1|t1", "y2|t2", "x1~1", "x2~1", "x1~~x1",
                      "x2~~x2", "y1~~x1", "y2~~x2", "x1~~x2") %in%
                      names(full$wls.obs[[g]])))
  }
})

test_that("mixed multiple-group models preserve invariance constraints", {
  data <- make_mixed_indicator_group_data()
  design <- make_mixed_indicator_design(data)
  fit <- fit_mixed_indicator_model(
    data,
    group="group",
    group.equal=c("loadings", "thresholds", "intercepts")
  )

  fit_survey <- suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit,
    survey.design=design,
    estimator="WLSMV",
    replicates=12
  ))

  partable <- lavaan::parTable(fit_survey)
  constraints <- partable[partable$op == "==", ]
  loadings <- partable[partable$op == "=~" & partable$rhs != "y1", ]
  thresholds <- partable[partable$op == "|", ]
  intercepts <- partable[partable$op == "~1" & partable$lhs %in% c("x1", "x2"), ]

  expect_true(lavaan::lavInspect(fit_survey, "converged"))
  expect_equal(lavaan::lavInspect(fit_survey, "ngroups"), 2)
  expect_gt(nrow(constraints), 0)
  expect_true(any(loadings$label != ""))
  expect_true(any(thresholds$label != ""))
  expect_true(any(intercepts$label != ""))
  expect_true(is.finite(lavaan::fitMeasures(fit_survey, "chisq.scaled")))
})

test_that("mixed ordered argument must be a subset of observed variables", {
  data <- make_mixed_indicator_data()
  design <- make_mixed_indicator_design(data)
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
