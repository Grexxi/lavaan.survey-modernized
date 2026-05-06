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

make_mixed_indicator_mi_data <- function() {
  data <- make_mixed_indicator_data()
  n <- nrow(data)

  set.seed(2603)
  missing_y2 <- sample(seq_len(n), 60)
  missing_x2 <- sample(setdiff(seq_len(n), missing_y2), 55)
  data$y2[missing_y2] <- NA
  data$x2[missing_x2] <- NA

  lapply(seq_len(3), function(i) {
    imputed <- data
    set.seed(9300 + i)

    observed_y2 <- imputed$y2[!is.na(imputed$y2)]
    imputed$y2[is.na(imputed$y2)] <- sample(
      observed_y2,
      sum(is.na(imputed$y2)),
      replace=TRUE
    )
    imputed$y2 <- ordered(imputed$y2, levels=levels(data$y2))

    observed_x2 <- imputed$x2[!is.na(imputed$x2)]
    n_missing_x2 <- sum(is.na(imputed$x2))
    imputed$x2[is.na(imputed$x2)] <- sample(
      observed_x2,
      n_missing_x2,
      replace=TRUE
    ) + stats::rnorm(n_missing_x2, sd=0.05)

    imputed
  })
}

make_mixed_indicator_group_mi_data <- function() {
  data <- make_mixed_indicator_group_data()
  n <- nrow(data)

  set.seed(2604)
  missing_y2 <- sample(seq_len(n), 75)
  missing_x2 <- sample(setdiff(seq_len(n), missing_y2), 70)
  data$y2[missing_y2] <- NA
  data$x2[missing_x2] <- NA

  lapply(seq_len(3), function(i) {
    imputed <- data
    set.seed(9400 + i)

    for(g in levels(data$group)) {
      missing_y2_g <- is.na(imputed$y2) & imputed$group == g
      observed_y2_g <- imputed$y2[!is.na(imputed$y2) & imputed$group == g]
      imputed$y2[missing_y2_g] <- sample(
        observed_y2_g,
        sum(missing_y2_g),
        replace=TRUE
      )

      missing_x2_g <- is.na(imputed$x2) & imputed$group == g
      observed_x2_g <- imputed$x2[!is.na(imputed$x2) & imputed$group == g]
      n_missing_x2_g <- sum(missing_x2_g)
      imputed$x2[missing_x2_g] <- sample(
        observed_x2_g,
        n_missing_x2_g,
        replace=TRUE
      ) + stats::rnorm(n_missing_x2_g, sd=0.05)
    }

    imputed$y2 <- ordered(imputed$y2, levels=levels(data$y2))
    imputed
  })
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
    estimator="WLSMV",
    point.wls="design",
    mi.pooling="sample.statistics"
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
    estimator="WLSMV",
    point.wls="design",
    mi.pooling="sample.statistics"
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

test_that("mixed MI survey models pool thresholds, means, and correlations", {
  skip_if_not_installed("mitools")

  imputed_data <- make_mixed_indicator_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=imputation_list,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=8
  )
  fit <- fit_mixed_indicator_model(imputed_data[[1]])

  fit_survey <- suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit,
    survey.design=rep_design,
    estimator="WLSMV",
    point.wls="design",
    mi.pooling="sample.statistics"
  ))
  survey_info <- attr(fit_survey, "lavaan.survey.info")

  per_imputation <- lapply(
    rep_design$designs,
    lavaan.survey:::get.ordinal.survey.stats,
    ov.names=c("y1", "y2", "x1", "x2"),
    ordered=c("y1", "y2"),
    ngroups=1
  )
  expected <- lavaan.survey:::pool.ordinal.mi.stats(
    per_imputation,
    ngroups=1
  )

  sampstat <- lavaan::lavInspect(fit_survey, "sampstat")
  Gamma <- lavaan::lavTech(fit_survey, "gamma")[[1]]

  expect_equal(survey_info$mode, "mixed ordinal/continuous")
  expect_equal(survey_info$mi.pooling, "sample.statistics")
  expect_equal(survey_info$point.wls, "design")
  expect_true(lavaan::lavInspect(fit_survey, "converged"))
  expect_equal(names(expected$point.stats$wls.obs),
               c("y1|t1", "y1|t2", "y2|t1", "y2|t2", "x1~1", "x2~1",
                 "x1~~x1", "x2~~x2", "y1~~y2", "y1~~x1", "y1~~x2",
                 "y2~~x1", "y2~~x2", "x1~~x2"))
  expect_equal(unclass(as.matrix(sampstat$cov)),
               unclass(as.matrix(expected$point.stats$sample.cov)),
               tolerance=1e-10, check.attributes=FALSE)
  expect_equal(as.numeric(sampstat$mean),
               as.numeric(expected$point.stats$sample.mean),
               tolerance=1e-10, check.attributes=FALSE)
  expect_equal(as.numeric(sampstat$th),
               as.numeric(expected$point.stats$sample.th),
               tolerance=1e-10, check.attributes=FALSE)
  expect_equal(Gamma, expected$Gamma,
               tolerance=1e-10, check.attributes=FALSE)
  expect_true(is.finite(lavaan::fitMeasures(fit_survey, "chisq.scaled")))
})

test_that("mixed multiple-group MI models preserve invariance constraints", {
  skip_if_not_installed("mitools")

  group_labels <- c("Girls", "Boys")
  imputed_data <- make_mixed_indicator_group_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=imputation_list,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=8
  )
  fit <- fit_mixed_indicator_model(
    imputed_data[[1]],
    group="group",
    group.equal=c("loadings", "thresholds", "intercepts")
  )

  fit_survey <- suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit,
    survey.design=rep_design,
    estimator="WLSMV",
    point.wls="design",
    mi.pooling="sample.statistics"
  ))

  per_imputation <- lapply(
    rep_design$designs,
    lavaan.survey:::get.ordinal.survey.stats,
    ov.names=c("y1", "y2", "x1", "x2"),
    ordered=c("y1", "y2"),
    ngroups=2,
    group.var="group",
    group.labels=group_labels
  )
  expected <- lavaan.survey:::pool.ordinal.mi.stats(
    per_imputation,
    ngroups=2,
    group.labels=group_labels
  )

  sampstat <- lavaan::lavInspect(fit_survey, "sampstat")
  Gamma <- lavaan::lavTech(fit_survey, "gamma")
  partable <- lavaan::parTable(fit_survey)
  constraints <- partable[partable$op == "==", ]
  loadings <- partable[partable$op == "=~" & partable$rhs != "y1", ]
  thresholds <- partable[partable$op == "|", ]
  intercepts <- partable[partable$op == "~1" & partable$lhs %in% c("x1", "x2"), ]

  expect_true(lavaan::lavInspect(fit_survey, "converged"))
  expect_equal(lavaan::lavInspect(fit_survey, "ngroups"), 2)
  expect_equal(lavaan::lavInspect(fit_survey, "nobs"),
               unname(expected$point.stats$nobs),
               check.attributes=FALSE)
  expect_equal(names(expected$point.stats$wls.obs), group_labels)
  expect_gt(nrow(constraints), 0)
  expect_true(any(loadings$label != ""))
  expect_true(any(thresholds$label != ""))
  expect_true(any(intercepts$label != ""))

  for(g in seq_along(group_labels)) {
    expect_equal(names(expected$point.stats$wls.obs[[g]]),
                 c("y1|t1", "y1|t2", "y2|t1", "y2|t2", "x1~1", "x2~1",
                   "x1~~x1", "x2~~x2", "y1~~y2", "y1~~x1", "y1~~x2",
                   "y2~~x1", "y2~~x2", "x1~~x2"))
    expect_equal(unclass(as.matrix(sampstat[[g]]$cov)),
                 unclass(as.matrix(expected$point.stats$sample.cov[[g]])),
                 tolerance=1e-10, check.attributes=FALSE)
    expect_equal(as.numeric(sampstat[[g]]$mean),
                 as.numeric(expected$point.stats$sample.mean[[g]]),
                 tolerance=1e-10, check.attributes=FALSE)
    expect_equal(as.numeric(sampstat[[g]]$th),
                 as.numeric(expected$point.stats$sample.th[[g]]),
                 tolerance=1e-10, check.attributes=FALSE)
    expect_equal(Gamma[[g]], expected$Gamma[[g]],
                 tolerance=1e-10, check.attributes=FALSE)
  }
  expect_true(is.finite(lavaan::fitMeasures(fit_survey, "chisq.scaled")))
})

test_that("mixed survey CFA can use lavaan WLS weights for point estimation", {
  data <- make_mixed_indicator_group_data()
  rep_design <- make_mixed_indicator_rep_design(data, replicates=8)

  stats_design <- lavaan.survey:::get.ordinal.survey.stats(
    rep.design=rep_design,
    ov.names=c("y1", "y2", "x1", "x2"),
    ordered=c("y1", "y2"),
    ngroups=2,
    group.var="group",
    group.labels=c("Girls", "Boys"),
    point.wls="design"
  )
  stats_lavaan <- lavaan.survey:::get.ordinal.survey.stats(
    rep.design=rep_design,
    ov.names=c("y1", "y2", "x1", "x2"),
    ordered=c("y1", "y2"),
    ngroups=2,
    group.var="group",
    group.labels=c("Girls", "Boys"),
    point.wls="lavaan"
  )

  for(g in c("Girls", "Boys")) {
    expect_equal(stats_lavaan$WLS.V[[g]],
                 stats_lavaan$point.stats$sample.wls.v[[g]],
                 tolerance=1e-10, check.attributes=FALSE)
    expect_false(isTRUE(all.equal(stats_design$WLS.V[[g]],
                                  stats_lavaan$WLS.V[[g]],
                                  tolerance=1e-7,
                                  check.attributes=FALSE)))
  }
})

test_that("mixed single-group MI models support Rubin parameter pooling", {
  skip_if_not_installed("mitools")

  imputed_data <- make_mixed_indicator_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=imputation_list,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=4
  )
  fit <- fit_mixed_indicator_model(imputed_data[[1]])

  fit_pooled <- suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit,
    survey.design=rep_design,
    estimator="WLSMV",
    point.wls="lavaan",
    mi.pooling="parameters",
    within.variance="naive"
  ))
  per_imputation <- lapply(rep_design$designs, function(design_i) {
    suppressWarnings(lavaan.survey:::fit.ordinal.weighted.lavaan(
      design=design_i,
      lavaan.fit=fit,
      ordered=c("y1", "y2"),
      estimator="WLSMV"
    ))
  })
  expected <- lavaan.survey:::pool.lavaan.mi.parameters(per_imputation)

  expect_s3_class(fit_pooled, "lavaan.survey.mi")
  expect_equal(fit_pooled$m, length(imputed_data))
  survey_info <- attr(fit_pooled, "lavaan.survey.info")
  expect_equal(survey_info$within.variance, "naive")
  expect_equal(coef(fit_pooled), coef(expected), tolerance=1e-10)
  expect_equal(vcov(fit_pooled), vcov(expected), tolerance=1e-10,
               check.attributes=FALSE)
  expect_true(all(is.finite(coef(fit_pooled))))
  expect_true(all(diag(vcov(fit_pooled)) > 0))
  expect_true(all(is.finite(fit_pooled$fit.measures)))
})

test_that("mixed MI parameter pooling reports the failing imputation", {
  skip_if_not_installed("mitools")

  imputed_data <- make_mixed_indicator_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=imputation_list,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=4
  )
  rep_design$designs[[2]]$variables$y1 <- NULL
  fit <- fit_mixed_indicator_model(imputed_data[[1]])

  expect_error(
    suppressWarnings(lavaan.survey:::pool.ordinal.mi.parameters(
      lavaan.fit=fit,
      rep.design=rep_design,
      ordered=c("y1", "y2"),
      estimator="WLSMV",
      point.wls="lavaan",
      rep.type="bootstrap",
      replicates=4,
      within.variance="naive"
    )),
    "Imputation 2 failed:"
  )
})

test_that("mixed MI parameter pooling can use replicate within variance", {
  skip_if_not_installed("mitools")

  imputed_data <- make_mixed_indicator_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=imputation_list,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=4
  )
  fit <- fit_mixed_indicator_model(imputed_data[[1]])

  progress_messages <- character()
  fit_pooled <- withCallingHandlers(
    suppressWarnings(lavaan.survey.ordinal(
      lavaan.fit=fit,
      survey.design=rep_design,
      estimator="WLSMV",
      point.wls="lavaan",
      mi.pooling="parameters",
      within.variance="replicate",
      verbose=TRUE
    )),
    message=function(m) {
      progress_messages <<- c(progress_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("Replicate fits for imputation 1/3: 1/4",
                        progress_messages,
                        fixed=TRUE)))
  per_imputation <- lapply(rep_design$designs, function(design_i) {
    suppressWarnings(lavaan.survey:::fit.ordinal.weighted.lavaan(
      design=design_i,
      lavaan.fit=fit,
      ordered=c("y1", "y2"),
      estimator="WLSMV"
    ))
  })
  replicate.vcov <- Map(
    lavaan.survey:::estimate.replicate.parameter.vcov,
    design=rep_design$designs,
    point.fit=per_imputation,
    MoreArgs=list(lavaan.fit=fit, ordered=c("y1", "y2"), estimator="WLSMV")
  )
  expected <- lavaan.survey:::pool.lavaan.mi.parameters(
    per_imputation,
    vcov.list=replicate.vcov
  )
  expect_error(
    lavaan.survey:::pool.lavaan.mi.parameters(
      per_imputation,
      vcov.list=replicate.vcov[-1]
    ),
    "vcov.list \\(.*\\) and fits \\(.*\\) must have the same length"
  )

  expect_s3_class(fit_pooled, "lavaan.survey.mi")
  survey_info <- attr(fit_pooled, "lavaan.survey.info")
  expect_equal(survey_info$within.variance, "replicate")
  expect_equal(fit_pooled$df.complete,
               min(vapply(rep_design$designs, survey::degf, numeric(1))))
  expect_equal(coef(fit_pooled), coef(expected), tolerance=1e-10)
  expect_equal(vcov(fit_pooled), vcov(expected), tolerance=1e-10,
               check.attributes=FALSE)
  expect_true(all(diag(fit_pooled$vcov.within) >= 0))

  summary_table <- capture.output(summary_out <- summary(fit_pooled))
  expect_true(length(summary_table) > 0)
  expect_true(any(grepl("lavaan.survey MI Summary", summary_table,
                        fixed=TRUE)))
  expect_true(any(grepl("Pooled Fit Measures", summary_table,
                        fixed=TRUE)))
  expect_true(any(grepl("Latent Variables", summary_table,
                        fixed=TRUE)))
  expect_true(any(grepl("P(>|t|)", summary_table,
                        fixed=TRUE)))
  expect_true(all(c("Estimate", "Std.Err", "t.value", "df", "P.value") %in%
                    names(summary_out)))
  expect_true(all(c("lhs", "op", "rhs") %in% names(summary_out)))
  expect_false("z.value" %in% names(summary_out))
  expected_df <- lavaan.survey:::barnard.rubin.df(fit_pooled)
  expected_p <- 2 * stats::pt(abs(summary_out$t.value),
                              df=expected_df,
                              lower.tail=FALSE)
  normal_ref <- is.infinite(expected_df)
  expected_p[normal_ref] <- 2 * stats::pnorm(abs(summary_out$t.value[normal_ref]),
                                             lower.tail=FALSE)
  expect_equal(summary_out$df, expected_df, tolerance=1e-10,
               check.attributes=FALSE)
  expect_equal(summary_out$P.value, expected_p, tolerance=1e-10)

  pe <- parameterEstimates(fit_pooled, remove.nonfree=TRUE)
  expect_true(all(c("lhs", "op", "rhs", "est", "se", "t", "df",
                    "pvalue") %in% names(pe)))
  expect_equal(pe$est, as.numeric(coef(fit_pooled)), tolerance=1e-10,
               check.attributes=FALSE)
  expect_equal(pe$df, expected_df, tolerance=1e-10,
               check.attributes=FALSE)
  expect_equal(pe$pvalue, expected_p, tolerance=1e-10)

  fm <- fitMeasures(fit_pooled, c("cfi.scaled", "rmsea.scaled"))
  expect_equal(fm,
               fit_pooled$fit.measures[c("cfi.scaled", "rmsea.scaled")],
               check.attributes=FALSE)
  fm.df <- fitMeasures(fit_pooled, c("cfi.scaled", "rmsea.scaled"),
                       output="data.frame")
  expect_equal(fm.df$measure, c("cfi.scaled", "rmsea.scaled"))
  expect_equal(fm.df$value, as.numeric(fm), check.attributes=FALSE)
})

test_that("mixed MI defaults to the Mplus-nearer parameter-pooling path", {
  skip_if_not_installed("mitools")

  imputed_data <- make_mixed_indicator_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=imputation_list,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=4
  )
  fit <- fit_mixed_indicator_model(imputed_data[[1]])

  fit_default <- expect_message(suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit,
    survey.design=rep_design,
    estimator="WLSMV"
  )), "lavaan\\.survey\\.ordinal mode: mixed ordinal/continuous")

  expect_s3_class(fit_default, "lavaan.survey.mi")
  survey_info <- attr(fit_default, "lavaan.survey.info")
  expect_null(fit_default$point.wls)
  expect_null(fit_default$mi.pooling)
  expect_null(fit_default$within.variance)
  expect_null(fit_default$survey.info)
  expect_equal(survey_info$mode, "mixed ordinal/continuous")
  expect_equal(survey_info$mi.pooling, "parameters")
  expect_equal(survey_info$point.wls, "lavaan")
  expect_equal(survey_info$within.variance, "replicate")
  printed <- utils::capture.output(print(fit_default))
  expect_true(any(grepl("lavaan.survey.ordinal mode: mixed ordinal/continuous",
                        printed, fixed=TRUE)))
  expect_true(any(grepl("MI pooling: parameters", printed, fixed=TRUE)))
  expect_true(any(grepl("Point WLS: lavaan", printed, fixed=TRUE)))
  expect_true(any(grepl("Within variance: replicate", printed, fixed=TRUE)))
  legacy <- fit_default
  attr(legacy, "lavaan.survey.info") <- NULL
  legacy$survey.info <- survey_info
  legacy_printed <- utils::capture.output(print(legacy))
  expect_true(any(grepl("MI pooling: parameters", legacy_printed, fixed=TRUE)))
  expect_true(all(is.finite(coef(fit_default))))
})

test_that("lavaan.survey dispatches mixed models to the ordinal wrapper", {
  data <- make_mixed_indicator_data()
  rep_design <- make_mixed_indicator_rep_design(data, replicates=6)
  fit <- fit_mixed_indicator_model(data, meanstructure=FALSE)

  fit_wrapper <- expect_message(suppressWarnings(lavaan.survey(
    lavaan.fit=fit,
    survey.design=rep_design,
    point.wls="design",
    mi.pooling="sample.statistics"
  )), "lavaan\\.survey\\.ordinal mode: mixed ordinal/continuous")
  survey_info <- attr(fit_wrapper, "lavaan.survey.info")

  expect_s4_class(fit_wrapper, "lavaan")
  expect_equal(survey_info$mode, "mixed ordinal/continuous")
  expect_equal(survey_info$mi.pooling, "none")
  expect_equal(survey_info$point.wls, "design")
  expect_equal(lavaan::lavInspect(fit_wrapper, "options")$estimator, "DWLS")
})

test_that("lavaan.survey dispatches mixed MI models to parameter pooling by default", {
  skip_if_not_installed("mitools")

  imputed_data <- make_mixed_indicator_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=imputation_list,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=4
  )
  fit <- fit_mixed_indicator_model(imputed_data[[1]])

  fit_wrapper <- expect_message(suppressWarnings(lavaan.survey(
    lavaan.fit=fit,
    survey.design=rep_design
  )), "MI pooling: parameters")

  expect_s3_class(fit_wrapper, "lavaan.survey.mi")
  survey_info <- attr(fit_wrapper, "lavaan.survey.info")
  expect_equal(survey_info$mode, "mixed ordinal/continuous")
  expect_equal(survey_info$mi.pooling, "parameters")
  expect_equal(survey_info$point.wls, "lavaan")
  expect_equal(survey_info$within.variance, "replicate")
})

test_that("mixed multiple-group MI models support Rubin parameter pooling", {
  skip_if_not_installed("mitools")

  imputed_data <- make_mixed_indicator_group_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=imputation_list,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=4
  )
  fit <- fit_mixed_indicator_model(
    imputed_data[[1]],
    group="group",
    group.equal=c("loadings", "thresholds", "intercepts")
  )

  fit_pooled <- suppressWarnings(lavaan.survey.ordinal(
    lavaan.fit=fit,
    survey.design=rep_design,
    estimator="WLSMV",
    point.wls="lavaan",
    mi.pooling="parameters",
    within.variance="naive"
  ))
  per_imputation <- lapply(rep_design$designs, function(design_i) {
    suppressWarnings(lavaan.survey:::fit.ordinal.weighted.lavaan(
      design=design_i,
      lavaan.fit=fit,
      ordered=c("y1", "y2"),
      estimator="WLSMV"
    ))
  })
  expected <- lavaan.survey:::pool.lavaan.mi.parameters(per_imputation)

  expect_s3_class(fit_pooled, "lavaan.survey.mi")
  expect_equal(fit_pooled$m, length(imputed_data))
  survey_info <- attr(fit_pooled, "lavaan.survey.info")
  expect_equal(survey_info$within.variance, "naive")
  expect_equal(coef(fit_pooled), coef(expected), tolerance=1e-10)
  expect_equal(vcov(fit_pooled), vcov(expected), tolerance=1e-10,
               check.attributes=FALSE)
  expect_true(any(grepl("\\.p", names(coef(fit_pooled)))))
  expect_true(all(diag(vcov(fit_pooled)) > 0))
})

test_that("lavaan robust within variance is guarded for ordered models", {
  skip_if_not_installed("mitools")

  imputed_data <- make_mixed_indicator_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    strata=~stratum,
    weights=~weight,
    data=imputation_list,
    nest=TRUE
  )
  rep_design <- survey::as.svrepdesign(
    design,
    type="bootstrap",
    replicates=4
  )
  fit <- fit_mixed_indicator_model(imputed_data[[1]])

  expect_error(
    lavaan.survey.ordinal(
      lavaan.fit=fit,
      survey.design=rep_design,
      estimator="WLSMV",
      point.wls="lavaan",
      mi.pooling="parameters",
      within.variance="lavaan.robust"
    ),
    "categorical \\+ clustered estimation is not supported"
  )
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
