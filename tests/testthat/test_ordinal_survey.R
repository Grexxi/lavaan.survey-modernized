suppressPackageStartupMessages(suppressWarnings(library(lavaan.survey)))

simulate_ordinal_survey_data <- function(n=600) {
  Sigma <- matrix(0, nrow=10, ncol=10)
  Sigma[1:5, 1:5] <- 0.6
  Sigma[6:10, 6:10] <- 0.6
  diag(Sigma) <- 1

  raw <- MASS::mvrnorm(n=n, mu=rep(0, 10), Sigma=Sigma)
  dat <- as.data.frame(raw)
  names(dat) <- paste0("v", 1:10)

  for(i in 1:5) {
    dat[[i]] <- cut(dat[[i]],
                    breaks=c(-Inf, -0.5, 0.5, Inf),
                    labels=c(1, 2, 3),
                    ordered_result=TRUE)
  }
  for(i in 6:10) {
    dat[[i]] <- cut(dat[[i]],
                    breaks=c(-Inf, -1, -0.5, 0.5, 1, Inf),
                    labels=c(1, 2, 3, 4, 5),
                    ordered_result=TRUE)
  }

  dat$cluster <- factor(rep(seq_len(30), each=n/30))
  dat$stratum <- factor(ifelse(as.numeric(dat$cluster) <= 15, 1, 2))
  dat$weight <- round(stats::runif(n, min=0.5, max=2.0), 2)
  dat$gender <- factor(sample(c("Girls", "Boys"), n, replace=TRUE))

  dat
}

.ordinal_survey_fits <- NULL

fit_ordinal_survey_models <- function() {
  if(!is.null(.ordinal_survey_fits)) return(.ordinal_survey_fits)

  set.seed(42)
  dat <- simulate_ordinal_survey_data()
  items <- paste0("v", 1:10)
  model <- "
    f1 =~ v1 + v2 + v3 + v4 + v5
    f2 =~ v6 + v7 + v8 + v9 + v10
  "
  design <- survey::svydesign(id=~cluster,
                              strata=~stratum,
                              weights=~weight,
                              data=dat,
                              nest=TRUE)

  fit_basis <- lavaan::cfa(model=model, data=dat, ordered=items)
  fit_group <- lavaan::cfa(model=model, data=dat, ordered=items,
                           group="gender")
  fit_scalar <- lavaan::cfa(model=model, data=dat, ordered=items,
                            group="gender",
                            group.equal=c("loadings", "thresholds"))

  .ordinal_survey_fits <<- suppressWarnings(list(
    base=lavaan.survey.ordinal(fit_basis, design, estimator="WLSMV"),
    group=lavaan.survey.ordinal(fit_group, design, estimator="WLSMV"),
    scalar=lavaan.survey.ordinal(fit_scalar, design, estimator="WLSMV")
  ))
  .ordinal_survey_fits
}

test_that("ordinal survey CFA works for single-group models", {
  fits <- fit_ordinal_survey_models()
  fit <- fits$base

  expect_true(lavaan::lavInspect(fit, "converged"))
  expect_equal(lavaan::lavInspect(fit, "ngroups"), 1)
  expect_true(is.finite(lavaan::fitMeasures(fit, "chisq.scaled")))
  expect_true(is.finite(lavaan::fitMeasures(fit, "cfi.robust")))
  expect_equal(lavaan::lavInspect(fit, "options")$se, "robust.sem")
  expect_true("scaled.shifted" %in% lavaan::lavInspect(fit, "options")$test)

  pe <- lavaan::parameterEstimates(fit)
  estimable <- !is.na(pe$se) & pe$se > 0 & pe$op != "|"
  expect_gt(sum(estimable), 0)
  expect_true(all(is.finite(pe$se[estimable])))
})

test_that("ordinal survey CFA works for multiple-group models", {
  fits <- fit_ordinal_survey_models()
  fit <- fits$group

  expect_true(lavaan::lavInspect(fit, "converged"))
  expect_equal(lavaan::lavInspect(fit, "ngroups"), 2)
  expect_equal(length(lavaan::lavTech(fit, "gamma")), 2)
  expect_true(is.finite(lavaan::fitMeasures(fit, "chisq.scaled")))
  expect_true("scaled.shifted" %in% lavaan::lavInspect(fit, "options")$test)

  pe <- lavaan::parameterEstimates(fit)
  loadings <- pe[pe$op == "=~", ]
  expect_equal(length(unique(loadings$group)), 2)
  free.loadings <- loadings[!is.na(loadings$se) & loadings$se > 0, ]
  expect_gt(nrow(free.loadings), 0)
  expect_true(all(is.finite(free.loadings$se)))
})

test_that("ordinal survey CFA preserves loading and threshold invariance constraints", {
  fits <- fit_ordinal_survey_models()
  fit <- fits$scalar

  expect_true(lavaan::lavInspect(fit, "converged"))
  expect_equal(lavaan::lavInspect(fit, "ngroups"), 2)
  expect_true(is.finite(lavaan::fitMeasures(fit, "chisq.scaled")))

  partable <- lavaan::parTable(fit)
  constraints <- partable[partable$op == "==", ]
  expect_gt(nrow(constraints), 0)

  loadings <- partable[partable$op == "=~" & partable$rhs != "v1" & partable$rhs != "v6", ]
  thresholds <- partable[partable$op == "|", ]
  expect_true(any(loadings$label != ""))
  expect_true(any(thresholds$label != ""))

  pe <- lavaan::parameterEstimates(fit)
  constrained.pe <- !is.na(pe$se) & pe$se > 0 & pe$op %in% c("=~", "|")
  expect_gt(sum(constrained.pe), 0)
  expect_true(all(is.finite(pe$se[constrained.pe])))
})
