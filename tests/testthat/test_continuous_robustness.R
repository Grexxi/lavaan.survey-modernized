suppressPackageStartupMessages(suppressWarnings(library(lavaan.survey)))

context("continuous survey robustness")

make_hs_design <- function() {
  data("HolzingerSwineford1939", package="lavaan")
  dat <- lavaan::HolzingerSwineford1939
  set.seed(2501)
  dat$cluster <- seq_len(nrow(dat))
  dat$weight <- stats::runif(nrow(dat), 0.8, 1.2)

  list(
    data=dat,
    design=survey::svydesign(ids=~cluster, weights=~weight, data=dat)
  )
}

test_that("Yuan-Bentler gamma works without a meanstructure", {
  hs <- make_hs_design()
  local_model <- "visual =~ x1 + x2 + x3"
  fit <- lavaan::cfa(
    model=local_model,
    data=hs$data,
    estimator="MLM",
    meanstructure=FALSE
  )

  fit_svy <- lavaan.survey(
    lavaan.fit=fit,
    survey.design=hs$design,
    estimator="MLM",
    estimator.gamma="Yuan-Bentler"
  )

  Gamma <- lavaan::lavTech(fit_svy, "gamma")
  if(is.list(Gamma)) Gamma <- Gamma[[1]]

  expect_true(lavaan::lavInspect(fit_svy, "converged"))
  expect_equal(dim(Gamma), c(6L, 6L))
})

test_that("group labels containing quotes subset safely", {
  hs <- make_hs_design()
  hs$data$grp_quote <- factor(ifelse(seq_len(nrow(hs$data)) <= 150,
                                     "won't", "ok"))
  hs$design <- survey::svydesign(ids=~cluster, weights=~weight, data=hs$data)
  local_model <- "visual =~ x1 + x2 + x3"

  fit <- lavaan::cfa(
    model=local_model,
    data=hs$data,
    group="grp_quote",
    estimator="MLM",
    meanstructure=TRUE
  )
  fit_svy <- lavaan.survey(fit, hs$design, estimator="MLM")

  expect_true(lavaan::lavInspect(fit_svy, "converged"))
  expect_equal(lavaan::lavInspect(fit_svy, "nobs"),
               lavaan::lavInspect(fit, "nobs"))
})

test_that("DWLS weights ignore nonpositive Gamma diagonal entries", {
  Gamma <- diag(c(4, -1, 0, .Machine$double.eps / 2, 2))
  W <- lavaan.survey:::get.dwls.weight(Gamma)

  expect_equal(diag(W), c(0.25, 0, 0, 0, 0.5))
})

test_that("continuous MI multiple-group models keep group sample sizes", {
  skip_if_not_installed("mitools")

  hs <- make_hs_design()
  imputed_data <- lapply(seq_len(2), function(i) hs$data)
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    weights=~weight,
    data=imputation_list
  )
  local_model <- "visual =~ x1 + x2 + x3"

  fit <- lavaan::cfa(
    model=local_model,
    data=hs$data,
    group="school",
    estimator="MLM",
    meanstructure=TRUE
  )
  fit_svy <- lavaan.survey(fit, design, estimator="MLM")

  expect_true(lavaan::lavInspect(fit_svy, "converged"))
  expect_equal(lavaan::lavInspect(fit_svy, "nobs"),
               lavaan::lavInspect(fit, "nobs"))
})

test_that("continuous MI models work without a meanstructure", {
  skip_if_not_installed("mitools")

  hs <- make_hs_design()
  imputed_data <- lapply(seq_len(2), function(i) hs$data)
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    weights=~weight,
    data=imputation_list
  )
  local_model <- "visual =~ x1 + x2 + x3"

  fit <- lavaan::cfa(
    model=local_model,
    data=hs$data,
    estimator="MLM",
    meanstructure=FALSE
  )
  fit_svy <- lavaan.survey(fit, design, estimator="MLM")

  Gamma <- lavaan::lavTech(fit_svy, "gamma")
  if(is.list(Gamma)) Gamma <- Gamma[[1]]

  expect_true(lavaan::lavInspect(fit_svy, "converged"))
  expect_equal(dim(Gamma), c(6L, 6L))
})

test_that("pval.pFsum accepts current robust lavaan test names", {
  hs <- make_hs_design()
  local_model <- "visual =~ x1 + x2 + x3"
  fit <- lavaan::cfa(
    model=local_model,
    data=hs$data,
    estimator="MLMV"
  )
  fit_svy <- lavaan.survey(fit, hs$design, estimator="MLMV")

  p_value <- pval.pFsum(fit_svy, hs$design)

  expect_true("scaled.shifted" %in%
                lavaan::lavInspect(fit_svy, "options")$test)
  expect_true(is.numeric(p_value))
  expect_equal(length(p_value), 1L)
  expect_true(is.na(p_value) || (p_value >= 0 && p_value <= 1))
})

test_that("pval.pFsum works for constrained multiple-group models", {
  hs <- make_hs_design()
  local_model <- "
    visual =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
  "
  fit <- lavaan::cfa(
    model=local_model,
    data=hs$data,
    group="school",
    group.equal=c("loadings", "intercepts"),
    estimator="MLMV",
    meanstructure=TRUE
  )
  fit_svy <- lavaan.survey(fit, hs$design, estimator="MLMV")

  UGamma <- lavaan.survey:::get.ugamma.matrix(fit_svy)
  p_value <- pval.pFsum(fit_svy, hs$design)

  expect_true(lavaan::lavInspect(fit_svy, "converged"))
  expect_equal(lavaan::lavInspect(fit_svy, "ngroups"), 2L)
  expect_equal(nrow(UGamma), ncol(UGamma))
  expect_true(is.numeric(p_value))
  expect_equal(length(p_value), 1L)
  expect_true(is.na(p_value) || (p_value >= 0 && p_value <= 1))
})

test_that("U.Gamma lists are combined block-diagonally", {
  UGamma <- list(diag(2), diag(c(3, 4)))
  UGamma.block <- lavaan.survey:::as.ugamma.matrix(UGamma)

  expect_equal(dim(UGamma.block), c(4L, 4L))
  expect_equal(diag(UGamma.block), c(1, 1, 3, 4))
  expect_equal(UGamma.block[1:2, 3:4], matrix(0, 2, 2))
  expect_equal(UGamma.block[3:4, 1:2], matrix(0, 2, 2))
})

test_that("U.Gamma extraction rejects malformed inputs", {
  expect_error(
    lavaan.survey:::as.ugamma.matrix(list()),
    "square matrix"
  )
  expect_error(
    lavaan.survey:::as.ugamma.matrix(matrix(1, 2, 3)),
    "square matrix"
  )
  expect_error(
    lavaan.survey:::as.ugamma.matrix(list(diag(2), matrix(1, 2, 3))),
    "square matrix"
  )
})
