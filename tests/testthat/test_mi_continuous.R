suppressPackageStartupMessages(suppressWarnings(library(lavaan.survey)))

context("continuous multiple-imputation survey models")

make_continuous_mi_data <- function() {
  set.seed(2401)
  n <- 240
  cluster <- rep(seq_len(40), each=6)
  weight <- round(stats::runif(n, 0.7, 1.5), 4)
  factor_score <- stats::rnorm(n)

  dat <- data.frame(
    y1=0.85 * factor_score + stats::rnorm(n, sd=0.55),
    y2=0.75 * factor_score + stats::rnorm(n, sd=0.65),
    y3=0.80 * factor_score + stats::rnorm(n, sd=0.60),
    cluster=cluster,
    weight=weight
  )

  missing_y2 <- sample(seq_len(n), 48)
  dat$y2[missing_y2] <- NA

  lapply(seq_len(5), function(i) {
    imp <- dat
    set.seed(9000 + i)
    imp$y2[is.na(imp$y2)] <- mean(imp$y2, na.rm=TRUE) +
      stats::rnorm(sum(is.na(imp$y2)), sd=0.8) +
      (i - 3) * 0.35
    imp
  })
}

pool_continuous_mi_sample_stats <- function(designs, ov.names, sample.nobs) {
  ov.formula <- stats::as.formula(paste("~", paste(ov.names, collapse="+")))
  Dplus <- lavaan::lav_matrix_duplication_ginv(length(ov.names))

  per_imputation <- lapply(designs, function(design) {
    sample.cov <- as.matrix(survey::svyvar(ov.formula, design=design, na.rm=TRUE))
    Gamma.cov <- attr(sample.cov, "var")
    Gamma.cov <- Dplus %*% Gamma.cov %*% t(Dplus)
    attr(sample.cov, "var") <- NULL

    sample.mean <- survey::svymean(ov.formula, design=design, na.rm=TRUE)
    Gamma.mean <- attr(sample.mean, "var")
    sample.mean <- structure(as.vector(sample.mean), names=names(sample.mean))

    Gamma <- lavaan::lav_matrix_bdiag(Gamma.mean, Gamma.cov) * sample.nobs

    list(sample.cov=sample.cov, sample.mean=sample.mean, Gamma=Gamma)
  })

  m <- length(per_imputation)
  sample.cov.list <- lapply(per_imputation, `[[`, "sample.cov")
  sample.mean.list <- lapply(per_imputation, `[[`, "sample.mean")

  sample.cov <- Reduce(`+`, sample.cov.list) / m
  sample.mean <- Reduce(`+`, sample.mean.list) / m

  mean.df <- do.call(rbind, sample.mean.list)
  cov.df <- do.call(rbind, lapply(sample.cov.list, lavaan::lav_matrix_vech))

  Gamma.within <- Reduce(`+`, lapply(per_imputation, `[[`, "Gamma")) / m
  Gamma.between <- stats::cov(cbind(mean.df, cov.df))
  dimnames(Gamma.between) <- dimnames(Gamma.within)
  Gamma <- Gamma.within + sample.nobs * ((m + 1) / m) * Gamma.between

  list(sample.cov=sample.cov, sample.mean=sample.mean, Gamma=Gamma)
}

fit_continuous_mi_lavaan_model <- function(data) {
  local_model <- "f =~ y1 + y2 + y3"
  lavaan::cfa(
    model=local_model,
    data=data,
    estimator="MLM",
    meanstructure=TRUE
  )
}

test_that("continuous MI survey models pool sample statistics with Rubin scaling", {
  skip_if_not_installed("mitools")

  imputed_data <- make_continuous_mi_data()
  imputation_list <- mitools::imputationList(imputed_data)
  design <- survey::svydesign(
    ids=~cluster,
    weights=~weight,
    data=imputation_list
  )

  fit_naive <- fit_continuous_mi_lavaan_model(imputed_data[[1]])

  fit_mi <- lavaan.survey(
    lavaan.fit=fit_naive,
    survey.design=design,
    estimator="MLM"
  )

  expected <- pool_continuous_mi_sample_stats(
    designs=design$designs,
    ov.names=c("y1", "y2", "y3"),
    sample.nobs=240
  )

  sampstat <- lavaan::lavInspect(fit_mi, "sampstat")
  Gamma <- lavaan::lavTech(fit_mi, "gamma")[[1]]
  expected.cov <- expected$sample.cov
  if(isTRUE(lavaan::lavInspect(fit_mi, "options")$sample.cov.rescale)) {
    expected.cov <- expected.cov * (240 - 1) / 240
  }

  expect_true(lavaan::lavInspect(fit_mi, "converged"))
  expect_equal(as.numeric(sampstat$mean), as.numeric(expected$sample.mean),
               tolerance=1e-10, check.attributes=FALSE)
  expect_equal(unclass(as.matrix(sampstat$cov)),
               unclass(as.matrix(expected.cov)),
               tolerance=1e-10, check.attributes=FALSE)
  expect_equal(Gamma, expected$Gamma,
               tolerance=1e-10, check.attributes=FALSE)
})
