validation_dir <- file.path("validation", "mplus-demo")
if (!dir.exists(validation_dir)) {
  stop("Run this script from the package root.")
}

required <- c("MASS", "lavaan", "survey", "mitools", "mice")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Install required packages first: ", paste(missing, collapse = ", "))
}

load_lavaan_survey <- function() {
  if (file.exists(file.path("R", "lavaan.survey.R"))) {
    library(lavaan)
    library(survey)
    source(file.path("R", "lavaan.survey.R"))
  } else if (requireNamespace("lavaan.survey", quietly = TRUE)) {
    library(lavaan.survey)
  } else {
    stop("Install lavaan.survey or run this script from the package root.")
  }
}

to_ordered4 <- function(x) {
  ordered(
    cut(x, breaks = c(-Inf, -0.8, 0, 0.8, Inf), labels = 1:4),
    levels = 1:4
  )
}

to_mplus_integer_data <- function(data, items) {
  out <- data
  for (item in items) {
    out[[item]] <- as.integer(as.character(out[[item]]))
  }
  out
}

load_lavaan_survey()

set.seed(20260506)

n_clusters <- 60L
cluster_size <- 10L
n <- n_clusters * cluster_size
items <- paste0("y", 1:6)
m <- 10L

latent_cor <- matrix(c(1.0, 0.35, 0.35, 1.0), 2, 2)
eta <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = latent_cor)
lambda <- c(0.85, 0.78, 0.72, 0.82, 0.76, 0.70)

y_star <- data.frame(
  y1 = lambda[1] * eta[, 1] + sqrt(1 - lambda[1]^2) * rnorm(n),
  y2 = lambda[2] * eta[, 1] + sqrt(1 - lambda[2]^2) * rnorm(n),
  y3 = lambda[3] * eta[, 1] + sqrt(1 - lambda[3]^2) * rnorm(n),
  y4 = lambda[4] * eta[, 2] + sqrt(1 - lambda[4]^2) * rnorm(n),
  y5 = lambda[5] * eta[, 2] + sqrt(1 - lambda[5]^2) * rnorm(n),
  y6 = lambda[6] * eta[, 2] + sqrt(1 - lambda[6]^2) * rnorm(n)
)

dat <- as.data.frame(lapply(y_star, to_ordered4))
names(dat) <- items

cluster <- rep(seq_len(n_clusters), each = cluster_size)
stratum_by_cluster <- rep(seq_len(6), each = n_clusters / 6)
stratum <- stratum_by_cluster[cluster]
base_weight <- c(0.75, 0.90, 1.05, 1.20, 1.35, 1.50)[stratum]
weight <- round(base_weight * runif(n, min = 0.85, max = 1.15), 6)

dat$wgt <- weight
dat$clu <- cluster
dat$str <- stratum

# Create moderate missingness in three ordered indicators. The imputation model
# includes survey design information so that the analysis comparison focuses on
# lavaan.survey versus Mplus rather than on different imputers.
dat_missing <- dat
p_y2 <- plogis(-1.45 + 0.45 * scale(y_star$y1)[, 1] + 0.20 * (dat$str >= 4))
p_y5 <- plogis(-1.35 + 0.45 * scale(y_star$y4)[, 1] - 0.15 * scale(dat$wgt)[, 1])
p_y6 <- rep(0.15, n)
dat_missing$y2[runif(n) < p_y2] <- NA
dat_missing$y5[runif(n) < p_y5] <- NA
dat_missing$y6[runif(n) < p_y6] <- NA

ini <- mice::mice(dat_missing, maxit = 0, printFlag = FALSE)
method <- ini$method
method[] <- ""
method[c("y2", "y5", "y6")] <- "polr"

predictor_matrix <- ini$predictorMatrix
predictor_matrix[,] <- 0
predictor_matrix[c("y2", "y5", "y6"), c(items, "wgt", "str")] <- 1
diag(predictor_matrix) <- 0

imp <- mice::mice(
  dat_missing,
  m = m,
  maxit = 10,
  method = method,
  predictorMatrix = predictor_matrix,
  seed = 20260507,
  printFlag = FALSE
)

imputed_data <- lapply(seq_len(m), function(i) {
  out <- mice::complete(imp, action = i)
  for (item in items) {
    out[[item]] <- ordered(out[[item]], levels = 1:4)
  }
  out
})

imputation_files <- sprintf("ordinal_mi_imp%02d.dat", seq_len(m))
for (i in seq_along(imputed_data)) {
  mplus_data <- to_mplus_integer_data(imputed_data[[i]], items)
  utils::write.table(
    mplus_data[, c(items, "wgt", "clu", "str")],
    file = file.path(validation_dir, imputation_files[[i]]),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

writeLines(
  imputation_files,
  con = file.path(validation_dir, "ordinal_mi_implist.dat")
)

writeLines(c(
  "TITLE: Ordinal MI complex CFA for lavaan.survey validation;",
  "DATA:",
  "  FILE IS ordinal_mi_implist.dat;",
  "  TYPE = IMPUTATION;",
  "VARIABLE:",
  "  NAMES ARE y1 y2 y3 y4 y5 y6 wgt clu str;",
  "  USEVARIABLES ARE y1 y2 y3 y4 y5 y6;",
  "  CATEGORICAL ARE y1 y2 y3 y4 y5 y6;",
  "  WEIGHT IS wgt;",
  "  CLUSTER IS clu;",
  "  STRATIFICATION IS str;",
  "ANALYSIS:",
  "  TYPE = COMPLEX;",
  "  ESTIMATOR = WLSMV;",
  "  PARAMETERIZATION = DELTA;",
  "MODEL:",
  "  f1 BY y1 y2 y3;",
  "  f2 BY y4 y5 y6;",
  "  f1 WITH f2;",
  "OUTPUT:",
  "  STDYX TECH1;"
), con = file.path(validation_dir, "ordinal_mi_complex.inp"))

cfa_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
"

design <- survey::svydesign(
  ids = ~clu,
  strata = ~str,
  weights = ~wgt,
  data = mitools::imputationList(imputed_data),
  nest = TRUE
)

fit_naive <- lavaan::cfa(
  model = cfa_model,
  data = imputed_data[[1]],
  ordered = items,
  estimator = "WLSMV"
)

fit_survey <- lavaan.survey.ordinal(
  lavaan.fit = fit_naive,
  survey.design = design,
  estimator = "WLSMV"
)

pe <- lavaan::parameterEstimates(fit_survey, standardized = TRUE)
pe <- pe[pe$op %in% c("=~", "~~", "|"), ]

utils::write.csv(
  pe[, c("lhs", "op", "rhs", "est", "se", "z", "pvalue", "std.all")],
  file = file.path(validation_dir, "ordinal_mi_lavaan_survey_parameters.csv"),
  row.names = FALSE
)

fm <- lavaan::fitMeasures(
  fit_survey,
  c("chisq", "df", "pvalue",
    "chisq.scaled", "df.scaled", "pvalue.scaled",
    "cfi.scaled", "tli.scaled", "rmsea.scaled",
    "cfi.robust", "tli.robust", "rmsea.robust", "srmr")
)

utils::write.csv(
  data.frame(measure = names(fm), value = unname(fm)),
  file = file.path(validation_dir, "ordinal_mi_lavaan_survey_fit.csv"),
  row.names = FALSE
)

missing_summary <- data.frame(
  variable = items,
  missing_proportion = vapply(
    dat_missing[items],
    function(x) mean(is.na(x)),
    numeric(1)
  )
)
utils::write.csv(
  missing_summary,
  file = file.path(validation_dir, "ordinal_mi_missing_summary.csv"),
  row.names = FALSE
)

message("Ordinal MI validation files written to: ", normalizePath(validation_dir))
message(
  "Next step: run `/Applications/MplusDemo/mpdemo ",
  "ordinal_mi_complex.inp` inside that folder."
)
