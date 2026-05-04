validation_dir <- file.path("validation", "mplus-demo")
if (!dir.exists(validation_dir)) {
  stop("Run this script from the package root.")
}

required <- c("MASS", "lavaan", "survey")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Install required packages first: ", paste(missing, collapse = ", "))
}

load_lavaan_survey <- function() {
  if (requireNamespace("lavaan.survey", quietly = TRUE) &&
      exists("lavaan.survey.ordinal", envir = asNamespace("lavaan.survey"))) {
    library(lavaan.survey)
  } else {
    source(file.path("R", "lavaan.survey.R"))
  }
}

download_file <- function(url, dest) {
  ok <- tryCatch({
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    status <- system2("curl", c("-L", "-k", "-s", "-o", dest, url))
    ok <- identical(status, 0L)
  }

  if (!ok || !file.exists(dest)) {
    stop("Could not download ", url)
  }
}

load_lavaan_survey()

# -------------------------------------------------------------------------
# 1. Official Mplus example 5.10 sanity check
# -------------------------------------------------------------------------

ex510_url <- "https://www.statmodel.com/usersguide/chap5/ex5.10.dat"
ex510_dat <- file.path(validation_dir, "ex5.10.dat")
download_file(ex510_url, ex510_dat)

writeLines(c(
  "TITLE: Official Mplus User Guide example 5.10;",
  "DATA: FILE IS ex5.10.dat;",
  "VARIABLE:",
  "  NAMES ARE u1a u1b u1c u2a u2b u2c;",
  "  CATEGORICAL ARE u1a u1b u1c u2a u2b u2c;",
  "MODEL:",
  "  f1 BY u1a u1b@1 u1c@1;",
  "  f2 BY u2a u2b@1 u2c@1;",
  "  [u1a$1 u1b$1 u1c$1] (1);",
  "  [u2a$1 u2b$1 u2c$1] (2);",
  "OUTPUT: STDYX TECH1;"
), con = file.path(validation_dir, "mplus_ex510.inp"))

# -------------------------------------------------------------------------
# 2. Simulated ordinal survey dataset for lavaan.survey vs Mplus
# -------------------------------------------------------------------------

set.seed(20260504)

n_clusters <- 60L
cluster_size <- 10L
n <- n_clusters * cluster_size
items <- paste0("y", 1:6)

latent_cor <- matrix(c(1.0, 0.35, 0.35, 1.0), 2, 2)
eta <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = latent_cor)
lambda <- c(0.85, 0.78, 0.72, 0.82, 0.76, 0.70)

y_star <- matrix(NA_real_, nrow = n, ncol = 6)
y_star[, 1] <- lambda[1] * eta[, 1] + sqrt(1 - lambda[1]^2) * rnorm(n)
y_star[, 2] <- lambda[2] * eta[, 1] + sqrt(1 - lambda[2]^2) * rnorm(n)
y_star[, 3] <- lambda[3] * eta[, 1] + sqrt(1 - lambda[3]^2) * rnorm(n)
y_star[, 4] <- lambda[4] * eta[, 2] + sqrt(1 - lambda[4]^2) * rnorm(n)
y_star[, 5] <- lambda[5] * eta[, 2] + sqrt(1 - lambda[5]^2) * rnorm(n)
y_star[, 6] <- lambda[6] * eta[, 2] + sqrt(1 - lambda[6]^2) * rnorm(n)

dat <- as.data.frame(y_star)
names(dat) <- items

for (item in items) {
  dat[[item]] <- cut(
    dat[[item]],
    breaks = c(-Inf, -0.8, 0, 0.8, Inf),
    labels = 1:4,
    ordered_result = TRUE
  )
}

cluster <- rep(seq_len(n_clusters), each = cluster_size)
stratum_by_cluster <- rep(seq_len(6), each = n_clusters / 6)
stratum <- stratum_by_cluster[cluster]
base_weight <- c(0.75, 0.90, 1.05, 1.20, 1.35, 1.50)[stratum]
weight <- round(base_weight * runif(n, min = 0.85, max = 1.15), 6)

dat$wgt <- weight
dat$clu <- cluster
dat$str <- stratum

mplus_dat <- dat
for (item in items) {
  mplus_dat[[item]] <- as.integer(as.character(mplus_dat[[item]]))
}

utils::write.table(
  mplus_dat[, c(items, "wgt", "clu", "str")],
  file = file.path(validation_dir, "ordinal_survey_complex.dat"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

writeLines(c(
  "TITLE: Demo-compatible ordinal survey CFA for lavaan.survey validation;",
  "DATA: FILE IS ordinal_survey_complex.dat;",
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
), con = file.path(validation_dir, "ordinal_survey_complex.inp"))

cfa_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
"

des <- survey::svydesign(
  ids = ~clu,
  strata = ~str,
  weights = ~wgt,
  data = dat,
  nest = TRUE
)

fit_naive <- lavaan::cfa(
  model = cfa_model,
  data = dat,
  ordered = items,
  estimator = "WLSMV"
)

fit_survey <- lavaan.survey.ordinal(
  lavaan.fit = fit_naive,
  survey.design = des,
  estimator = "WLSMV"
)

pe <- lavaan::parameterEstimates(fit_survey, standardized = TRUE)
pe <- pe[pe$op %in% c("=~", "|", "~~"), ]

utils::write.csv(
  pe[, c("lhs", "op", "rhs", "est", "se", "z", "pvalue", "std.all")],
  file = file.path(validation_dir, "lavaan_survey_complex_parameters.csv"),
  row.names = FALSE
)

fm <- lavaan::fitMeasures(
  fit_survey,
  c("chisq", "df", "pvalue",
    "chisq.scaled", "df.scaled", "pvalue.scaled",
    "baseline.chisq.scaled", "baseline.df.scaled",
    "cfi.scaled", "tli.scaled", "rmsea.scaled",
    "cfi.robust", "tli.robust", "rmsea.robust", "srmr")
)

utils::write.csv(
  data.frame(measure = names(fm), value = unname(fm)),
  file = file.path(validation_dir, "lavaan_survey_complex_fit.csv"),
  row.names = FALSE
)

message("Validation files written to: ", normalizePath(validation_dir))
message("Next step: run `mpdemo ordinal_survey_complex.inp` inside that folder.")
