validation_dir <- file.path("validation", "mplus-demo")
if (!dir.exists(validation_dir)) {
  stop("Run this script from the package root.")
}

required <- c("lavaan", "survey")
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
    cut(x, breaks = c(-Inf, -0.7, 0.1, 0.9, Inf), labels = 1:4),
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

set.seed(20260510)

n_clusters <- 48L
cluster_size <- 12L
n <- n_clusters * cluster_size
items <- paste0("y", 1:4)

cluster <- rep(seq_len(n_clusters), each = cluster_size)
group_by_cluster <- rep(c(1L, 2L), each = n_clusters / 2L)
group_id <- group_by_cluster[cluster]
group <- factor(group_id, levels = c(1L, 2L), labels = c("boys", "girls"))

stratum_by_cluster <- rep(rep(seq_len(4), each = 6L), times = 2L)
stratum <- stratum_by_cluster[cluster]
base_weight <- c(0.80, 0.95, 1.15, 1.35)[stratum]
group_weight <- ifelse(group_id == 1L, 1.05, 0.95)
weight <- round(base_weight * group_weight * runif(n, min = 0.85, max = 1.15), 6)

eta <- rnorm(n, mean = ifelse(group_id == 1L, -0.15, 0.20), sd = 1)
lambda_boys <- c(0.82, 0.76, 0.71, 0.66)
lambda_girls <- c(0.86, 0.73, 0.74, 0.69)
lambda <- rbind(lambda_boys, lambda_girls)

y_star <- matrix(NA_real_, nrow = n, ncol = length(items))
for (j in seq_along(items)) {
  loading <- lambda[group_id, j]
  y_star[, j] <- loading * eta + sqrt(1 - loading^2) * rnorm(n)
}

dat <- as.data.frame(y_star)
names(dat) <- items
for (item in items) {
  dat[[item]] <- to_ordered4(dat[[item]])
}
dat$grp <- group
dat$grp_id <- group_id
dat$wgt <- weight
dat$clu <- cluster
dat$str <- stratum

mplus_dat <- to_mplus_integer_data(dat, items)
utils::write.table(
  mplus_dat[, c(items, "grp_id", "wgt", "clu", "str")],
  file = file.path(validation_dir, "ordinal_group_complex.dat"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

writeLines(c(
  "TITLE: Multiple-group ordinal complex CFA for lavaan.survey validation;",
  "DATA: FILE IS ordinal_group_complex.dat;",
  "VARIABLE:",
  "  NAMES ARE y1 y2 y3 y4 grp wgt clu str;",
  "  USEVARIABLES ARE y1 y2 y3 y4;",
  "  CATEGORICAL ARE y1 y2 y3 y4;",
  "  GROUPING IS grp (1 = boys 2 = girls);",
  "  WEIGHT IS wgt;",
  "  CLUSTER IS clu;",
  "  STRATIFICATION IS str;",
  "ANALYSIS:",
  "  TYPE = COMPLEX;",
  "  ESTIMATOR = WLSMV;",
  "  PARAMETERIZATION = DELTA;",
  "MODEL:",
  "  f BY y1 y2 y3 y4;",
  "OUTPUT:",
  "  STDYX TECH1;"
), con = file.path(validation_dir, "ordinal_group_complex.inp"))

cfa_model <- "
  f =~ y1 + y2 + y3 + y4
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
  group = "grp",
  group.equal = c("loadings", "thresholds", "intercepts"),
  estimator = "WLSMV"
)

fit_survey <- lavaan.survey.ordinal(
  lavaan.fit = fit_naive,
  survey.design = des,
  estimator = "WLSMV"
)

pe <- lavaan::parameterEstimates(fit_survey, standardized = TRUE)
pe <- pe[
  pe$op %in% c("=~", "|") |
    (pe$op == "~~" & pe$lhs == "f" & pe$rhs == "f"),
]
group_labels <- lavaan::lavInspect(fit_survey, "group.label")
pe$group_label <- group_labels[pe$group]

utils::write.csv(
  pe[, c("group", "group_label", "lhs", "op", "rhs", "est", "se", "z",
         "pvalue", "std.all")],
  file = file.path(validation_dir, "ordinal_group_lavaan_survey_parameters.csv"),
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
  file = file.path(validation_dir, "ordinal_group_lavaan_survey_fit.csv"),
  row.names = FALSE
)

message("Ordinal multiple-group validation files written to: ",
        normalizePath(validation_dir))
message(
  "Next step: run `/Applications/MplusDemo/mpdemo ",
  "ordinal_group_complex.inp` inside that folder."
)
