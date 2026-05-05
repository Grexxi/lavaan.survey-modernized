validation_dir <- file.path("validation", "mplus-demo")
if (!dir.exists(validation_dir)) {
  stop("Run this script from the package root.")
}

required <- c("lavaan", "survey", "mitools", "mice")
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

to_ordered3 <- function(x, breaks) {
  ordered(cut(x, breaks = breaks, labels = 1:3), levels = 1:3)
}

to_mplus_integer_data <- function(data, items) {
  out <- data
  for (item in items) {
    out[[item]] <- as.integer(as.character(out[[item]]))
  }
  out
}

load_lavaan_survey()

set.seed(20260518)

n_clusters <- 48L
cluster_size <- 12L
n <- n_clusters * cluster_size
m <- 10L
items_ord <- c("y1", "y2")
items_cont <- c("x1", "x2")
items <- c(items_ord, items_cont)

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
lambda_boys <- c(0.86, 0.78, 0.82, 0.70)
lambda_girls <- c(0.88, 0.76, 0.84, 0.72)
lambda <- rbind(lambda_boys, lambda_girls)

y1_star <- lambda[group_id, 1] * eta +
  sqrt(1 - lambda[group_id, 1]^2) * rnorm(n)
y2_star <- lambda[group_id, 2] * eta +
  sqrt(1 - lambda[group_id, 2]^2) * rnorm(n)
x1 <- 0.25 + lambda[group_id, 3] * eta +
  sqrt(1 - lambda[group_id, 3]^2) * rnorm(n)
x2 <- -0.20 + lambda[group_id, 4] * eta +
  sqrt(1 - lambda[group_id, 4]^2) * rnorm(n)

dat <- data.frame(
  y1 = to_ordered3(y1_star, breaks = c(-Inf, -0.65, 0.25, Inf)),
  y2 = to_ordered3(y2_star, breaks = c(-Inf, -0.45, 0.55, Inf)),
  x1 = x1,
  x2 = x2,
  grp = group,
  grp_id = group_id,
  wgt = weight,
  clu = cluster,
  str = stratum
)

# Create moderate missingness in one ordered and one continuous indicator. The
# imputation model includes group and design information; Mplus and
# lavaan.survey.ordinal then receive the same completed datasets.
dat_missing <- dat
p_y2 <- plogis(
  -1.35 + 0.45 * scale(y1_star)[, 1] + 0.25 * (group_id == 2L)
)
p_x2 <- plogis(
  -1.45 + 0.35 * scale(x1)[, 1] - 0.20 * scale(weight)[, 1]
)
dat_missing$y2[runif(n) < p_y2] <- NA
dat_missing$x2[runif(n) < p_x2] <- NA

ini <- mice::mice(dat_missing, maxit = 0, printFlag = FALSE)
method <- ini$method
method[] <- ""
method["y2"] <- "polr"
method["x2"] <- "pmm"

predictor_matrix <- ini$predictorMatrix
predictor_matrix[,] <- 0
predictor_matrix[c("y2", "x2"), c(items, "grp_id", "wgt", "str")] <- 1
diag(predictor_matrix) <- 0

imp <- mice::mice(
  dat_missing,
  m = m,
  maxit = 10,
  method = method,
  predictorMatrix = predictor_matrix,
  seed = 20260519,
  printFlag = FALSE
)

imputed_data <- lapply(seq_len(m), function(i) {
  out <- mice::complete(imp, action = i)
  for (item in items_ord) {
    out[[item]] <- ordered(out[[item]], levels = 1:3)
  }
  out$grp <- factor(out$grp, levels = c("boys", "girls"))
  out$grp_id <- as.integer(out$grp_id)
  out
})

imputation_files <- sprintf("mixed_group_mi_imp%02d.dat", seq_len(m))
for (i in seq_along(imputed_data)) {
  mplus_data <- to_mplus_integer_data(imputed_data[[i]], items_ord)
  utils::write.table(
    mplus_data[, c(items, "grp_id", "wgt", "clu", "str")],
    file = file.path(validation_dir, imputation_files[[i]]),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

writeLines(
  imputation_files,
  con = file.path(validation_dir, "mixed_group_mi_implist.dat")
)

writeLines(c(
  "TITLE: Mixed group MI complex CFA validation;",
  "DATA:",
  "  FILE IS mixed_group_mi_implist.dat;",
  "  TYPE = IMPUTATION;",
  "VARIABLE:",
  "  NAMES ARE y1 y2 x1 x2 grp wgt clu str;",
  "  USEVARIABLES ARE y1 y2 x1 x2;",
  "  CATEGORICAL ARE y1 y2;",
  "  GROUPING IS grp (1 = boys 2 = girls);",
  "  WEIGHT IS wgt;",
  "  CLUSTER IS clu;",
  "  STRATIFICATION IS str;",
  "ANALYSIS:",
  "  TYPE = COMPLEX;",
  "  ESTIMATOR = WLSMV;",
  "  PARAMETERIZATION = DELTA;",
  "MODEL:",
  "  f BY y1 y2 x1 x2;",
  "OUTPUT:",
  "  STDYX TECH1;"
), con = file.path(validation_dir, "mixed_group_mi_complex.inp"))

cfa_model <- "
  f =~ y1 + y2 + x1 + x2
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
  ordered = items_ord,
  group = "grp",
  group.equal = c("loadings", "thresholds", "intercepts"),
  estimator = "WLSMV"
)

fit_survey <- lavaan.survey.ordinal(
  lavaan.fit = fit_naive,
  survey.design = design,
  estimator = "WLSMV"
)

fit_survey_parameter_pooling <- lavaan.survey.ordinal(
  lavaan.fit = fit_naive,
  survey.design = design,
  estimator = "WLSMV",
  point.wls = "lavaan",
  mi.pooling = "parameters"
)

keep_parameter_rows <- function(pe) {
  pe[
    pe$op %in% c("=~", "|") |
      (pe$op == "~1" & pe$lhs %in% items_cont) |
      (pe$op == "~~" & pe$lhs == "f" & pe$rhs == "f"),
  ]
}

write_lavaan_outputs <- function(fit, parameter_file, fit_file) {
  pe <- lavaan::parameterEstimates(fit, standardized = TRUE)
  pe <- keep_parameter_rows(pe)
  group_labels <- lavaan::lavInspect(fit, "group.label")
  pe$group_label <- group_labels[pe$group]

  utils::write.csv(
    pe[, c("group", "group_label", "lhs", "op", "rhs", "est", "se", "z",
           "pvalue", "std.all")],
    file = file.path(validation_dir, parameter_file),
    row.names = FALSE
  )

  fm <- lavaan::fitMeasures(
    fit,
    c("chisq", "df", "pvalue",
      "chisq.scaled", "df.scaled", "pvalue.scaled",
      "cfi.scaled", "tli.scaled", "rmsea.scaled",
      "cfi.robust", "tli.robust", "rmsea.robust", "srmr")
  )
  utils::write.csv(
    data.frame(measure = names(fm), value = unname(fm)),
    file = file.path(validation_dir, fit_file),
    row.names = FALSE
  )
}

write_parameter_pooled_outputs <- function(fit, parameter_file, fit_file) {
  template <- fit$fits[[1]]
  pt <- lavaan::parTable(template)
  free_rows <- pt[pt$free > 0 & pt$op != "==", ]
  free_rows <- free_rows[order(free_rows$free), ]

  coef_names <- names(coef(fit))
  if (nrow(free_rows) != length(coef_names)) {
    stop("Could not map pooled free parameters back to lavaan parameter table.")
  }

  se <- sqrt(diag(vcov(fit)))
  z <- coef(fit) / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  group_labels <- lavaan::lavInspect(template, "group.label")

  pe <- data.frame(
    group = free_rows$group,
    group_label = group_labels[free_rows$group],
    lhs = free_rows$lhs,
    op = free_rows$op,
    rhs = free_rows$rhs,
    est = unname(coef(fit)),
    se = unname(se),
    z = unname(z),
    pvalue = unname(p),
    std.all = NA_real_,
    stringsAsFactors = FALSE
  )
  pe <- keep_parameter_rows(pe)

  utils::write.csv(
    pe,
    file = file.path(validation_dir, parameter_file),
    row.names = FALSE
  )

  utils::write.csv(
    data.frame(measure = names(fit$fit.measures),
               value = unname(fit$fit.measures)),
    file = file.path(validation_dir, fit_file),
    row.names = FALSE
  )
}

write_lavaan_outputs(
  fit_survey,
  parameter_file = "mixed_group_mi_lavaan_survey_parameters.csv",
  fit_file = "mixed_group_mi_lavaan_survey_fit.csv"
)
write_parameter_pooled_outputs(
  fit_survey_parameter_pooling,
  parameter_file = "mixed_group_mi_lavaan_survey_parameters_parameter_pooling.csv",
  fit_file = "mixed_group_mi_lavaan_survey_fit_parameter_pooling.csv"
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
  file = file.path(validation_dir, "mixed_group_mi_missing_summary.csv"),
  row.names = FALSE
)

message("Mixed multiple-group MI validation files written to: ",
        normalizePath(validation_dir))
message(
  "Next step: run `/Applications/MplusDemo/mpdemo ",
  "mixed_group_mi_complex.inp` inside that folder."
)
