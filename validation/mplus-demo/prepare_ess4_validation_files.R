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
  if (requireNamespace("lavaan.survey", quietly = TRUE) &&
      exists("lavaan.survey.ordinal", envir = asNamespace("lavaan.survey"))) {
    library(lavaan.survey)
  } else {
    source(file.path("R", "lavaan.survey.R"))
  }
}

load_lavaan_survey()

data_path <- file.path("data", "ess4.gb.rda")
if (!file.exists(data_path)) {
  stop("Could not find ", data_path)
}

load(data_path)

items <- c("gvjbevn", "gvhlthc", "gvslvol", "gvslvue", "gvcldcr", "gvpdlwk")
design_vars <- c("psu", "stratval", "dweight")
keep <- stats::complete.cases(ess4.gb[, c(items, design_vars)])
dat <- ess4.gb[keep, c(items, design_vars)]

# Collapse 0-10 ESS response scales into four ordered categories. This keeps
# the original variables and one-factor structure, but avoids very sparse
# thresholds in replicate-weight polychoric correlations.
collapse_ess <- function(x) {
  cut(
    x,
    breaks = c(-Inf, 4, 6, 8, Inf),
    labels = 1:4,
    ordered_result = TRUE
  )
}

for (item in items) {
  dat[[item]] <- collapse_ess(dat[[item]])
}

mplus_dat <- dat
for (item in items) {
  mplus_dat[[item]] <- as.integer(as.character(mplus_dat[[item]]))
}
mplus_dat$psu <- as.integer(droplevels(mplus_dat$psu))
mplus_dat$strat <- as.integer(droplevels(mplus_dat$stratval))
mplus_dat$stratval <- NULL

utils::write.table(
  mplus_dat[, c(items, "dweight", "psu", "strat")],
  file = file.path(validation_dir, "ess4_range_complex.dat"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

writeLines(c(
  "TITLE: ESS4 GB ordinal survey CFA for lavaan.survey validation;",
  "DATA: FILE IS ess4_range_complex.dat;",
  "VARIABLE:",
  "  NAMES ARE gvjbevn gvhlthc gvslvol gvslvue gvcldcr gvpdlwk dweight psu strat;",
  "  USEVARIABLES ARE gvjbevn gvhlthc gvslvol gvslvue gvcldcr gvpdlwk;",
  "  CATEGORICAL ARE gvjbevn gvhlthc gvslvol gvslvue gvcldcr gvpdlwk;",
  "  WEIGHT IS dweight;",
  "  CLUSTER IS psu;",
  "  STRATIFICATION IS strat;",
  "ANALYSIS:",
  "  TYPE = COMPLEX;",
  "  ESTIMATOR = WLSMV;",
  "  PARAMETERIZATION = DELTA;",
  "MODEL:",
  "  range BY gvjbevn gvhlthc gvslvol gvslvue gvcldcr gvpdlwk;",
  "OUTPUT:",
  "  STDYX TECH1;"
), con = file.path(validation_dir, "ess4_range_complex.inp"))

cfa_model <- "
  range =~ gvjbevn + gvhlthc + gvslvol + gvslvue + gvcldcr + gvpdlwk
"

options(survey.lonely.psu = "adjust")

des <- survey::svydesign(
  ids = ~psu,
  strata = ~stratval,
  weights = ~dweight,
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
  file = file.path(validation_dir, "ess4_lavaan_survey_parameters.csv"),
  row.names = FALSE
)

fm <- lavaan::fitMeasures(
  fit_survey,
  c("chisq.scaled", "df.scaled", "pvalue.scaled",
    "cfi.robust", "tli.robust", "rmsea.robust", "srmr")
)

utils::write.csv(
  data.frame(measure = names(fm), value = unname(fm)),
  file = file.path(validation_dir, "ess4_lavaan_survey_fit.csv"),
  row.names = FALSE
)

category_counts <- as.data.frame.matrix(sapply(dat[items], table))
utils::write.csv(
  category_counts,
  file = file.path(validation_dir, "ess4_collapsed_category_counts.csv")
)

message("ESS4 validation files written to: ", normalizePath(validation_dir))
message("Next step: run `/Applications/MplusDemo/mpdemo ess4_range_complex.inp` inside that folder.")

