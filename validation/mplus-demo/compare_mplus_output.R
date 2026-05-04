validation_dir <- file.path("validation", "mplus-demo")
if (!dir.exists(validation_dir)) {
  stop("Run this script from the package root.")
}

out_file <- file.path(validation_dir, "ordinal_survey_complex.out")
if (!file.exists(out_file)) {
  stop(
    "Mplus output not found. Run `mpdemo ordinal_survey_complex.inp` ",
    "inside validation/mplus-demo first."
  )
}

if (!requireNamespace("MplusAutomation", quietly = TRUE)) {
  stop("Install MplusAutomation first: install.packages('MplusAutomation')")
}

model <- MplusAutomation::readModels(out_file, quiet = TRUE)

if (!is.null(model$summaries)) {
  utils::write.csv(
    as.data.frame(model$summaries),
    file = file.path(validation_dir, "mplus_complex_fit_summary.csv"),
    row.names = FALSE
  )
}

params <- model$parameters$unstandardized
if (is.null(params)) {
  stop("MplusAutomation did not find unstandardized parameter estimates.")
}

utils::write.csv(
  params,
  file = file.path(validation_dir, "mplus_complex_parameters_raw.csv"),
  row.names = FALSE
)

lavaan_params <- utils::read.csv(
  file.path(validation_dir, "lavaan_survey_complex_parameters.csv"),
  stringsAsFactors = FALSE
)

lavaan_params$key <- with(
  lavaan_params,
  paste(tolower(lhs), op, tolower(rhs), sep = "|")
)

normalize_mplus <- function(params) {
  out <- data.frame(
    lhs = character(),
    op = character(),
    rhs = character(),
    mplus_est = numeric(),
    mplus_se = numeric(),
    stringsAsFactors = FALSE
  )

  add <- function(lhs, op, rhs, est, se) {
    data.frame(
      lhs = tolower(lhs),
      op = op,
      rhs = tolower(rhs),
      mplus_est = est,
      mplus_se = se,
      stringsAsFactors = FALSE
    )
  }

  for (i in seq_len(nrow(params))) {
    header <- params$paramHeader[i]
    par <- params$param[i]
    est <- params$est[i]
    se <- params$se[i]

    if (grepl("\\.BY$", header)) {
      lhs <- sub("\\.BY$", "", header)
      out <- rbind(out, add(lhs, "=~", par, est, se))
    } else if (identical(header, "Thresholds") && grepl("\\$", par)) {
      bits <- strsplit(par, "\\$")[[1]]
      out <- rbind(out, add(bits[1], "|", paste0("t", bits[2]), est, se))
    } else if (grepl("\\.WITH$", header)) {
      lhs <- sub("\\.WITH$", "", header)
      out <- rbind(out, add(lhs, "~~", par, est, se))
    } else if (identical(header, "Variances")) {
      out <- rbind(out, add(par, "~~", par, est, se))
    }
  }

  out$key <- with(out, paste(lhs, op, rhs, sep = "|"))
  out
}

mplus_params <- normalize_mplus(params)

comparison <- merge(
  lavaan_params,
  mplus_params[, c("key", "mplus_est", "mplus_se")],
  by = "key",
  all.x = TRUE
)

comparison$est_diff_mplus_minus_lavaan <- comparison$mplus_est - comparison$est
comparison$se_diff_mplus_minus_lavaan <- comparison$mplus_se - comparison$se

utils::write.csv(
  comparison,
  file = file.path(validation_dir, "mplus_lavaan_parameter_comparison.csv"),
  row.names = FALSE
)

message("Comparison written to: ",
        normalizePath(file.path(validation_dir, "mplus_lavaan_parameter_comparison.csv")))

