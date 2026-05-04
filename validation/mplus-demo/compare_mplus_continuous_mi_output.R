validation_dir <- file.path("validation", "mplus-demo")
if (!dir.exists(validation_dir)) {
  stop("Run this script from the package root.")
}

out_file <- file.path(validation_dir, "continuous_mi_complex.out")
if (!file.exists(out_file)) {
  stop(
    "Mplus output not found. Run `/Applications/MplusDemo/mpdemo ",
    "continuous_mi_complex.inp` inside validation/mplus-demo first."
  )
}

if (!requireNamespace("MplusAutomation", quietly = TRUE)) {
  stop("Install MplusAutomation first: install.packages('MplusAutomation')")
}

model <- MplusAutomation::readModels(out_file, quiet = TRUE)

pick_one <- function(data, candidates) {
  hit <- intersect(candidates, names(data))
  if (length(hit) == 0L) return(NA_real_)
  data[[hit[[1L]]]][1]
}

if (!is.null(model$summaries)) {
  mplus_summary <- as.data.frame(model$summaries)
  utils::write.csv(
    mplus_summary,
    file = file.path(validation_dir, "continuous_mi_mplus_fit_summary.csv"),
    row.names = FALSE
  )

  lavaan_fit <- utils::read.csv(
    file.path(validation_dir, "continuous_mi_lavaan_survey_fit.csv"),
    stringsAsFactors = FALSE
  )
  get_lavaan_fit <- function(measure) {
    value <- lavaan_fit$value[lavaan_fit$measure == measure]
    if (length(value) == 0L) NA_real_ else value[[1L]]
  }

  fit_comparison <- data.frame(
    measure = c("chi_square", "df", "p_value", "cfi", "tli", "rmsea", "srmr"),
    lavaan_scaled = c(
      get_lavaan_fit("chisq.scaled"),
      get_lavaan_fit("df.scaled"),
      get_lavaan_fit("pvalue.scaled"),
      get_lavaan_fit("cfi.scaled"),
      get_lavaan_fit("tli.scaled"),
      get_lavaan_fit("rmsea.scaled"),
      get_lavaan_fit("srmr")
    ),
    mplus_mlr = c(
      pick_one(
        mplus_summary,
        c("ChiSqM_Value", "ChiSqM_Mean", "ChiSqM_Value_Mean", "ChiSqM_Value_MLR")
      ),
      pick_one(mplus_summary, c("ChiSqM_DF", "ChiSqM_DF_Mean")),
      pick_one(mplus_summary, c("ChiSqM_PValue", "ChiSqM_PValue_Mean")),
      pick_one(mplus_summary, c("CFI", "CFI_Mean")),
      pick_one(mplus_summary, c("TLI", "TLI_Mean")),
      pick_one(
        mplus_summary,
        c("RMSEA_Estimate", "RMSEA_Mean", "RMSEA_Estimate_Mean")
      ),
      pick_one(mplus_summary, c("SRMR", "SRMR_Mean"))
    )
  )

  utils::write.csv(
    fit_comparison,
    file = file.path(validation_dir, "continuous_mi_mplus_lavaan_fit_comparison.csv"),
    row.names = FALSE
  )
}

params <- model$parameters$unstandardized
if (is.null(params)) {
  stop("MplusAutomation did not find unstandardized parameter estimates.")
}

utils::write.csv(
  params,
  file = file.path(validation_dir, "continuous_mi_mplus_parameters_raw.csv"),
  row.names = FALSE
)

lavaan_params <- utils::read.csv(
  file.path(validation_dir, "continuous_mi_lavaan_survey_parameters.csv"),
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
    } else if (grepl("\\.WITH$", header)) {
      lhs <- sub("\\.WITH$", "", header)
      out <- rbind(out, add(lhs, "~~", par, est, se))
    } else if (identical(header, "Variances")) {
      out <- rbind(out, add(par, "~~", par, est, se))
    } else if (identical(header, "Residual.Variances") ||
               identical(header, "Residual Variances")) {
      out <- rbind(out, add(par, "~~", par, est, se))
    } else if (identical(header, "Intercepts")) {
      out <- rbind(out, add(par, "~1", "", est, se))
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
  file = file.path(validation_dir, "continuous_mi_mplus_lavaan_parameter_comparison.csv"),
  row.names = FALSE
)

message("Continuous MI comparison written to: ",
        normalizePath(file.path(
          validation_dir,
          "continuous_mi_mplus_lavaan_parameter_comparison.csv"
        )))
message("Continuous MI fit comparison written to: ",
        normalizePath(file.path(
          validation_dir,
          "continuous_mi_mplus_lavaan_fit_comparison.csv"
        )))
