validation_dir <- file.path("validation", "mplus-demo")
if (!dir.exists(validation_dir)) {
  stop("Run this script from the package root.")
}

out_file <- file.path(validation_dir, "mixed_group_mi_complex.out")
if (!file.exists(out_file)) {
  stop(
    "Mplus output not found. Run `/Applications/MplusDemo/mpdemo ",
    "mixed_group_mi_complex.inp` inside validation/mplus-demo first."
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
    file = file.path(validation_dir, "mixed_group_mi_mplus_fit_summary.csv"),
    row.names = FALSE
  )

  read_lavaan_fit <- function(file) {
    path <- file.path(validation_dir, file)
    if (!file.exists(path)) {
      return(data.frame(measure = character(), value = numeric()))
    }
    utils::read.csv(path, stringsAsFactors = FALSE)
  }
  get_lavaan_fit <- function(lavaan_fit, measure) {
    value <- lavaan_fit$value[lavaan_fit$measure == measure]
    if (length(value) == 0L) NA_real_ else value[[1L]]
  }
  lavaan_fit_stats <- read_lavaan_fit("mixed_group_mi_lavaan_survey_fit.csv")
  lavaan_fit_params <- read_lavaan_fit(
    "mixed_group_mi_lavaan_survey_fit_parameter_pooling.csv"
  )

  fit_comparison <- data.frame(
    measure = c("chi_square", "df", "p_value", "cfi", "tli", "rmsea", "srmr"),
    lavaan_sample_stat_scaled = c(
      get_lavaan_fit(lavaan_fit_stats, "chisq.scaled"),
      get_lavaan_fit(lavaan_fit_stats, "df.scaled"),
      get_lavaan_fit(lavaan_fit_stats, "pvalue.scaled"),
      get_lavaan_fit(lavaan_fit_stats, "cfi.scaled"),
      get_lavaan_fit(lavaan_fit_stats, "tli.scaled"),
      get_lavaan_fit(lavaan_fit_stats, "rmsea.scaled"),
      get_lavaan_fit(lavaan_fit_stats, "srmr")
    ),
    lavaan_sample_stat_robust = c(
      NA_real_,
      NA_real_,
      NA_real_,
      get_lavaan_fit(lavaan_fit_stats, "cfi.robust"),
      get_lavaan_fit(lavaan_fit_stats, "tli.robust"),
      get_lavaan_fit(lavaan_fit_stats, "rmsea.robust"),
      get_lavaan_fit(lavaan_fit_stats, "srmr")
    ),
    lavaan_parameter_pooling_mean = c(
      get_lavaan_fit(lavaan_fit_params, "chisq.scaled"),
      get_lavaan_fit(lavaan_fit_params, "df.scaled"),
      NA_real_,
      get_lavaan_fit(lavaan_fit_params, "cfi.scaled"),
      NA_real_,
      get_lavaan_fit(lavaan_fit_params, "rmsea.scaled"),
      get_lavaan_fit(lavaan_fit_params, "srmr")
    ),
    mplus_wlsmv = c(
      pick_one(
        mplus_summary,
        c("ChiSqM_Value", "ChiSqM_Mean", "ChiSqM_Value_Mean")
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
    file = file.path(validation_dir,
                     "mixed_group_mi_mplus_lavaan_fit_comparison.csv"),
    row.names = FALSE
  )
}

params <- model$parameters$unstandardized
if (is.null(params)) {
  stop("MplusAutomation did not find unstandardized parameter estimates.")
}

utils::write.csv(
  params,
  file = file.path(validation_dir, "mixed_group_mi_mplus_parameters_raw.csv"),
  row.names = FALSE
)

normalize_group <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "", x)
  ifelse(x %in% c("boy", "boys"), "boys",
         ifelse(x %in% c("girl", "girls"), "girls", x))
}

detect_group <- function(row) {
  candidates <- c("Group", "group", "BetweenWithin", "LatentClass")
  hit <- intersect(candidates, names(row))
  for (name in hit) {
    value <- row[[name]][1]
    if (!is.na(value) && nzchar(as.character(value))) {
      group <- normalize_group(as.character(value))
      if (group %in% c("boys", "girls")) return(group)
    }
  }

  header <- as.character(row$paramHeader[1])
  header_group <- sub("^.*\\.(BOYS|GIRLS)$", "\\1", header, ignore.case = TRUE)
  header_group <- normalize_group(header_group)
  if (header_group %in% c("boys", "girls")) return(header_group)

  NA_character_
}

normalize_mplus <- function(params) {
  out <- data.frame(
    group_label = character(),
    lhs = character(),
    op = character(),
    rhs = character(),
    mplus_est = numeric(),
    mplus_se = numeric(),
    stringsAsFactors = FALSE
  )

  add <- function(group_label, lhs, op, rhs, est, se) {
    data.frame(
      group_label = group_label,
      lhs = tolower(lhs),
      op = op,
      rhs = tolower(rhs),
      mplus_est = est,
      mplus_se = se,
      stringsAsFactors = FALSE
    )
  }

  current_group <- NA_character_
  for (i in seq_len(nrow(params))) {
    row <- params[i, , drop = FALSE]
    group <- detect_group(row)
    if (!is.na(group)) {
      current_group <- group
    }
    if (is.na(current_group)) {
      current_group <- "boys"
    }

    header <- as.character(params$paramHeader[i])
    header_clean <- sub("\\.(BOYS|GIRLS)$", "", header, ignore.case = TRUE)
    par <- as.character(params$param[i])
    if (identical(header_clean, "Variances") && par == "FALSE") {
      par <- "F"
    }
    est <- params$est[i]
    se <- params$se[i]

    if (grepl("\\.BY$", header_clean)) {
      lhs <- sub("\\.BY$", "", header_clean)
      out <- rbind(out, add(current_group, lhs, "=~", par, est, se))
    } else if (identical(header_clean, "Thresholds") && grepl("\\$", par)) {
      bits <- strsplit(par, "\\$")[[1]]
      out <- rbind(out, add(current_group, bits[1], "|", paste0("t", bits[2]),
                            est, se))
    } else if (grepl("\\.WITH$", header_clean)) {
      lhs <- sub("\\.WITH$", "", header_clean)
      out <- rbind(out, add(current_group, lhs, "~~", par, est, se))
    } else if (identical(header_clean, "Variances")) {
      out <- rbind(out, add(current_group, par, "~~", par, est, se))
    } else if (identical(header_clean, "Intercepts") ||
               identical(header_clean, "Means")) {
      out <- rbind(out, add(current_group, par, "~1", "", est, se))
    }
  }

  out$key <- with(out, paste(group_label, lhs, op, rhs, sep = "|"))
  out
}

mplus_params <- normalize_mplus(params)

compare_lavaan_parameters <- function(input_file, output_file, algorithm) {
  path <- file.path(validation_dir, input_file)
  if (!file.exists(path)) {
    warning("lavaan parameter file not found: ", path)
    return(invisible(NULL))
  }

  lavaan_params <- utils::read.csv(path, stringsAsFactors = FALSE)
  lavaan_params$key <- with(
    lavaan_params,
    paste(tolower(group_label), tolower(lhs), op, tolower(rhs), sep = "|")
  )

  comparison <- merge(
    lavaan_params,
    mplus_params[, c("key", "mplus_est", "mplus_se")],
    by = "key",
    all.x = TRUE
  )
  comparison$algorithm <- algorithm
  comparison$est_diff_mplus_minus_lavaan <- comparison$mplus_est - comparison$est
  comparison$se_diff_mplus_minus_lavaan <- comparison$mplus_se - comparison$se

  utils::write.csv(
    comparison,
    file = file.path(validation_dir, output_file),
    row.names = FALSE
  )
  invisible(comparison)
}

comparison_sample_stats <- compare_lavaan_parameters(
  input_file = "mixed_group_mi_lavaan_survey_parameters.csv",
  output_file = "mixed_group_mi_mplus_lavaan_parameter_comparison.csv",
  algorithm = "sample_statistics"
)
comparison_parameter_pooling <- compare_lavaan_parameters(
  input_file = "mixed_group_mi_lavaan_survey_parameters_parameter_pooling.csv",
  output_file = "mixed_group_mi_mplus_lavaan_parameter_comparison_parameter_pooling.csv",
  algorithm = "parameter_pooling"
)

combined <- do.call(rbind, Filter(Negate(is.null), list(
  comparison_sample_stats,
  comparison_parameter_pooling
)))
if (!is.null(combined) && nrow(combined) > 0L) {
  utils::write.csv(
    combined,
    file = file.path(validation_dir,
                     "mixed_group_mi_mplus_lavaan_parameter_comparison_all_algorithms.csv"),
    row.names = FALSE
  )
}

message("Mixed multiple-group MI parameter comparison written to: ",
        normalizePath(file.path(
          validation_dir,
          "mixed_group_mi_mplus_lavaan_parameter_comparison.csv"
        )))
message("Mixed multiple-group MI parameter-pooling comparison written to: ",
        normalizePath(file.path(
          validation_dir,
          "mixed_group_mi_mplus_lavaan_parameter_comparison_parameter_pooling.csv"
        )))
message("Mixed multiple-group MI fit comparison written to: ",
        normalizePath(file.path(
          validation_dir,
          "mixed_group_mi_mplus_lavaan_fit_comparison.csv"
        )))
