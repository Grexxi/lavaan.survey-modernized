# Complex sampling analysis of SEM models
# Daniel Oberski, 2015-11-03

lavaan.survey <- 
  function(lavaan.fit, survey.design, 
           estimator=c("MLM", "MLMV", "MLMVS", "WLS", "DWLS", "ML"),
           estimator.gamma=c("default","Yuan-Bentler"),
           ordered=NULL, verbose=interactive(), ...) {

  verbose <- isTRUE(verbose)
  estimator.supplied <- !missing(estimator)
  ordered.from.fit <- lavNames(lavaan.fit, type="ov.ord", group=1)
  ordinal.requested <- isTRUE(ordered) || is.character(ordered) ||
    (is.null(ordered) && length(ordered.from.fit) > 0L)
  if(identical(ordered, FALSE)) ordinal.requested <- FALSE

  if(ordinal.requested) {
    if(isTRUE(ordered)) ordered <- ordered.from.fit
    if(estimator.supplied) {
      estimator <- match.arg(estimator, choices=c("WLSMV", "DWLS"))
    }
    else {
      estimator <- "WLSMV"
    }

    return(lavaan.survey.ordinal(lavaan.fit=lavaan.fit,
                                 survey.design=survey.design,
                                 ordered=ordered,
                                 estimator=estimator,
                                 verbose=verbose, ...))
  }
  
  # Not all estimators in lavaan make sense to use here, therefore matching args
  estimator <- match.arg(estimator) 
  if(estimator=="ML") warning("Estimator 'ML' will not correct standard errors and chi-square statistic.")
  estimator.gamma <- match.arg(estimator.gamma) # Smoothing Gamma or not
  meanstructure <- isTRUE(lavInspect(lavaan.fit, "options")$meanstructure)
  
  # Names of the observed variables (same for each group)
  ov.names <- lavNames(lavaan.fit, type="ov", group=1)
  
  # The MP-inverse duplication matrix is handy for removing redundancy
  Dplus <- lavaan::lav_matrix_duplication_ginv( length(ov.names) )
  # Create a formula that includes all observed variables for svymean
  ov.formula <- as.formula(paste("~", paste(ov.names, collapse="+")))
  
  # <no. group>-sized lists that will contain the asy. covariance matrix,
  #  and sample covariance matrix and mean vector
  ngroups <- lavInspect(lavaan.fit, "ngroups")
  
  Gamma <- vector("list", ngroups)
  sample.cov <- vector("list", ngroups)
  sample.mean <- vector("list", ngroups)
  sample.nobs <- vector("list", ngroups)
  
  for(g in seq(ngroups)) {
    if(ngroups > 1) {
      # Use survey::subset to create data groups
      survey.design.g <- 
        survey_subset_by_group(survey.design,
                               lavInspect(lavaan.fit, "group"),
                               lavInspect(lavaan.fit, "group.label")[[g]])
    } 
    else { # In case of no groups, just use the original survey design object.
      survey.design.g <- survey.design  
    }
    
    # Function that takes survey design and returns the Gamma & observed moments
    get.stats.design <- function(survey.design.g, sample.nobs.g) {
      sample.cov.g <- as.matrix(svyvar(ov.formula, design=survey.design.g, na.rm=TRUE))  
      # survey package returns the variance matrix of the (co)variances as attr:
      Gamma.cov.g <- attr(sample.cov.g, "var")
      # Remove (co)variances wrt redundant elements of covs; not used by lavaan. 
      Gamma.cov.g <- Dplus %*% Gamma.cov.g %*% t(Dplus)
      
      if(meanstructure) {
        # Same for mean vector
        sample.mean.g <- svymean(ov.formula, design=survey.design.g, na.rm=TRUE)
        Gamma.mean.g <- attr(sample.mean.g, "var")

        # Join asy. variance matrices for means and covariances
        # TODO add offdiag
        Gamma.g <- lavaan::lav_matrix_bdiag(Gamma.mean.g, Gamma.cov.g)
      }
      else {
        sample.mean.g <- NULL
        Gamma.g <- Gamma.cov.g
      }
      
      Gamma.g <- Gamma.g * sample.nobs.g # lavaan wants nobs * Gamma.
      
      # Since the above nonparametric estimate of Gamma can be unstable, Yuan
      # and Bentler suggested a model-smoothed estimate of it, optional here:
      if(estimator.gamma == "Yuan-Bentler") {
        r <- get.residuals(lavaan.fit, group=g,
                           meanstructure=meanstructure) # Iff these asy = 0, all will be well...
        Gamma.g <- Gamma.g + (sample.nobs.g/(sample.nobs.g - 1)) * (r %*% t(r))
      }
      # Get rid of attribute, preventing errors caused by lazy evaluation
      # (This has to be at the end or lazy evaluation mayhem will ensue)
      attr(sample.cov.g, "var") <- NULL
      if(meanstructure) {
        tmp  <- as.vector(sample.mean.g)
        names(tmp) <- names(sample.mean.g)
        sample.mean.g <- tmp
      }
      
      list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
    }
    # The data may be a list of multiply imputed datasets
    if(!inherits(survey.design.g, "svyimputationList")) {
      # If no imputations, just use usual no. observations and asy variance
      sample.nobs.g <- lavInspect(lavaan.fit, "nobs")[[g]] 

      stats <- get.stats.design(survey.design.g, sample.nobs.g)
    } 
    else { # In case of multiply imputed data
      # Not only can nobs differ from lavaan.fit, but also per imputation
      sample.nobs.g <- get.sample.nobs(survey.design.g)
      
      # Retrieve point and variance estimates per imputation.
      stats.list <- lapply(seq_along(survey.design.g$designs), function(i) {
        if(verbose) {
          if(ngroups > 1) {
            message("Computing sample statistics for group ", g, "/", ngroups,
                    ", imputation ", i, "/", length(survey.design.g$designs),
                    ".")
          }
          else {
            message("Computing sample statistics for imputation ", i, "/",
                    length(survey.design.g$designs), ".")
          }
        }
        stats <- get.stats.design(survey.design.g$designs[[i]],
                                  sample.nobs.g=sample.nobs.g)
        if(verbose) {
          if(ngroups > 1) {
            message("Computing sample statistics for group ", g, "/", ngroups,
                    ", imputation ", i, "/", length(survey.design.g$designs),
                    " complete.")
          }
          else {
            message("Computing sample statistics for imputation ", i, "/",
                    length(survey.design.g$designs), " complete.")
          }
        }
        stats
      })
      m  <- length(stats.list) # no. imputation
      
      # Point estimates are average over imputations
      sample.cov.list <- lapply(stats.list, `[[`, 'sample.cov.g')
      sample.cov.g <- Reduce(`+`, sample.cov.list) / m
      cov.df <- do.call(rbind, lapply(sample.cov.list, lavaan::lav_matrix_vech))
      if(meanstructure) {
        sample.mean.list <- lapply(stats.list, `[[`, 'sample.mean.g')
        sample.mean.g <- Reduce(`+`, sample.mean.list) / m
        mean.df <- do.call(rbind, sample.mean.list)
      }
      else {
        sample.mean.g <- NULL
        mean.df <- NULL
      }
      
      # Variance estimates depend on within- and between-imputation variance:
      Gamma.within  <- Reduce(`+`, lapply(stats.list, `[[`, 'Gamma.g')) / m
      if(m > 1L) {
        between.df <- if(meanstructure) cbind(mean.df, cov.df) else cov.df
        Gamma.between <- cov(between.df)
      }
      else {
        Gamma.between <- matrix(0, nrow=nrow(Gamma.within),
                                ncol=ncol(Gamma.within))
      }
      dimnames(Gamma.between) <- dimnames(Gamma.within)
      Gamma.g <- Gamma.within + sample.nobs.g * ((m + 1)/m) * Gamma.between
      
      # set stats with multiple imputation point and variance estimates
      stats <- list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
    }

    # Augment the list for this group
    Gamma[[g]] <- stats$Gamma.g
    sample.cov[[g]] <- stats$sample.cov.g
    sample.mean[g] <- list(stats$sample.mean.g)
    sample.nobs[[g]] <- sample.nobs.g
  } # End of loop over groups

  new.call <- lavInspect(lavaan.fit, "call")
  new.call$data <- NULL                # Remove any data argument
  new.call$sample.cov <- sample.cov    # Set survey covariances
  if(meanstructure) new.call$sample.mean <- sample.mean  # Set survey means
  else new.call$sample.mean <- NULL
  new.call$sample.nobs <- sample.nobs  
  new.call$estimator <- estimator  # Always use Satorra-Bentler or WLS estimator

  if(substr(estimator, 1, 2) == "ML") { # ML, robust options
    # Set asymptotic covariance matrix of sample means and covariances
    new.call$NACOV <- Gamma      
  }
  if(estimator %in% c("WLS", "DWLS")) {
    # Weighted Least Squares, adjust the weight matrix.
    # Note that Gamma may be singular.
    if(estimator == "WLS") {
      # WLS uses the full Moore-Penrose inverse of Gamma.
      new.call$WLS.V <- lapply(Gamma, ginv)
    }
    else {
      # DWLS uses only the inverse diagonal of Gamma.
      new.call$WLS.V <- lapply(Gamma, get.dwls.weight)
      new.call$NACOV <- Gamma
    }
  }
  new.fit <- tryCatch(
    eval(as.call(new.call), envir=parent.frame()),
    error=function(e) e
  )
  if(inherits(new.fit, "error")) {
    first.error <- new.fit
    new.call$model <- lavaan::parTable(lavaan.fit)
    new.fit <- tryCatch(
      eval(as.call(new.call), envir=parent.frame()),
      error=function(e) e
    )
    if(inherits(new.fit, "error")) stop(first.error)
  }
  
  if(estimator %in% c("WLS", "DWLS")) return(new.fit) # We are done for WLS

  # For ML with robust se's, check that a possibly singular Gamma has not
  # created dependencies in the parameter estimates.
  # (Code below should really be implemented in lavaan...)
  evs.too.small <- sapply(Gamma, function(Gamma.g) {
    any(eigen(Gamma.g, only.values=TRUE)$values < (.Machine$double.eps*100))
  })
  if(any(evs.too.small)) {
    V.est <- lavaan::vcov(new.fit)
    if(any(Re(eigen(V.est, only.values=TRUE)$values) < (.Machine$double.eps*100))) {
      long.string  <- sprintf("Some of the standard errors may not be trustworthy.
        Some of the observed covariances or means are
        collinear, and this has generated collinearity in your
        parameter estimates.  This may be a sample size issue,
        missing data problem, or due to having too few
        clusters relative to the number of parameters. Problem
        encountered in group(s) %s",
        paste(which(evs.too.small), collapse=", "))
  
      warning(strwrap(long.string, width=9999, simplify=TRUE))#gotta love it
    }
  }

  new.fit
}

# Complex sampling analysis of ordinal SEM models.
# This first implementation is best validated for models where all observed
# model variables are ordinal. It also includes an experimental mixed
# ordinal/continuous path that delegates WLS statistic ordering to lavaan.
# It estimates the design-based covariance matrix from replicate weights
# generated by the survey package.
lavaan.survey.ordinal <-
  function(lavaan.fit, survey.design, ordered=NULL,
           estimator=c("WLSMV", "DWLS"),
           rep.type="auto", replicates=NULL,
           point.wls=c("auto", "lavaan", "design"),
           mi.pooling=c("auto", "parameters", "sample.statistics"),
           within.variance=c("auto", "replicate", "lavaan.robust", "naive"),
           standardized.se=c("lavaan", "replicate"),
           verbose=interactive()) {

  estimator <- match.arg(estimator)
  point.wls <- match.arg(point.wls)
  mi.pooling <- match.arg(mi.pooling)
  within.variance <- match.arg(within.variance)
  standardized.se <- match.arg(standardized.se)
  verbose <- isTRUE(verbose)

  ngroups <- lavInspect(lavaan.fit, "ngroups")
  meanstructure <- isTRUE(lavInspect(lavaan.fit, "options")$meanstructure)
  group.var <- NULL
  group.labels <- NULL
  if(ngroups > 1) {
    group.var <- lavInspect(lavaan.fit, "group")
    group.labels <- lavInspect(lavaan.fit, "group.label")
  }
  ov.names <- lavNames(lavaan.fit, type="ov", group=1)
  if(is.null(ordered)) ordered <- lavNames(lavaan.fit, type="ov.ord", group=1)
  ordered <- unique(ordered)
  if(length(ordered) == 0) {
    stop("No ordered variables were found. Please fit lavaan with ordered= or pass ordered=.")
  }
  if(!all(ordered %in% ov.names)) {
    stop("Ordered variables are not observed model variables: ",
         paste(setdiff(ordered, ov.names), collapse=", "))
  }
  all.ordinal <- all(ov.names %in% ordered)

  if(inherits(survey.design, "svyimputationList") &&
     all(vapply(survey.design$designs, inherits, logical(1), "svyrep.design"))) {
    rep.design <- survey.design
  }
  else if(inherits(survey.design, "svyrep.design")) {
    rep.design <- survey.design
  }
  else {
    rep.args <- list(design=survey.design, type=rep.type)
    if(!is.null(replicates)) rep.args$replicates <- replicates
    rep.design <- do.call(survey::as.svrepdesign, rep.args)
  }

  data <- if(inherits(rep.design, "svyimputationList")) {
    rep.design$designs[[1]]$variables
  }
  else {
    rep.design$variables
  }
  if(!all(ov.names %in% names(data))) {
    stop("The survey design is missing observed model variables: ",
         paste(setdiff(ov.names, names(data)), collapse=", "))
  }
  if(!is.null(group.var) && !group.var %in% names(data)) {
    stop("The survey design is missing the lavaan group variable: ", group.var)
  }
  is.mi.design <- inherits(rep.design, "svyimputationList")
  if(point.wls == "auto") {
    point.wls <- if(all.ordinal) "design" else "lavaan"
  }
  if(mi.pooling == "auto") {
    mi.pooling <- if(is.mi.design && !all.ordinal) "parameters" else "sample.statistics"
  }
  if(mi.pooling == "parameters") {
    within.variance <- resolve.parameter.mi.within.variance(within.variance,
                                                            rep.design,
                                                            point.wls)
    if(standardized.se == "replicate" && within.variance != "replicate") {
      stop("standardized.se = \"replicate\" requires within.variance = \"replicate\".",
           call.=FALSE)
    }
  }
  else {
    within.variance <- "none"
    if(standardized.se == "replicate") {
      warning("standardized.se = \"replicate\" is only used with mi.pooling = \"parameters\".",
              call.=FALSE)
    }
  }
  survey.info <- make.ordinal.survey.info(mode=if(all.ordinal) "ordinal" else "mixed ordinal/continuous",
                                          point.wls=point.wls,
                                          mi.pooling=if(is.mi.design) mi.pooling else "none",
                                          estimator=estimator,
                                          multiple.imputation=is.mi.design,
                                          within.variance=within.variance,
                                          standardized.se=standardized.se)

  if(mi.pooling == "parameters") {
    if(!is.mi.design) {
      stop("mi.pooling = \"parameters\" requires a svyimputationList survey design.")
    }
    fit <- pool.ordinal.mi.parameters(
      lavaan.fit=lavaan.fit, rep.design=rep.design, ordered=ordered,
      estimator=estimator, point.wls=point.wls, rep.type=rep.type,
      replicates=replicates, within.variance=within.variance,
      survey.info=survey.info, standardized.se=standardized.se,
      verbose=verbose
    )
    inform.ordinal.survey.info(survey.info)
    return(fit)
  }

  if(is.mi.design) {
    ordinal.stats <- lapply(rep.design$designs, get.ordinal.survey.stats,
                            ov.names=ov.names, ordered=ordered,
                            ngroups=ngroups, group.var=group.var,
                            group.labels=group.labels,
                            meanstructure=meanstructure,
                            estimator=estimator,
                            point.wls=point.wls)
    ordinal.stats <- pool.ordinal.mi.stats(ordinal.stats, ngroups=ngroups,
                                           group.labels=group.labels,
                                           point.wls=point.wls)
  }
  else {
    ordinal.stats <- get.ordinal.survey.stats(rep.design=rep.design,
                                              ov.names=ov.names,
                                              ordered=ordered,
                                              ngroups=ngroups,
                                              group.var=group.var,
                                              group.labels=group.labels,
                                              meanstructure=meanstructure,
                                              estimator=estimator,
                                              point.wls=point.wls)
  }

  point.stats <- ordinal.stats$point.stats
  Gamma <- ordinal.stats$Gamma
  WLS.V <- ordinal.stats$WLS.V

  new.call <- lavInspect(lavaan.fit, "call")
  new.call$data <- NULL
  new.call$model <- lavaan::parTable(lavaan.fit)
  new.call$sample.cov <- point.stats$sample.cov
  new.call$sample.mean <- point.stats$sample.mean
  new.call$sample.th <- point.stats$sample.th
  new.call$sample.nobs <- point.stats$nobs
  new.call$estimator <- estimator
  new.call$ordered <- ordered
  if(!is.null(group.var)) {
    new.call$group <- NULL
    new.call$group.label <- group.labels
  }
  new.call$sampling.weights <- NULL
  new.call$WLS.V <- WLS.V
  new.call$NACOV <- Gamma

  fit <- tryCatch(
    eval(as.call(new.call), envir=parent.frame()),
    error=function(e) {
      stop("Ordinal survey refit failed: ", conditionMessage(e), call.=FALSE)
    }
  )
  fit <- attach.ordinal.survey.info(fit, survey.info)
  inform.ordinal.survey.info(survey.info)
  fit
}

get.ordinal.survey.stats <- function(rep.design, ov.names, ordered,
                                     ngroups, group.var=NULL,
                                     group.labels=NULL,
                                     meanstructure=TRUE,
                                     estimator="WLSMV",
                                     point.wls="design") {
  data <- rep.design$variables
  point.weights <- stats::weights(rep.design, type="sampling")
  point.stats <- get.ordinal.stats(data=data, ov.names=ov.names,
                                   ordered=ordered, weights=point.weights,
                                   group=group.var, group.labels=group.labels,
                                   full=TRUE,
                                   meanstructure=meanstructure,
                                   estimator=estimator)

  rep.weights <- stats::weights(rep.design, type="analysis")
  if(ngroups == 1) {
    rep.stats <- matrix(NA_real_, nrow=ncol(rep.weights),
                        ncol=length(point.stats$wls.obs))
    colnames(rep.stats) <- names(point.stats$wls.obs)

    for(r in seq_len(ncol(rep.weights))) {
      rep.stats[r, ] <- get.ordinal.stats(data=data, ov.names=ov.names,
                                          ordered=ordered,
                                          weights=rep.weights[, r],
                                          wls.names=names(point.stats$wls.obs),
                                          meanstructure=meanstructure,
                                          estimator=estimator)$wls.obs
    }

    Gamma <- survey::svrVar(rep.stats, scale=rep.design$scale,
                            rscales=rep.design$rscales, mse=rep.design$mse,
                            coef=point.stats$wls.obs)
    Gamma <- as.matrix(Gamma) * point.stats$nobs
    dimnames(Gamma) <- list(names(point.stats$wls.obs), names(point.stats$wls.obs))
    WLS.V <- get.ordinal.point.wls(Gamma, point.stats, point.wls)

    return(list(point.stats=point.stats, Gamma=Gamma, WLS.V=WLS.V))
  }

  rep.stats <- vector("list", ngroups)
  names(rep.stats) <- group.labels
  for(g in seq_len(ngroups)) {
    rep.stats[[g]] <- matrix(NA_real_, nrow=ncol(rep.weights),
                             ncol=length(point.stats$wls.obs[[g]]))
    colnames(rep.stats[[g]]) <- names(point.stats$wls.obs[[g]])
  }

  for(r in seq_len(ncol(rep.weights))) {
    rep.wls <- get.ordinal.stats(data=data, ov.names=ov.names,
                                 ordered=ordered,
                                 weights=rep.weights[, r],
                                 group=group.var,
                                 group.labels=group.labels,
                                 wls.names=lapply(point.stats$wls.obs, names),
                                 meanstructure=meanstructure,
                                 estimator=estimator)$wls.obs
    for(g in seq_len(ngroups)) rep.stats[[g]][r, ] <- rep.wls[[g]]
  }

  Gamma <- vector("list", ngroups)
  WLS.V <- vector("list", ngroups)
  names(Gamma) <- names(WLS.V) <- group.labels
  for(g in seq_len(ngroups)) {
    Gamma[[g]] <- survey::svrVar(rep.stats[[g]], scale=rep.design$scale,
                                 rscales=rep.design$rscales, mse=rep.design$mse,
                                 coef=point.stats$wls.obs[[g]])
    Gamma[[g]] <- as.matrix(Gamma[[g]]) * point.stats$nobs[[g]]
    dimnames(Gamma[[g]]) <- list(names(point.stats$wls.obs[[g]]),
                                 names(point.stats$wls.obs[[g]]))
    WLS.V[[g]] <- get.ordinal.point.wls(Gamma[[g]], point.stats, point.wls, group=g)
  }

  list(point.stats=point.stats, Gamma=Gamma, WLS.V=WLS.V)
}

pool.ordinal.mi.stats <- function(ordinal.stats, ngroups, group.labels=NULL,
                                  point.wls="design") {
  if(length(ordinal.stats) == 0L) {
    stop("At least one imputed ordinal survey design is required.")
  }

  if(ngroups == 1) {
    point.stats <- pool.ordinal.mi.point.stats(
      lapply(ordinal.stats, `[[`, "point.stats")
    )
    Gamma <- pool.ordinal.mi.gamma(
      gamma.list=lapply(ordinal.stats, `[[`, "Gamma"),
      wls.list=lapply(ordinal.stats, function(x) x$point.stats$wls.obs),
      nobs.list=vapply(ordinal.stats, function(x) x$point.stats$nobs, numeric(1)),
      sample.nobs=point.stats$nobs
    )
    WLS.V <- get.ordinal.point.wls(Gamma, point.stats, point.wls)

    return(list(point.stats=point.stats, Gamma=Gamma, WLS.V=WLS.V))
  }

  point.stats <- pool.ordinal.mi.point.stats.by.group(
    lapply(ordinal.stats, `[[`, "point.stats"),
    group.labels=group.labels
  )
  Gamma <- WLS.V <- vector("list", ngroups)
  names(Gamma) <- names(WLS.V) <- group.labels
  for(g in seq_len(ngroups)) {
    Gamma[[g]] <- pool.ordinal.mi.gamma(
      gamma.list=lapply(ordinal.stats, function(x) x$Gamma[[g]]),
      wls.list=lapply(ordinal.stats, function(x) x$point.stats$wls.obs[[g]]),
      nobs.list=vapply(ordinal.stats, function(x) x$point.stats$nobs[[g]], numeric(1)),
      sample.nobs=point.stats$nobs[[g]]
    )
    WLS.V[[g]] <- get.ordinal.point.wls(Gamma[[g]], point.stats, point.wls, group=g)
  }

  list(point.stats=point.stats, Gamma=Gamma, WLS.V=WLS.V)
}

pool.ordinal.mi.point.stats <- function(point.stats.list) {
  first <- point.stats.list[[1]]
  sample.th <- average.named.vectors(lapply(point.stats.list, `[[`, "sample.th"))
  attr(sample.th, "th.idx") <- attr(first$sample.th, "th.idx")

  list(sample.cov=average.matrices(lapply(point.stats.list, `[[`, "sample.cov")),
       sample.mean=average.named.vectors(lapply(point.stats.list, `[[`, "sample.mean")),
       sample.th=sample.th,
       wls.obs=average.named.vectors(lapply(point.stats.list, `[[`, "wls.obs")),
       sample.wls.v=average.matrices(lapply(point.stats.list, `[[`, "sample.wls.v")),
       nobs=stats::median(vapply(point.stats.list, `[[`, numeric(1), "nobs")))
}

pool.ordinal.mi.point.stats.by.group <- function(point.stats.list, group.labels) {
  ngroups <- length(group.labels)
  first <- point.stats.list[[1]]

  sample.cov <- sample.mean <- sample.th <- wls.obs <- sample.wls.v <-
    vector("list", ngroups)
  names(sample.cov) <- names(sample.mean) <- names(sample.th) <-
    names(wls.obs) <- names(sample.wls.v) <- group.labels

  for(g in seq_len(ngroups)) {
    sample.cov[[g]] <- average.matrices(
      lapply(point.stats.list, function(x) x$sample.cov[[g]])
    )
    sample.mean[[g]] <- average.named.vectors(
      lapply(point.stats.list, function(x) x$sample.mean[[g]])
    )
    sample.th[[g]] <- average.named.vectors(
      lapply(point.stats.list, function(x) x$sample.th[[g]])
    )
    wls.obs[[g]] <- average.named.vectors(
      lapply(point.stats.list, function(x) x$wls.obs[[g]])
    )
    sample.wls.v[[g]] <- average.matrices(
      lapply(point.stats.list, function(x) x$sample.wls.v[[g]])
    )
  }
  attr(sample.th, "th.idx") <- attr(first$sample.th, "th.idx")

  nobs <- vapply(seq_len(ngroups), function(g) {
    stats::median(vapply(point.stats.list, function(x) x$nobs[[g]], numeric(1)))
  }, numeric(1))
  names(nobs) <- group.labels

  list(sample.cov=sample.cov,
       sample.mean=sample.mean,
       sample.th=sample.th,
       wls.obs=wls.obs,
       sample.wls.v=sample.wls.v,
       nobs=nobs)
}

pool.ordinal.mi.gamma <- function(gamma.list, wls.list, nobs.list, sample.nobs) {
  m <- length(gamma.list)
  ref.names <- names(wls.list[[1]])
  if(!all(vapply(wls.list, function(x) identical(names(x), ref.names), logical(1)))) {
    stop("Ordinal sample statistics do not match across imputations.")
  }

  gamma.list <- Map(function(Gamma, nobs) {
    Gamma / nobs * sample.nobs
  }, gamma.list, nobs.list)
  Gamma.within <- Reduce(`+`, gamma.list) / m

  if(m > 1L) {
    wls.df <- do.call(rbind, lapply(wls.list, `[`, ref.names))
    Gamma.between <- stats::cov(wls.df)
  }
  else {
    Gamma.between <- matrix(0, nrow=nrow(Gamma.within), ncol=ncol(Gamma.within))
  }
  dimnames(Gamma.between) <- dimnames(Gamma.within)

  Gamma.within + sample.nobs * ((m + 1) / m) * Gamma.between
}

get.ordinal.point.wls <- function(Gamma, point.stats, point.wls, group=NULL) {
  if(point.wls == "design") {
    return(get.dwls.weight(Gamma))
  }

  WLS.V <- if(is.null(group)) point.stats$sample.wls.v else point.stats$sample.wls.v[[group]]
  if(is.null(WLS.V)) {
    stop("lavaan WLS weights are not available for point.wls = \"lavaan\".")
  }
  dimnames(WLS.V) <- dimnames(Gamma)
  WLS.V
}

make.ordinal.survey.info <- function(mode, point.wls, mi.pooling, estimator,
                                     multiple.imputation,
                                     within.variance="none",
                                     standardized.se="lavaan") {
  list(function.name="lavaan.survey.ordinal",
       mode=mode,
       mi.pooling=mi.pooling,
       point.wls=point.wls,
       estimator=estimator,
       multiple.imputation=multiple.imputation,
       within.variance=within.variance,
       standardized.se=standardized.se)
}

attach.ordinal.survey.info <- function(fit, survey.info) {
  attr(fit, "lavaan.survey.info") <- survey.info
  fit
}

get.ordinal.survey.info <- function(fit) {
  survey.info <- attr(fit, "lavaan.survey.info")
  if(!is.null(survey.info)) return(survey.info)

  if(is.list(fit) && !is.null(fit$survey.info)) {
    return(fit$survey.info)
  }

  if(is.list(fit)) {
    fields <- c("mode", "mi.pooling", "point.wls", "estimator",
                "multiple.imputation", "within.variance", "standardized.se")
    survey.info <- list()
    for(field in fields) {
      value <- fit[[field]]
      if(!is.null(value)) survey.info[[field]] <- value
    }
    if(length(survey.info)) return(survey.info)
  }

  NULL
}

format.ordinal.survey.info <- function(survey.info) {
  out <- c(paste0("lavaan.survey.ordinal mode: ", survey.info$mode),
           paste0("MI pooling: ", survey.info$mi.pooling),
           paste0("Point WLS: ", survey.info$point.wls))
  if(!is.null(survey.info$within.variance) &&
     survey.info$mi.pooling == "parameters") {
    out <- c(out, paste0("Within variance: ", survey.info$within.variance))
  }
  if(!is.null(survey.info$standardized.se) &&
     survey.info$mi.pooling == "parameters") {
    out <- c(out, paste0("Standardized SE: ", survey.info$standardized.se))
  }
  out
}

inform.ordinal.survey.info <- function(survey.info) {
  message(paste(format.ordinal.survey.info(survey.info), collapse="\n"))
}

resolve.parameter.mi.within.variance <- function(within.variance, rep.design,
                                                 point.wls) {
  if(within.variance != "auto") return(within.variance)
  if(point.wls == "design") return("replicate")
  if(inherits(rep.design, "svyimputationList") &&
     all(vapply(rep.design$designs, inherits, logical(1), "svyrep.design"))) {
    return("replicate")
  }
  "naive"
}

pool.ordinal.mi.parameters <- function(lavaan.fit, rep.design, ordered,
                                       estimator, point.wls, rep.type,
                                       replicates,
                                       within.variance="naive",
                                       survey.info=NULL,
                                       standardized.se="lavaan",
                                       verbose=FALSE) {
  verbose <- isTRUE(verbose)
  if(point.wls != "lavaan" && within.variance != "replicate") {
    stop("within.variance = \"", within.variance, "\" is only available with ",
         "point.wls = \"lavaan\". The point.wls = \"design\" parameter-pooling ",
         "path uses the design-based replicate covariance from each ",
         "imputation.")
  }
  if(within.variance == "lavaan.robust") {
    stop("within.variance = \"lavaan.robust\" is not available for ordered ",
         "lavaan models in current lavaan versions because categorical + ",
         "clustered estimation is not supported. Use \"replicate\" or ",
         "\"naive\".")
  }
  fits <- if(point.wls == "lavaan") {
    if(verbose) {
      message("Fitting imputed datasets with lavaan sampling weights.")
    }
    lapply(seq_along(rep.design$designs), function(i) {
      if(verbose) {
        message("Fitting imputation ", i, "/", length(rep.design$designs), ".")
      }
      fit <- tryCatch(
        fit.ordinal.weighted.lavaan(
          design=rep.design$designs[[i]],
          lavaan.fit=lavaan.fit,
          ordered=ordered,
          estimator=estimator
        ),
        error=function(e) {
          stop("Imputation ", i, " failed: ", conditionMessage(e), call.=FALSE)
        }
      )
      if(!isTRUE(try(lavaan::lavInspect(fit, "converged"), silent=TRUE))) {
        stop("Imputation ", i, " did not converge.", call.=FALSE)
      }
      if(verbose) {
        message("Fitting imputation ", i, "/", length(rep.design$designs),
                " complete.")
      }
      fit
    })
  }
  else {
    if(verbose) {
      message("Fitting survey-adjusted imputed datasets.")
    }
    lapply(seq_along(rep.design$designs), function(i) {
      if(verbose) {
        message("Fitting imputation ", i, "/", length(rep.design$designs), ".")
      }
      fit <- tryCatch(
        suppressWarnings(lavaan.survey.ordinal(
          lavaan.fit=lavaan.fit,
          survey.design=rep.design$designs[[i]],
          ordered=ordered,
          estimator=estimator,
          rep.type=rep.type,
          replicates=replicates,
          point.wls=point.wls,
          mi.pooling="sample.statistics",
          verbose=FALSE
        )),
        error=function(e) {
          stop("Imputation ", i, " failed: ", conditionMessage(e), call.=FALSE)
        }
      )
      if(!isTRUE(try(lavaan::lavInspect(fit, "converged"), silent=TRUE))) {
        stop("Imputation ", i, " did not converge.", call.=FALSE)
      }
      if(verbose) {
        message("Fitting imputation ", i, "/", length(rep.design$designs),
                " complete.")
      }
      fit
    })
  }

  vcov.list <- NULL
  standardized.vcov.list <- NULL
  if(point.wls == "lavaan" && within.variance == "replicate") {
    if(verbose) {
      message("Estimating replicate within-imputation variance.")
    }
    standardized.types <- if(standardized.se == "replicate") {
      c("std.lv", "std.all", "std.nox")
    }
    else {
      NULL
    }
    vcov.list <- Map(estimate.replicate.parameter.vcov,
                     design=rep.design$designs,
                     point.fit=fits,
                     imputation.index=seq_along(rep.design$designs),
                     MoreArgs=list(lavaan.fit=lavaan.fit,
                                   ordered=ordered,
                                   estimator=estimator,
                                   imputations=length(rep.design$designs),
                                   standardized.types=standardized.types,
                                   verbose=verbose))
    if(standardized.se == "replicate") {
      standardized.vcov.list <- lapply(vcov.list, attr, "standardized.vcov")
    }
  }

  pooled <- pool.lavaan.mi.parameters(
    fits,
    vcov.list=vcov.list,
    standardized.vcov.list=standardized.vcov.list,
    df.complete=get.mi.complete.df(rep.design)
  )
  pooled$call <- lavInspect(lavaan.fit, "call")
  pooled$ordered <- ordered
  pooled <- attach.ordinal.survey.info(pooled, survey.info)
  pooled
}

fit.ordinal.weighted.lavaan <- function(design, lavaan.fit, ordered, estimator,
                                        weights=NULL) {
  data <- design$variables
  weight.name <- NULL
  if(is.null(weights)) {
    weight.name <- infer.sampling.weight.name(design, data)
    weights <- as.numeric(stats::weights(design, type="sampling"))
  }
  if(is.null(weight.name)) {
    weight.name <- ".lavaan.survey.sampling.weight"
    while(weight.name %in% names(data)) weight.name <- paste0(".", weight.name)
    data[[weight.name]] <- as.numeric(weights)
  }

  new.call <- lavInspect(lavaan.fit, "call")
  new.call$data <- data
  new.call$model <- lavaan::parTable(lavaan.fit)
  new.call$sample.cov <- NULL
  new.call$sample.mean <- NULL
  new.call$sample.th <- NULL
  new.call$sample.nobs <- NULL
  new.call$WLS.V <- NULL
  new.call$NACOV <- NULL
  new.call$estimator <- estimator
  new.call$ordered <- ordered
  new.call$sampling.weights <- weight.name

  eval(as.call(new.call), envir=parent.frame())
}

estimate.replicate.parameter.vcov <- function(design, point.fit, lavaan.fit,
                                              ordered, estimator,
                                              imputation.index=NULL,
                                              imputations=NULL,
                                              standardized.types=NULL,
                                              verbose=FALSE) {
  verbose <- isTRUE(verbose)
  ref.names <- names(lavaan::coef(point.fit))
  point.coef <- lavaan::coef(point.fit)[ref.names]
  standardized.types <- unique(standardized.types)
  if(length(standardized.types)) {
    point.std <- lapply(standardized.types, function(type) {
      x <- lavaan::standardizedSolution(point.fit, type=type, se=FALSE,
                                        output="data.frame")
      key <- lavaan.survey.mi.parameter.key(x)
      out <- x$est.std
      names(out) <- key
      out
    })
    names(point.std) <- standardized.types
  }
  else {
    point.std <- list()
  }
  rep.weights <- stats::weights(design, type="analysis")
  n.rep <- ncol(rep.weights)
  rep.coef <- matrix(NA_real_, nrow=n.rep, ncol=length(ref.names),
                     dimnames=list(NULL, ref.names))
  rep.std <- lapply(point.std, function(x) {
    matrix(NA_real_, nrow=n.rep, ncol=length(x), dimnames=list(NULL, names(x)))
  })

  report.every <- max(1L, floor(n.rep / 20L))
  imputation.label <- if(!is.null(imputation.index) && !is.null(imputations)) {
    paste0(" for imputation ", imputation.index, "/", imputations)
  }
  else {
    ""
  }
  if(verbose) {
    message("Replicate fits", imputation.label, ": 0/", n.rep)
  }

  for(r in seq_len(n.rep)) {
    fit.r <- try(suppressWarnings(fit.ordinal.weighted.lavaan(
      design=design,
      lavaan.fit=lavaan.fit,
      ordered=ordered,
      estimator=estimator,
      weights=rep.weights[, r]
    )), silent=TRUE)
    if(!inherits(fit.r, "try-error") &&
       isTRUE(try(lavaan::lavInspect(fit.r, "converged"), silent=TRUE))) {
      rep.coef[r, ] <- lavaan::coef(fit.r)[ref.names]
      for(type in names(rep.std)) {
        std.r <- try(lavaan::standardizedSolution(fit.r, type=type, se=FALSE,
                                                  output="data.frame"),
                     silent=TRUE)
        if(!inherits(std.r, "try-error")) {
          std.r.key <- lavaan.survey.mi.parameter.key(std.r)
          rep.std[[type]][r, ] <- std.r$est.std[match(colnames(rep.std[[type]]),
                                                      std.r.key)]
        }
      }
    }

    if(verbose && (r == n.rep || r %% report.every == 0L)) {
      message("Replicate fits", imputation.label, ": ", r, "/", n.rep)
    }
  }

  keep <- stats::complete.cases(rep.coef)
  n.total <- nrow(rep.coef)
  n.keep <- sum(keep)
  min.keep <- min(10L, max(2L, ceiling(n.total / 2)))
  if(n.keep < min.keep) {
    stop("Only ", n.keep, " of ", n.total,
         " replicate lavaan fits converged; cannot estimate a stable ",
         "within-imputation replicate variance.", call.=FALSE)
  }
  failure.rate <- 1 - n.keep / n.total
  if(failure.rate > 0.10) {
    warning(sprintf(
      "%.0f%% of replicate lavaan fits did not converge. The replicate-based within-imputation variance estimate may underestimate uncertainty.",
      100 * failure.rate
    ), call.=FALSE)
  }
  if(verbose) {
    message("Replicate fits", imputation.label, " complete: ", n.keep, "/",
            n.total, " converged.")
  }

  V <- survey::svrVar(rep.coef[keep, , drop=FALSE],
                      scale=design$scale,
                      rscales=design$rscales[keep],
                      mse=design$mse,
                      coef=point.coef)
  V <- as.matrix(V)
  dimnames(V) <- list(ref.names, ref.names)
  if(length(rep.std)) {
    standardized.vcov <- lapply(names(rep.std), function(type) {
      keep.std <- keep & stats::complete.cases(rep.std[[type]])
      n.keep.std <- sum(keep.std)
      if(n.keep.std < min.keep) {
        stop("Only ", n.keep.std, " of ", n.total,
             " replicate standardized lavaan fits converged for ", type,
             "; cannot estimate a stable standardized replicate variance.",
             call.=FALSE)
      }
      failure.rate.std <- 1 - n.keep.std / n.total
      if(failure.rate.std > 0.10) {
        warning(sprintf(
          "%.0f%% of replicate standardized lavaan fits did not converge for %s. The replicate-based standardized within-imputation variance estimate may underestimate uncertainty.",
          100 * failure.rate.std, type
        ), call.=FALSE)
      }
      V.std <- survey::svrVar(rep.std[[type]][keep.std, , drop=FALSE],
                              scale=design$scale,
                              rscales=design$rscales[keep.std],
                              mse=design$mse,
                              coef=point.std[[type]])
      V.std <- as.matrix(V.std)
      dimnames(V.std) <- list(names(point.std[[type]]), names(point.std[[type]]))
      V.std
    })
    names(standardized.vcov) <- names(rep.std)
    attr(V, "standardized.vcov") <- standardized.vcov
  }
  V
}

infer.sampling.weight.name <- function(design, data) {
  sampling.weights <- as.numeric(stats::weights(design, type="sampling"))
  if(length(sampling.weights) != nrow(data)) return(NULL)

  numeric.names <- names(data)[vapply(data, is.numeric, logical(1))]
  for(nm in numeric.names) {
    candidate <- as.numeric(data[[nm]])
    if(length(candidate) == length(sampling.weights) &&
       isTRUE(all.equal(candidate, sampling.weights,
                        tolerance=sqrt(.Machine$double.eps),
                        check.attributes=FALSE))) {
      return(nm)
    }
  }

  NULL
}

get.mi.complete.df <- function(rep.design) {
  designs <- if(inherits(rep.design, "svyimputationList")) {
    rep.design$designs
  }
  else {
    list(rep.design)
  }
  dfs <- vapply(designs, function(design) {
    out <- try(survey::degf(design), silent=TRUE)
    if(inherits(out, "try-error") || length(out) != 1L || is.na(out)) {
      Inf
    }
    else {
      as.numeric(out)
    }
  }, numeric(1))
  finite.dfs <- dfs[is.finite(dfs)]
  if(length(finite.dfs)) min(finite.dfs) else Inf
}

pool.lavaan.mi.parameters <- function(fits, vcov.list=NULL,
                                      standardized.vcov.list=NULL,
                                      df.complete=Inf) {
  if(length(fits) == 0L) {
    stop("At least one fitted lavaan object is required for MI parameter pooling.")
  }
  if(!is.null(vcov.list) && length(vcov.list) != length(fits)) {
    stop("vcov.list (", length(vcov.list), ") and fits (", length(fits),
         ") must have the same length.")
  }
  if(!is.null(standardized.vcov.list) &&
     length(standardized.vcov.list) != length(fits)) {
    stop("standardized.vcov.list (", length(standardized.vcov.list),
         ") and fits (", length(fits), ") must have the same length.")
  }

  ref.names <- names(lavaan::coef(fits[[1]]))
  coef.list <- lapply(fits, function(fit) lavaan::coef(fit)[ref.names])
  if(!all(vapply(coef.list, function(x) identical(names(x), ref.names), logical(1)))) {
    stop("Free parameter names do not match across imputations.")
  }
  coef.matrix <- do.call(rbind, coef.list)
  coef.pooled <- colMeans(coef.matrix)
  names(coef.pooled) <- ref.names

  if(is.null(vcov.list)) {
    vcov.list <- lapply(fits, function(fit) {
      V <- as.matrix(lavaan::vcov(fit))
      V[ref.names, ref.names, drop=FALSE]
    })
  }
  else {
    vcov.list <- lapply(vcov.list, function(V) {
      V <- as.matrix(V)
      V[ref.names, ref.names, drop=FALSE]
    })
  }
  vcov.within <- Reduce(`+`, vcov.list) / length(vcov.list)
  if(length(fits) > 1L) {
    vcov.between <- stats::cov(coef.matrix)
  }
  else {
    vcov.between <- matrix(0, nrow=length(ref.names), ncol=length(ref.names),
                           dimnames=list(ref.names, ref.names))
  }
  dimnames(vcov.between) <- dimnames(vcov.within)
  vcov.total <- vcov.within + ((length(fits) + 1) / length(fits)) * vcov.between

  fit.measures <- pool.lavaan.mi.fit.measures(fits)
  out <- list(coef=coef.pooled,
              vcov=vcov.total,
              vcov.within=vcov.within,
              vcov.between=vcov.between,
              fits=fits,
              fit.measures=fit.measures,
              m=length(fits),
              df.complete=df.complete)
  if(!is.null(standardized.vcov.list)) {
    out$standardized.vcov.within <-
      pool.lavaan.mi.standardized.within.vcov(standardized.vcov.list)
  }
  class(out) <- "lavaan.survey.mi"
  out
}

pool.lavaan.mi.standardized.within.vcov <- function(standardized.vcov.list) {
  types <- Reduce(intersect, lapply(standardized.vcov.list, names))
  out <- lapply(types, function(type) {
    mats <- lapply(standardized.vcov.list, `[[`, type)
    ref.dimnames <- dimnames(mats[[1]])
    if(!all(vapply(mats, function(m) identical(dimnames(m), ref.dimnames),
                   logical(1)))) {
      stop("Standardized replicate covariance names do not match across imputations for ",
           type, ".", call.=FALSE)
    }
    Reduce(`+`, mats) / length(mats)
  })
  names(out) <- types
  out
}

pool.lavaan.mi.fit.measures <- function(fits) {
  measure.names <- c("chisq", "chisq.scaled", "df", "df.scaled",
                     "cfi", "cfi.scaled", "rmsea", "rmsea.scaled", "srmr")
  fit.list <- lapply(fits, function(fit) {
    out <- try(lavaan::fitMeasures(fit, measure.names), silent=TRUE)
    if(inherits(out, "try-error")) {
      stats::setNames(rep(NA_real_, length(measure.names)), measure.names)
    }
    else {
      out
    }
  })
  fit.matrix <- do.call(rbind, fit.list)
  colMeans(fit.matrix, na.rm=TRUE)
}

.onLoad <- function(libname, pkgname) {
  register.semPlotModel.method <- function(...) {
    if(isNamespaceLoaded("semPlot")) {
      registerS3method(
        "semPlotModel",
        "lavaan.survey.mi",
        semPlotModel.lavaan.survey.mi,
        envir=asNamespace("semPlot")
      )
    }
  }

  register.semPlotModel.method()
  setHook(packageEvent("semPlot", "onLoad"), register.semPlotModel.method,
          action="append")
}

coef.lavaan.survey.mi <- function(object, ...) {
  object$coef
}

vcov.lavaan.survey.mi <- function(object, ...) {
  object$vcov
}

fitMeasures.lavaan.survey.mi <- function(object, fit.measures="all",
                                         baseline.model=NULL,
                                         h1.model=NULL, fm.args=list(),
                                         output=c("vector", "data.frame",
                                                  "list"),
                                         ...) {
  output <- match.arg(output)
  available <- object$fit.measures %||% numeric(0)
  if(length(fit.measures) == 1L && identical(tolower(fit.measures), "all")) {
    out <- available
  }
  else {
    out <- available[fit.measures]
    missing <- setdiff(fit.measures, names(available))
    if(length(missing)) {
      out[missing] <- NA_real_
      out <- out[fit.measures]
    }
  }

  if(output == "list") {
    return(as.list(out))
  }
  if(output == "data.frame") {
    return(data.frame(measure=names(out), value=as.numeric(out),
                      row.names=NULL, check.names=FALSE))
  }
  out
}

semPlotModel.lavaan.survey.mi <- function(object, ...) {
  if(!requireNamespace("semPlot", quietly=TRUE)) {
    stop("The semPlot package is required for semPlotModel() or semPaths() ",
         "with lavaan.survey.mi objects.", call.=FALSE)
  }

  slots <- lavaan.survey.mi.semPlotModel.slots(object)
  new("semPlotModel",
      Pars=slots$Pars,
      Vars=slots$Vars,
      Thresholds=slots$Thresholds,
      Computed=slots$Computed,
      ObsCovs=slots$ObsCovs,
      ImpCovs=slots$ImpCovs,
      Original=slots$Original)
}

lavaan.survey.mi.semPlotModel.slots <- function(object) {
  if(!inherits(object, "lavaan.survey.mi")) {
    stop("A lavaan.survey.mi object is required.", call.=FALSE)
  }
  if(length(object$fits) == 0L) {
    stop("The lavaan.survey.mi object does not contain fitted lavaan models.",
         call.=FALSE)
  }

  ref.fit <- object$fits[[1]]
  pt <- as.data.frame(lavaan::parTable(ref.fit), stringsAsFactors=FALSE)
  pars <- parameterEstimates.lavaan.survey.mi(
    object,
    se=FALSE,
    zstat=FALSE,
    pvalue=FALSE,
    ci=FALSE,
    standardized=TRUE,
    remove.system.eq=FALSE,
    remove.eq=FALSE,
    remove.ineq=FALSE,
    remove.def=FALSE
  )

  if(!"label" %in% names(pt)) pt$label <- rep("", nrow(pt))
  if(!"group" %in% names(pt)) pt$group <- ""
  if(!"free" %in% names(pt)) pt$free <- 0L
  if(!"est" %in% names(pt)) pt$est <- NA_real_

  idx <- lavaan.survey.mi.match.parameter.rows(pt, pars)

  est <- pt$est
  std <- rep(NA_real_, nrow(pt))
  has.match <- !is.na(idx)
  est[has.match] <- pars$est[idx[has.match]]
  if("std.all" %in% names(pars)) {
    std[has.match] <- pars$std.all[idx[has.match]]
  }

  sem.pars <- data.frame(
    label=pt$label,
    lhs=ifelse(pt$op %in% c("~", "~1"), pt$rhs, pt$lhs),
    edge="--",
    rhs=ifelse(pt$op %in% c("~", "~1"), pt$lhs, pt$rhs),
    est=est,
    std=std,
    group=pt$group,
    fixed=pt$free == 0L,
    par=pt$free,
    stringsAsFactors=FALSE
  )

  sem.pars$edge[pt$op %in% c("~~", "~*~")] <- "<->"
  sem.pars$edge[pt$op == "~"] <- "~>"
  sem.pars$edge[pt$op == "=~"] <- "->"
  sem.pars$edge[pt$op == "~1"] <- "int"
  sem.pars$edge[grepl("\\|", pt$op)] <- "|"

  thresholds <- sem.pars[grepl("\\|", sem.pars$edge), -(3:4), drop=FALSE]
  keep <- !pt$op %in% c(":=", "<", ">", "==", "|")
  sem.pars <- sem.pars[keep, , drop=FALSE]

  var.names <- lavaan::lavNames(ref.fit, type="ov")
  fact.names <- lavaan::lavNames(ref.fit, type="lv")
  fact.names <- fact.names[!fact.names %in% var.names]
  vars <- data.frame(
    name=c(var.names, fact.names),
    manifest=c(var.names, fact.names) %in% var.names,
    exogenous=NA,
    stringsAsFactors=FALSE
  )

  list(
    Pars=sem.pars,
    Vars=vars,
    Thresholds=thresholds,
    Computed=TRUE,
    ObsCovs=lavaan.survey.mi.semPlotModel.covs(object$fits, "sampstat"),
    ImpCovs=lavaan.survey.mi.semPlotModel.covs(object$fits, "implied"),
    Original=list(object)
  )
}

lavaan.survey.mi.semPlotModel.covs <- function(fits, what) {
  covs <- lapply(fits, lavaan.survey.mi.extract.semPlotModel.covs, what=what)
  ref.length <- length(covs[[1]])
  if(!all(vapply(covs, length, integer(1)) == ref.length)) {
    return(covs[[1]])
  }

  out <- lapply(seq_len(ref.length), function(g) {
    mats <- lapply(covs, `[[`, g)
    if(!all(vapply(mats, is.matrix, logical(1)))) return(mats[[1]])
    average.matrices(mats)
  })
  group.labels <- lavaan::lavInspect(fits[[1]], "group.label")
  if(length(group.labels) == length(out)) names(out) <- group.labels
  out
}

lavaan.survey.mi.extract.semPlotModel.covs <- function(fit, what) {
  element <- if(isTRUE(lavaan::lavInspect(fit, "options")$conditional.x)) {
    "res.cov"
  }
  else {
    "cov"
  }
  stats <- lavaan::lavTech(fit, what, add.labels=TRUE)
  if(is.matrix(stats)) return(list(stats))

  lapply(stats, function(x) {
    if(is.matrix(x)) return(x)
    x[[element]] %||% x$cov
  })
}

parameterEstimates <- function(object, ...) {
  if(inherits(object, "lavaan.survey.mi")) {
    return(parameterEstimates.lavaan.survey.mi(object, ...))
  }
  lavaan::parameterEstimates(object, ...)
}

standardizedSolution <- function(object, ...) {
  if(inherits(object, "lavaan.survey.mi")) {
    return(standardizedSolution.lavaan.survey.mi(object, ...))
  }
  lavaan::standardizedSolution(object, ...)
}

parameterEstimates.lavaan.survey.mi <- function(object, se=TRUE, zstat=TRUE,
                                                pvalue=TRUE, ci=TRUE,
                                                standardized=FALSE,
                                                standardized.se=c("lavaan",
                                                                  "replicate"),
                                                level=0.95, plabel=FALSE,
                                                cov.std=TRUE,
                                                remove.system.eq=TRUE,
                                                remove.eq=TRUE,
                                                remove.ineq=TRUE,
                                                remove.def=FALSE,
                                                remove.nonfree=FALSE,
                                                output=c("data.frame", "table"),
                                                header=FALSE, ...) {
  output <- match.arg(output)
  standardized.se <- match.arg(standardized.se)
  pt <- as.data.frame(lavaan::parTable(object$fits[[1]]),
                      stringsAsFactors=FALSE)
  ref.names <- names(object$coef)
  free.order <- unique(pt$free[pt$free > 0L])
  if(length(free.order) != length(ref.names)) {
    stop("Could not align pooled free parameters with the lavaan parameter table.",
         call.=FALSE)
  }
  free.map <- seq_along(free.order)
  names(free.map) <- as.character(free.order)
  free.idx <- ifelse(pt$free > 0L, free.map[as.character(pt$free)], NA_integer_)

  est <- pt$est
  has.free <- !is.na(free.idx)
  est[has.free] <- object$coef[free.idx[has.free]]

  out <- data.frame(lhs=pt$lhs,
                    op=pt$op,
                    rhs=pt$rhs,
                    stringsAsFactors=FALSE,
                    check.names=FALSE)
  if("block" %in% names(pt)) out$block <- pt$block
  if("group" %in% names(pt)) out$group <- pt$group
  if("label" %in% names(pt) && any(nchar(pt$label) > 0L)) out$label <- pt$label
  if(isTRUE(plabel) && "plabel" %in% names(pt)) out$plabel <- pt$plabel
  out$est <- est

  df <- barnard.rubin.df(object)
  if(isTRUE(se)) {
    se.out <- rep(0, nrow(pt))
    se.out[has.free] <- sqrt(diag(object$vcov))[free.idx[has.free]]
    out$se <- se.out

    if(isTRUE(zstat)) {
      t.out <- rep(NA_real_, nrow(pt))
      t.out[has.free] <- out$est[has.free] / out$se[has.free]
      out$t <- t.out
      df.out <- rep(NA_real_, nrow(pt))
      df.out[has.free] <- df[free.idx[has.free]]
      out$df <- df.out

      if(isTRUE(pvalue)) {
        p.out <- rep(NA_real_, nrow(pt))
        p.out[has.free] <- 2 * stats::pt(abs(t.out[has.free]),
                                          df=df.out[has.free],
                                          lower.tail=FALSE)
        normal.ref <- has.free & is.infinite(df.out)
        p.out[normal.ref] <- 2 * stats::pnorm(abs(t.out[normal.ref]),
                                              lower.tail=FALSE)
        out$pvalue <- p.out
      }
    }

    if(isTRUE(ci)) {
      ci.lower <- ci.upper <- rep(NA_real_, nrow(pt))
      alpha <- 1 - level
      df.out <- rep(NA_real_, nrow(pt))
      df.out[has.free] <- df[free.idx[has.free]]
      crit <- stats::qt(1 - alpha / 2, df=df.out)
      normal.ref <- is.infinite(df.out)
      crit[normal.ref] <- stats::qnorm(1 - alpha / 2)
      ci.lower[has.free] <- out$est[has.free] - crit[has.free] * out$se[has.free]
      ci.upper[has.free] <- out$est[has.free] + crit[has.free] * out$se[has.free]
      fixed <- !has.free
      ci.lower[fixed] <- out$est[fixed]
      ci.upper[fixed] <- out$est[fixed]
      out$ci.lower <- ci.lower
      out$ci.upper <- ci.upper
    }
  }

  keep <- rep(TRUE, nrow(out))
  if(isTRUE(remove.system.eq)) keep <- keep & !(out$op == "==" & (pt$user %||% 0L) != 1L)
  if(isTRUE(remove.eq)) keep <- keep & out$op != "=="
  if(isTRUE(remove.ineq)) keep <- keep & !out$op %in% c("<", ">")
  if(isTRUE(remove.def)) keep <- keep & out$op != ":="
  if(isTRUE(remove.nonfree)) {
    keep <- keep & !(pt$free == 0L & !out$op %in% c("==", ">", "<", ":="))
  }
  out <- out[keep, , drop=FALSE]

  if(isTRUE(standardized)) {
    for(type in c("std.lv", "std.all", "std.nox")) {
      std <- lavaan.survey.mi.standardized.table(
        object=object,
        type=type,
        se=FALSE,
        zstat=FALSE,
        pvalue=FALSE,
        ci=FALSE,
        standardized.se=standardized.se,
        level=level,
        cov.std=cov.std,
        remove.eq=remove.eq,
        remove.ineq=remove.ineq,
        remove.def=remove.def
      )
      out[[type]] <- std$est.std[lavaan.survey.mi.match.parameter.rows(out, std)]
    }
  }

  if(identical(output, "table")) {
    class(out) <- c("lavaan.data.frame", "data.frame")
  }
  out
}

standardizedSolution.lavaan.survey.mi <- function(object, type="std.all",
                                                  se=TRUE, zstat=TRUE,
                                                  pvalue=TRUE, ci=TRUE,
                                                  standardized.se=c("lavaan",
                                                                    "replicate"),
                                                  level=0.95, cov.std=TRUE,
                                                  remove.eq=TRUE,
                                                  remove.ineq=TRUE,
                                                  remove.def=FALSE,
                                                  partable=NULL, GLIST=NULL,
                                                  est=NULL,
                                                  output="data.frame") {
  standardized.se <- match.arg(standardized.se)
  if(!is.null(partable) || !is.null(GLIST) || !is.null(est)) {
    warning("`partable`, `GLIST`, and `est` are ignored for lavaan.survey.mi objects.",
            call.=FALSE)
  }
  out <- lavaan.survey.mi.standardized.table(
    object=object,
    type=type,
    se=se,
    zstat=zstat,
    pvalue=pvalue,
    ci=ci,
    standardized.se=standardized.se,
    level=level,
    cov.std=cov.std,
    remove.eq=remove.eq,
    remove.ineq=remove.ineq,
    remove.def=remove.def
  )
  attr(out, "standardized.type") <- type
  attr(out, "standardized.se") <- standardized.se
  attr(out, "se.method") <- if(standardized.se == "replicate") {
    "Rubin pooling with replicate-based standardized within-imputation variance"
  }
  else {
    "Rubin pooling of per-imputation lavaan standardized-solution SEs"
  }
  if(identical(output, "table")) {
    class(out) <- c("lavaan.data.frame", "data.frame")
  }
  out
}

lavaan.survey.mi.standardized.table <- function(object, type="std.all",
                                                se=TRUE, zstat=TRUE,
                                                pvalue=TRUE, ci=TRUE,
                                                standardized.se="lavaan",
                                                level=0.95, cov.std=TRUE,
                                                remove.eq=TRUE,
                                                remove.ineq=TRUE,
                                                remove.def=FALSE) {
  need.se <- isTRUE(se) || isTRUE(zstat) || isTRUE(pvalue) || isTRUE(ci)
  standardized.se <- match.arg(standardized.se, c("lavaan", "replicate"))
  std.list <- lapply(object$fits, function(fit) {
    lavaan::standardizedSolution(
      fit,
      type=type,
      se=need.se,
      zstat=FALSE,
      pvalue=FALSE,
      ci=FALSE,
      level=level,
      cov.std=cov.std,
      remove.eq=remove.eq,
      remove.ineq=remove.ineq,
      remove.def=remove.def,
      output="data.frame"
    )
  })

  ref <- std.list[[1]]
  ref.key <- lavaan.survey.mi.parameter.key(ref)
  std.list <- lapply(std.list, function(x) {
    x.key <- lavaan.survey.mi.parameter.key(x)
    idx <- match(ref.key, x.key)
    if(anyNA(idx)) {
      stop("Could not align standardized solution rows across imputations.",
           call.=FALSE)
    }
    x[idx, , drop=FALSE]
  })

  est.matrix <- do.call(rbind, lapply(std.list, function(x) x$est.std))
  est.pooled <- colMeans(est.matrix, na.rm=TRUE)
  out <- ref
  out$est.std <- as.numeric(est.pooled)

  if(need.se) {
    if(standardized.se == "replicate") {
      within.vcov <- object$standardized.vcov.within[[type]] %||% NULL
      if(is.null(within.vcov)) {
        stop("Replicate standardized SEs are not available for this object. ",
             "Refit with standardized.se = \"replicate\" and ",
             "within.variance = \"replicate\".", call.=FALSE)
      }
      within.vcov <- within.vcov[ref.key, ref.key, drop=FALSE]
    }
    else {
      if(!all(vapply(std.list, function(x) "se" %in% names(x), logical(1)))) {
        stop("Per-imputation lavaan standardized SEs are not available.",
             call.=FALSE)
      }
      se.matrix <- do.call(rbind, lapply(std.list, function(x) x$se))
      within.var <- colMeans(se.matrix^2, na.rm=TRUE)
      within.var[!is.finite(within.var)] <- 0
      within.vcov <- diag(within.var, nrow=length(within.var))
      dimnames(within.vcov) <- list(ref.key, ref.key)
    }
    if(length(std.list) > 1L) {
      between.vcov <- stats::cov(est.matrix, use="pairwise.complete.obs")
    }
    else {
      between.vcov <- matrix(0, nrow=length(est.pooled), ncol=length(est.pooled),
                             dimnames=list(ref.key, ref.key))
    }
    between.vcov[!is.finite(between.vcov)] <- 0
    dimnames(between.vcov) <- list(ref.key, ref.key)
    total.vcov <- within.vcov +
      ((length(std.list) + 1) / length(std.list)) * between.vcov
    total.var <- diag(total.vcov)
    total.var[total.var < 0] <- NA_real_
    std.se <- sqrt(total.var)
    df <- barnard.rubin.df(list(
      coef=stats::setNames(est.pooled, ref.key),
      vcov=total.vcov,
      vcov.within=within.vcov,
      vcov.between=between.vcov,
      m=object$m,
      df.complete=object$df.complete
    ))

    if(isTRUE(se)) {
      out$se <- std.se
    }
    else if("se" %in% names(out)) {
      out$se <- NULL
    }
    if(isTRUE(zstat)) {
      out$t <- out$est.std / std.se
      out$df <- df
    }
    if(isTRUE(pvalue)) {
      t.value <- out$est.std / std.se
      p <- 2 * stats::pt(abs(t.value), df=df, lower.tail=FALSE)
      normal.ref <- is.infinite(df)
      p[normal.ref] <- 2 * stats::pnorm(abs(t.value[normal.ref]),
                                        lower.tail=FALSE)
      out$pvalue <- p
    }
    if(isTRUE(ci)) {
      alpha <- 1 - level
      crit <- stats::qt(1 - alpha / 2, df=df)
      normal.ref <- is.infinite(df)
      crit[normal.ref] <- stats::qnorm(1 - alpha / 2)
      out$ci.lower <- out$est.std - crit * std.se
      out$ci.upper <- out$est.std + crit * std.se
    }
  }
  else {
    drop.cols <- intersect(c("se", "z", "pvalue", "ci.lower", "ci.upper"),
                           names(out))
    out[drop.cols] <- NULL
  }
  if("z" %in% names(out)) {
    out$z <- NULL
  }
  rownames(out) <- NULL
  out
}

lavaan.survey.mi.parameter.key <- function(x) {
  key.cols <- intersect(c("lhs", "op", "rhs", "group"), names(x))
  do.call(paste, c(x[key.cols], sep="\r"))
}

lavaan.survey.mi.match.parameter.rows <- function(x, table) {
  key.cols <- c("lhs", "op", "rhs")
  if("group" %in% names(x) && "group" %in% names(table)) {
    key.cols <- c(key.cols, "group")
  }
  key.cols <- intersect(key.cols, intersect(names(x), names(table)))
  match(lavaan.survey.mi.parameter.key.with.cols(x, key.cols),
        lavaan.survey.mi.parameter.key.with.cols(table, key.cols))
}

lavaan.survey.mi.parameter.key.with.cols <- function(x, key.cols) {
  do.call(paste, c(x[key.cols], sep="\r"))
}

print.lavaan.survey.mi <- function(x, ...) {
  survey.info <- get.ordinal.survey.info(x)
  mode <- if(!is.null(survey.info)) survey.info$mode %||% "survey SEM" else "survey SEM"
  cat("lavaan.survey multiply imputed ", mode, " fit\n", sep="")
  if(!is.null(survey.info)) {
    cat(paste0("  ", format.ordinal.survey.info(survey.info)), sep="\n")
    cat("\n")
  }
  cat("  Imputations:", x$m, "\n")
  cat("  Parameter pooling: Rubin\n")
  cat("  Point WLS:", survey.info$point.wls %||% "unknown", "\n")
  cat("  Within variance:", survey.info$within.variance %||% "unknown", "\n")
  cat("  Standardized SE:", survey.info$standardized.se %||% "lavaan", "\n")
  invisible(x)
}

summary.lavaan.survey.mi <- function(object, fit.measures=TRUE,
                                     standardized=FALSE, ci=FALSE,
                                     level=0.95, ...) {
  out <- lavaan.survey.mi.parameter.table(object, ci=ci, level=level,
                                          standardized=standardized)
  attr(out, "fit.measures") <- if(isTRUE(fit.measures)) object$fit.measures else numeric(0)
  attr(out, "survey.info") <- get.ordinal.survey.info(object)
  attr(out, "m") <- object$m
  attr(out, "df.complete") <- object$df.complete
  attr(out, "ci") <- isTRUE(ci)
  attr(out, "level") <- level
  attr(out, "standardized.requested") <- isTRUE(standardized)
  class(out) <- c("summary.lavaan.survey.mi", "data.frame")
  print(out)
  invisible(out)
}

lavaan.survey.mi.parameter.table <- function(object, ci=FALSE, level=0.95,
                                             standardized=FALSE) {
  se <- sqrt(diag(object$vcov))
  t.value <- object$coef / se
  df <- barnard.rubin.df(object)
  p <- 2 * stats::pt(abs(t.value), df=df, lower.tail=FALSE)
  normal.ref <- is.infinite(df)
  p[normal.ref] <- 2 * stats::pnorm(abs(t.value[normal.ref]),
                                    lower.tail=FALSE)
  parsed <- parse.lavaan.parameter.names(names(object$coef))
  out <- data.frame(lhs=parsed$lhs,
                    op=parsed$op,
                    rhs=parsed$rhs,
                    Estimate=object$coef,
                    Std.Err=se,
                    t.value=t.value,
                    df=df,
                    P.value=p,
                    check.names=FALSE)
  if(isTRUE(ci)) {
    alpha <- 1 - level
    crit <- stats::qt(1 - alpha / 2, df=df)
    normal.ref <- is.infinite(df)
    crit[normal.ref] <- stats::qnorm(1 - alpha / 2)
    out$ci.lower <- out$Estimate - crit * out$Std.Err
    out$ci.upper <- out$Estimate + crit * out$Std.Err
  }
  if(isTRUE(standardized)) {
    pe.std <- parameterEstimates.lavaan.survey.mi(
      object,
      se=FALSE,
      zstat=FALSE,
      pvalue=FALSE,
      ci=FALSE,
      standardized=TRUE,
      remove.nonfree=TRUE
    )
    out$Std.all <- pe.std$std.all[
      lavaan.survey.mi.match.parameter.rows(out, pe.std)
    ]
  }
  rownames(out) <- names(object$coef)
  out
}

parse.lavaan.parameter.names <- function(x) {
  out <- data.frame(lhs=x, op="", rhs="", stringsAsFactors=FALSE)
  for(i in seq_along(x)) {
    nm <- x[[i]]
    if(grepl("~1$", nm)) {
      out$lhs[[i]] <- sub("~1$", "", nm)
      out$op[[i]] <- "~1"
      out$rhs[[i]] <- ""
    }
    else if(grepl("=~", nm, fixed=TRUE)) {
      parts <- strsplit(nm, "=~", fixed=TRUE)[[1]]
      out$lhs[[i]] <- parts[[1]]
      out$op[[i]] <- "=~"
      out$rhs[[i]] <- parts[[2]]
    }
    else if(grepl("~~", nm, fixed=TRUE)) {
      parts <- strsplit(nm, "~~", fixed=TRUE)[[1]]
      out$lhs[[i]] <- parts[[1]]
      out$op[[i]] <- "~~"
      out$rhs[[i]] <- parts[[2]]
    }
    else if(grepl("|", nm, fixed=TRUE)) {
      parts <- strsplit(nm, "|", fixed=TRUE)[[1]]
      out$lhs[[i]] <- parts[[1]]
      out$op[[i]] <- "|"
      out$rhs[[i]] <- parts[[2]]
    }
    else if(grepl("~", nm, fixed=TRUE)) {
      parts <- strsplit(nm, "~", fixed=TRUE)[[1]]
      out$lhs[[i]] <- parts[[1]]
      out$op[[i]] <- "~"
      out$rhs[[i]] <- parts[[2]]
    }
  }
  out
}

lavaan.survey.mi.parameter.blocks <- function(x) {
  ifelse(x$op == "=~", "Latent Variables",
         ifelse(x$op == "~", "Regressions",
                ifelse(x$op == "~~" & x$lhs == x$rhs, "Variances",
                       ifelse(x$op == "~~", "Covariances",
                              ifelse(x$op == "|", "Thresholds",
                                     ifelse(x$op == "~1", "Intercepts", "Other"))))))
}

format.mi.number <- function(x, digits=3) {
  out <- formatC(x, format="f", digits=digits)
  out[is.na(x)] <- ""
  out
}

format.mi.df <- function(x, digits=2) {
  out <- formatC(x, format="f", digits=digits)
  out[is.infinite(x)] <- "Inf"
  out[is.na(x)] <- ""
  out
}

format.mi.p <- function(x) {
  out <- ifelse(is.na(x), "",
                ifelse(x < .001, "<.001",
                       sub("^0", "", formatC(x, format="f", digits=3))))
  out
}

format.mi.table.cell <- function(x, width, align="left") {
  x <- as.character(x)
  x[is.na(x)] <- ""
  pad <- pmax(0L, width - nchar(x, type="width", allowNA=FALSE))
  if(align == "right") {
    return(paste0(strrep(" ", pad), x))
  }
  if(align == "center") {
    left <- floor(pad / 2)
    right <- pad - left
    return(paste0(strrep(" ", left), x, strrep(" ", right)))
  }
  paste0(x, strrep(" ", pad))
}

print.mi.table <- function(x, align=NULL, spacing=2L) {
  x <- as.data.frame(x, stringsAsFactors=FALSE, check.names=FALSE)
  if(nrow(x) == 0L || ncol(x) == 0L) return(invisible(x))
  x[] <- lapply(x, as.character)
  if(is.null(align)) {
    align <- rep("left", ncol(x))
    names(align) <- names(x)
  }
  align <- align[names(x)]
  align[is.na(align)] <- "left"
  widths <- vapply(seq_along(x), function(i) {
    max(nchar(c(names(x)[[i]], x[[i]]), type="width", allowNA=FALSE))
  }, integer(1))

  rows <- rbind(names(x), as.matrix(x))
  formatted <- do.call(cbind, lapply(seq_along(x), function(i) {
    format.mi.table.cell(rows[, i], widths[[i]], align[[i]])
  }))
  line.sep <- strrep(" ", spacing)
  cat(apply(formatted, 1L, paste, collapse=line.sep), sep="\n")
  cat("\n")
  invisible(x)
}

print.summary.lavaan.survey.mi <- function(x, digits=3, ...) {
  survey.info <- attr(x, "survey.info")
  fit.measures <- attr(x, "fit.measures")
  m <- attr(x, "m")
  df.complete <- attr(x, "df.complete")
  ci <- isTRUE(attr(x, "ci"))
  level <- attr(x, "level") %||% 0.95
  standardized.requested <- isTRUE(attr(x, "standardized.requested"))

  cat("\nlavaan.survey MI Summary\n")
  cat(strrep("=", 72), "\n", sep="")
  if(!is.null(survey.info)) {
    cat("Mode:            ", survey.info$mode %||% "unknown", "\n", sep="")
    cat("Estimator:       ", survey.info$estimator %||% "unknown", "\n", sep="")
    cat("MI pooling:      ", survey.info$mi.pooling %||% "unknown", "\n", sep="")
    cat("Point WLS:       ", survey.info$point.wls %||% "unknown", "\n", sep="")
    cat("Within variance: ", survey.info$within.variance %||% "unknown", "\n", sep="")
    cat("Standardized SE: ", survey.info$standardized.se %||% "lavaan", "\n", sep="")
  }
  cat("Imputations:     ", m %||% NA_integer_, "\n", sep="")
  if(!is.null(df.complete)) {
    cat("Survey df:       ", format.mi.df(df.complete), "\n", sep="")
  }
  cat("Tests:           t statistics with Barnard-Rubin degrees of freedom\n")
  if(standardized.requested) {
    cat("Standardized:    pooled std.all column shown\n")
  }

  if(length(fit.measures)) {
    cat("\nPooled Fit Measures (averaged over imputations)\n")
    fit.order <- c("chisq.scaled", "df.scaled", "cfi.scaled",
                   "rmsea.scaled", "srmr", "chisq", "df", "cfi", "rmsea")
    fit.names <- intersect(fit.order, names(fit.measures))
    fit.out <- data.frame(
      Measure=fit.names,
      Value=format.mi.number(fit.measures[fit.names], digits=3),
      row.names=NULL,
      check.names=FALSE
    )
    print.mi.table(fit.out, align=c(Measure="left", Value="right"))
  }

  blocks <- lavaan.survey.mi.parameter.blocks(x)
  block.order <- c("Latent Variables", "Regressions", "Covariances",
                   "Intercepts", "Thresholds", "Variances", "Other")
  for(block in block.order[block.order %in% blocks]) {
    cat("\n", block, "\n", sep="")
    cat(strrep("-", nchar(block)), "\n", sep="")
    rows <- x[blocks == block, , drop=FALSE]
    out <- data.frame(
      lhs=rows$lhs,
      op=rows$op,
      rhs=rows$rhs,
      Estimate=format.mi.number(rows$Estimate, digits=digits),
      Std.Err=format.mi.number(rows$Std.Err, digits=digits),
      t=format.mi.number(rows$t.value, digits=digits),
      df=format.mi.df(rows$df),
      `P(>|t|)`=format.mi.p(rows$P.value),
      check.names=FALSE
    )
    if("Std.all" %in% names(rows)) {
      out$Std.all <- format.mi.number(rows$Std.all, digits=digits)
    }
    if(ci) {
      out[[paste0("ci.lower.", round(100 * level))]] <-
        format.mi.number(rows$ci.lower, digits=digits)
      out[[paste0("ci.upper.", round(100 * level))]] <-
        format.mi.number(rows$ci.upper, digits=digits)
    }
    align <- rep("right", ncol(out))
    names(align) <- names(out)
    align[c("lhs", "rhs")] <- "left"
    align["op"] <- "center"
    print.mi.table(out, align=align)
  }
  invisible(x)
}

barnard.rubin.df <- function(object) {
  npar <- length(object$coef)
  ref.names <- names(object$coef)
  if(is.null(object$vcov.within) || is.null(object$vcov.between) ||
     is.null(object$vcov) || is.null(object$m) || object$m <= 1L) {
    out <- rep(object$df.complete %||% Inf, npar)
    names(out) <- ref.names
    return(out)
  }

  m <- object$m
  total.var <- diag(object$vcov)
  between.var <- diag(object$vcov.between)
  lambda <- ((m + 1) / m) * between.var / total.var
  lambda[!is.finite(lambda)] <- NA_real_
  lambda <- pmin(pmax(lambda, 0), 1)

  df.old <- rep(Inf, length(lambda))
  has.missing.info <- !is.na(lambda) & lambda > .Machine$double.eps
  df.old[has.missing.info] <- (m - 1) / (lambda[has.missing.info]^2)

  df.complete <- object$df.complete %||% Inf
  if(is.finite(df.complete)) {
    df.obs <- ((df.complete + 1) / (df.complete + 3)) *
      df.complete * (1 - lambda)
    out <- combine.df(df.old, df.obs)
  }
  else {
    out <- df.old
  }

  out[is.na(lambda) | out <= .Machine$double.eps] <- NA_real_
  names(out) <- ref.names
  out
}

combine.df <- function(df1, df2) {
  out <- (df1 * df2) / (df1 + df2)
  old.inf <- is.infinite(df1) & is.finite(df2)
  obs.inf <- is.finite(df1) & is.infinite(df2)
  both.inf <- is.infinite(df1) & is.infinite(df2)
  out[old.inf] <- df2[old.inf]
  out[obs.inf] <- df1[obs.inf]
  out[both.inf] <- Inf
  out
}

`%||%` <- function(x, y) {
  if(is.null(x)) y else x
}

average.named.vectors <- function(x) {
  ref.names <- names(x[[1]])
  if(!all(vapply(x, function(v) identical(names(v), ref.names), logical(1)))) {
    stop("Ordinal sample statistic names do not match across imputations.")
  }
  out <- Reduce(`+`, lapply(x, `[`, ref.names)) / length(x)
  names(out) <- ref.names
  out
}

average.matrices <- function(x) {
  ref.dimnames <- dimnames(x[[1]])
  if(!all(vapply(x, function(m) identical(dimnames(m), ref.dimnames), logical(1)))) {
    stop("Ordinal covariance names do not match across imputations.")
  }
  out <- Reduce(`+`, x) / length(x)
  dimnames(out) <- ref.dimnames
  out
}

get.ordinal.stats <- function(data, ov.names, ordered, weights,
                              group=NULL, group.labels=NULL,
                              wls.names=NULL, full=FALSE,
                              meanstructure=TRUE,
                              estimator="WLSMV") {
  weight.name <- ".lavaan.survey.weight"
  while(weight.name %in% names(data)) weight.name <- paste0(".", weight.name)

  keep.names <- ov.names
  if(!is.null(group)) keep.names <- c(keep.names, group)
  data <- data[, keep.names, drop=FALSE]
  for(v in ordered) {
    if(!is.ordered(data[[v]])) data[[v]] <- base::ordered(data[[v]])
  }
  if(!is.null(group)) data[[group]] <- factor(data[[group]], levels=group.labels)
  data[[weight.name]] <- as.numeric(weights)

  keep <- stats::complete.cases(data[, c(keep.names, weight.name), drop=FALSE])
  data <- data[keep, , drop=FALSE]
  if(nrow(data) == 0) stop("No complete observations are available for the ordinal model variables.")
  all.ordinal <- all(ov.names %in% ordered)

  if(full) {
    unrestricted <- lavaan::lavCor(data, ordered=ordered, sampling.weights=weight.name,
                                   group=group,
                                   output="lavaan", missing="listwise",
                                   meanstructure=meanstructure,
                                   estimator=estimator)
    sample.stats <- lavInspect(unrestricted, "sampstat")
    wls.obs <- lavInspect(unrestricted, "wls.obs")
    sample.wls.v <- lavInspect(unrestricted, "wls.v")

    if(is.null(group)) {
      th.idx <- unrestricted@SampleStats@th.idx[[1]]
      names(th.idx) <- names(sample.stats$th)
      sample.th <- sample.stats$th
      attr(sample.th, "th.idx") <- th.idx

      return(list(sample.cov=sample.stats$cov,
                  sample.mean=sample.stats$mean,
                  sample.th=sample.th,
                  wls.obs=wls.obs,
                  sample.wls.v=sample.wls.v,
                  nobs=lavInspect(unrestricted, "nobs")))
    }

    sample.th <- lapply(sample.stats, `[[`, "th")
    th.idx <- unrestricted@SampleStats@th.idx
    for(g in seq_along(th.idx)) names(th.idx[[g]]) <- names(sample.th[[g]])
    attr(sample.th, "th.idx") <- th.idx

    return(list(sample.cov=lapply(sample.stats, `[[`, "cov"),
                sample.mean=lapply(sample.stats, `[[`, "mean"),
                sample.th=sample.th,
                wls.obs=wls.obs,
                sample.wls.v=sample.wls.v,
                nobs=lavInspect(unrestricted, "nobs")))
  }

  if(!all.ordinal) {
    unrestricted <- lavaan::lavCor(data, ordered=ordered, sampling.weights=weight.name,
                                   group=group,
                                   output="lavaan", missing="listwise",
                                   meanstructure=meanstructure,
                                   estimator=estimator)
    wls.obs <- lavInspect(unrestricted, "wls.obs")

    if(is.null(group)) {
      if(!is.null(wls.names)) wls.obs <- wls.obs[wls.names]
      if(anyNA(wls.obs)) stop("Could not match ordinal sample statistics across replicate weights.")

      return(list(wls.obs=wls.obs))
    }

    if(!is.null(wls.names)) {
      for(g in seq_along(wls.obs)) wls.obs[[g]] <- wls.obs[[g]][wls.names[[g]]]
    }
    if(any(unlist(lapply(wls.obs, anyNA)))) {
      stop("Could not match ordinal sample statistics across replicate weights.")
    }

    return(list(wls.obs=wls.obs))
  }

  sample.stats <- lavaan::lavCor(data, ordered=ordered, sampling.weights=weight.name,
                                 group=group,
                                 output="sampstat", missing="listwise",
                                 meanstructure=meanstructure,
                                 estimator=estimator)
  if(is.null(group)) {
    wls.obs <- get.ordinal.wls.obs(sample.stats)
    if(!is.null(wls.names)) wls.obs <- wls.obs[wls.names]
    if(anyNA(wls.obs)) stop("Could not match ordinal sample statistics across replicate weights.")

    return(list(wls.obs=wls.obs))
  }

  wls.obs <- lapply(sample.stats, get.ordinal.wls.obs)
  if(!is.null(wls.names)) {
    for(g in seq_along(wls.obs)) wls.obs[[g]] <- wls.obs[[g]][wls.names[[g]]]
  }
  if(any(unlist(lapply(wls.obs, anyNA)))) {
    stop("Could not match ordinal sample statistics across replicate weights.")
  }

  list(wls.obs=wls.obs)
}

get.ordinal.wls.obs <- function(sample.stats) {
  COV <- sample.stats$cov
  cov.idx <- which(upper.tri(COV), arr.ind=TRUE)
  cov.names <- paste(rownames(COV)[cov.idx[, 1]], colnames(COV)[cov.idx[, 2]], sep="~~")
  cov.obs <- COV[upper.tri(COV)]
  names(cov.obs) <- cov.names
  c(sample.stats$th, cov.obs)
}

# Safely subset survey designs by lavaan group labels. Avoid parse(text=...),
# because group labels can contain quotes or other non-syntactic characters.
survey_subset_by_group <- function(survey.design, group.var, group.label) {
  eval(substitute(
    subset(survey.design, get(GROUPVAR) == GROUPLABEL),
    list(GROUPVAR=group.var, GROUPLABEL=group.label)
  ))
}

# Obtain residuals from a lavaan fit object, concatenating means w/ covariances
# (used in Yuan-Bentler correction)
get.residuals <- function(fit, group=1, meanstructure=TRUE) {
    r <- lavaan::residuals(fit)
    
    if(!is.null(r$cov)) {
      r.g <- r
    }
    else {
      r.g <- r[[group]]
      if(is.null(r.g)) stop("Could not find residuals for group ", group, ".")
    }

    if(meanstructure) c(r.g$mean, lavaan::lav_matrix_vech(r.g$cov))
    else lavaan::lav_matrix_vech(r.g$cov)
}

# Obtain the diagonal DWLS weight matrix from Gamma.
get.dwls.weight <- function(Gamma) {
    gamma.diag <- diag(Gamma)
    w <- rep(0, length(gamma.diag))
    use <- gamma.diag > .Machine$double.eps
    w[use] <- 1 / gamma.diag[use]
    W <- diag(w, nrow=length(w))
    dimnames(W) <- dimnames(Gamma)
    W
}

# Obtain sample size from multiply imputed svydesign object.
# In case sample size differs over imputations, takes median over imputations.
# For multiple-group models, lavaan.survey subsets the svyimputationList by
# group before calling this helper.
get.sample.nobs  <- function(svy.imp.design) {
  if(!inherits(svy.imp.design, "svyimputationList")) {
    stop("Expected a svyimputationList object.")
  }
  nobs.imp <- vapply(svy.imp.design$designs, function(des) {
    nrow(des$variables)
  }, numeric(1))
  median(nobs.imp)
}

# Use the pFsum function from the survey package to obtain p value for the 
#   overall model fit using an F reference distribution where the 
#   denominator degrees of freedom is the design degrees of freedom.  
# An anonymous reviewer for J Stat Software suggested that 
#  "in surveys with small numbers of primary sampling units this sort of 
#   correction has often improved the 
#   behavior of tests in other contexts."
# The eigenvalues of the U.Gamma matrix will be the coefficients in the 
#   mixture of F's distribution (Skinner, Holt & Smith, pp. 86-87).
pval.pFsum <- function(lavaan.fit, survey.design, method = "saddlepoint") {
  # Check that a robust lavaan test with U.Gamma information is present.
  test.options <- lavInspect(lavaan.fit, "options")$test
  robust.tests <- c("satorra.bentler", "scaled.shifted", "mean.var.adjusted")
  if(!any(test.options %in% robust.tests)) {
    stop("Please refit the model with a robust lavaan test (MLM, MLMV, or MLMVS).")
  }
  
  UGamma <- get.ugamma.matrix(lavaan.fit)
  real.eigen.values <- Re(eigen(UGamma, only.values = TRUE)$values)

  return(survey::pFsum(x=fitMeasures(lavaan.fit, "chisq"), df=rep(1, length(real.eigen.values)), 
                a=real.eigen.values, ddf=survey::degf(survey.design), lower.tail=FALSE,
                method=method))
}

get.ugamma.matrix <- function(lavaan.fit) {
  as.ugamma.matrix(lavTech(lavaan.fit, "ugamma"))
}

as.ugamma.matrix <- function(UGamma) {
  if(is.list(UGamma)) {
    if(length(UGamma) == 0L ||
       !all(vapply(UGamma, is.square.matrix, logical(1)))) {
      stop("Could not extract U.Gamma as a square matrix or list of square matrices.")
    }
    UGamma <- lavaan::lav_matrix_bdiag(UGamma)
  }

  if(!is.square.matrix(UGamma)) {
    stop("Could not extract U.Gamma as a square matrix.")
  }

  UGamma
}

is.square.matrix <- function(x) {
  is.matrix(x) && length(dim(x)) == 2L && nrow(x) == ncol(x)
}
