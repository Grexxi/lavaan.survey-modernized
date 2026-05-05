# Complex sampling analysis of SEM models
# Daniel Oberski, 2015-11-03

lavaan.survey <- 
  function(lavaan.fit, survey.design, 
           estimator=c("MLM", "MLMV", "MLMVS", "WLS", "DWLS", "ML"),
           estimator.gamma=c("default","Yuan-Bentler")) {
  
  # Not all estimators in lavaan make sense to use here, therefore matching args
  estimator <- match.arg(estimator) 
  if(estimator=="ML") warning("Estimator 'ML' will not correct standard errors and chi-square statistic.")
  estimator.gamma <- match.arg(estimator.gamma) # Smoothing Gamma or not
  
  # Names of the observed variables (same for each group)
  ov.names <- lavaanNames(lavaan.fit, type="ov", group=1)
  
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
        subset(survey.design, eval(parse(text=sprintf("%s == '%s'", 
                                                      lavInspect(lavaan.fit, "group"), 
                                                      lavInspect(lavaan.fit, "group.label")[[g]]))))
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
      
      # Same for mean vector
      sample.mean.g <- svymean(ov.formula, design=survey.design.g, na.rm=TRUE)  
      Gamma.mean.g <- attr(sample.mean.g, "var")
      
      # Join asy. variance matrices for means and covariances
      # TODO add offdiag
      Gamma.g <- lavaan::lav_matrix_bdiag(Gamma.mean.g, Gamma.cov.g)
      
      Gamma.g <- Gamma.g * sample.nobs.g # lavaan wants nobs * Gamma.
      
      # Since the above nonparametric estimate of Gamma can be unstable, Yuan
      # and Bentler suggested a model-smoothed estimate of it, optional here:
      if(estimator.gamma == "Yuan-Bentler") {
        r <- get.residuals(lavaan.fit, group=g) # Iff these asy = 0, all will be well...
        Gamma.g <- Gamma.g + (sample.nobs.g/(sample.nobs.g - 1)) * (r %*% t(r))
      }
      # Get rid of attribute, preventing errors caused by lazy evaluation
      # (This has to be at the end or lazy evaluation mayhem will ensue)
      attr(sample.cov.g, "var") <- NULL
      tmp  <- as.vector(sample.mean.g)
      names(tmp) <- names(sample.mean.g)
      sample.mean.g <- tmp	
      
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
      sample.nobs.g <- get.sample.nobs(survey.design.g, lavInspect(lavaan.fit, "group"))
      
      # Retrieve point and variance estimates per imputation.
      stats.list <- lapply(survey.design.g$designs, get.stats.design,
                           sample.nobs.g=sample.nobs.g)
      m  <- length(stats.list) # no. imputation
      
      # Point estimates are average over imputations
      sample.cov.list <- lapply(stats.list, `[[`, 'sample.cov.g')
      sample.cov.g <- Reduce(`+`, sample.cov.list) / m
      cov.df <- Reduce(`rbind`, lapply(sample.cov.list, lavaan::lav_matrix_vech))
      sample.mean.list <- lapply(stats.list, `[[`, 'sample.mean.g')
      sample.mean.g <- Reduce(`+`, sample.mean.list) / m
      mean.df <- Reduce(`rbind`, sample.mean.list)
      
      # Variance estimates depend on within- and between-imputation variance:
      Gamma.within  <- Reduce(`+`, lapply(stats.list, `[[`, 'Gamma.g')) / m
      if(m > 1L) {
        Gamma.between <- cov(cbind(mean.df, cov.df))
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
    sample.mean[[g]] <- stats$sample.mean.g
    sample.nobs[[g]] <- sample.nobs.g
  } # End of loop over groups

  new.call <- lavInspect(lavaan.fit, "call")
  new.call$data <- NULL                # Remove any data argument
  new.call$sample.cov <- sample.cov    # Set survey covariances
  new.call$sample.mean <- sample.mean  # Set survey means
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
    if(grepl("object .* not found", conditionMessage(new.fit))) {
      new.call$model <- lavaan::parTable(lavaan.fit)
      new.fit <- eval(as.call(new.call))
    }
    else {
      stop(new.fit)
    }
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
# This first implementation supports models where all observed model variables
# are ordinal. It estimates thresholds and polychoric correlations using lavaan,
# and estimates their design-based covariance matrix from replicate weights
# generated by the survey package.
lavaan.survey.ordinal <-
  function(lavaan.fit, survey.design, ordered=NULL,
           estimator=c("WLSMV", "DWLS"),
           rep.type="auto", replicates=NULL) {

  estimator <- match.arg(estimator)

  ngroups <- lavInspect(lavaan.fit, "ngroups")
  group.var <- NULL
  group.labels <- NULL
  if(ngroups > 1) {
    group.var <- lavInspect(lavaan.fit, "group")
    group.labels <- lavInspect(lavaan.fit, "group.label")
  }
  ov.names <- lavaan::lavNames(lavaan.fit, type="ov", group=1)
  if(is.null(ordered)) ordered <- lavaan::lavNames(lavaan.fit, type="ov.ord", group=1)
  if(length(ordered) == 0) {
    stop("No ordered variables were found. Please fit lavaan with ordered= or pass ordered=.")
  }
  if(!all(ov.names %in% ordered)) {
    stop("lavaan.survey.ordinal currently requires all observed model variables to be ordered.")
  }

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

  if(inherits(rep.design, "svyimputationList")) {
    ordinal.stats <- lapply(rep.design$designs, get.ordinal.survey.stats,
                            ov.names=ov.names, ordered=ordered,
                            ngroups=ngroups, group.var=group.var,
                            group.labels=group.labels)
    ordinal.stats <- pool.ordinal.mi.stats(ordinal.stats, ngroups=ngroups,
                                           group.labels=group.labels)
  }
  else {
    ordinal.stats <- get.ordinal.survey.stats(rep.design=rep.design,
                                              ov.names=ov.names,
                                              ordered=ordered,
                                              ngroups=ngroups,
                                              group.var=group.var,
                                              group.labels=group.labels)
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

  eval(as.call(new.call))
}

get.ordinal.survey.stats <- function(rep.design, ov.names, ordered,
                                     ngroups, group.var=NULL,
                                     group.labels=NULL) {
  data <- rep.design$variables
  point.weights <- stats::weights(rep.design, type="sampling")
  point.stats <- get.ordinal.stats(data=data, ov.names=ov.names,
                                   ordered=ordered, weights=point.weights,
                                   group=group.var, group.labels=group.labels,
                                   full=TRUE)

  rep.weights <- stats::weights(rep.design, type="analysis")
  if(ngroups == 1) {
    rep.stats <- matrix(NA_real_, nrow=ncol(rep.weights),
                        ncol=length(point.stats$wls.obs))
    colnames(rep.stats) <- names(point.stats$wls.obs)

    for(r in seq_len(ncol(rep.weights))) {
      rep.stats[r, ] <- get.ordinal.stats(data=data, ov.names=ov.names,
                                          ordered=ordered,
                                          weights=rep.weights[, r],
                                          wls.names=names(point.stats$wls.obs))$wls.obs
    }

    Gamma <- survey::svrVar(rep.stats, scale=rep.design$scale,
                            rscales=rep.design$rscales, mse=rep.design$mse,
                            coef=point.stats$wls.obs)
    Gamma <- as.matrix(Gamma) * point.stats$nobs
    dimnames(Gamma) <- list(names(point.stats$wls.obs), names(point.stats$wls.obs))
    WLS.V <- get.dwls.weight(Gamma)

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
                                 wls.names=lapply(point.stats$wls.obs, names))$wls.obs
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
    WLS.V[[g]] <- get.dwls.weight(Gamma[[g]])
  }

  list(point.stats=point.stats, Gamma=Gamma, WLS.V=WLS.V)
}

pool.ordinal.mi.stats <- function(ordinal.stats, ngroups, group.labels=NULL) {
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
    WLS.V <- get.dwls.weight(Gamma)

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
    WLS.V[[g]] <- get.dwls.weight(Gamma[[g]])
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
       nobs=stats::median(vapply(point.stats.list, `[[`, numeric(1), "nobs")))
}

pool.ordinal.mi.point.stats.by.group <- function(point.stats.list, group.labels) {
  ngroups <- length(group.labels)
  first <- point.stats.list[[1]]

  sample.cov <- sample.mean <- sample.th <- wls.obs <- vector("list", ngroups)
  names(sample.cov) <- names(sample.mean) <- names(sample.th) <-
    names(wls.obs) <- group.labels

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
                              wls.names=NULL, full=FALSE) {
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

  if(full) {
    unrestricted <- lavaan::lavCor(data, ordered=ordered, sampling.weights=weight.name,
                                   group=group,
                                   output="lavaan", missing="listwise")
    sample.stats <- lavInspect(unrestricted, "sampstat")
    wls.obs <- lavInspect(unrestricted, "wls.obs")

    if(is.null(group)) {
      th.idx <- unrestricted@SampleStats@th.idx[[1]]
      names(th.idx) <- names(sample.stats$th)
      sample.th <- sample.stats$th
      attr(sample.th, "th.idx") <- th.idx

      return(list(sample.cov=sample.stats$cov,
                  sample.mean=sample.stats$mean,
                  sample.th=sample.th,
                  wls.obs=wls.obs,
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
                nobs=lavInspect(unrestricted, "nobs")))
  }

  sample.stats <- lavaan::lavCor(data, ordered=ordered, sampling.weights=weight.name,
                                 group=group,
                                 output="sampstat", missing="listwise")
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

# Obtain residuals from a lavaan fit object, concatenating means w/ covariances
# (used in Yuan-Bentler correction)
get.residuals <- function(fit, group=1) {
    r <- lavaan::residuals(fit)
    
    if(!is.null(r$cov)) {
      r.g <- r
    }
    else {
      r.g <- r[[group]]
      if(is.null(r.g)) stop("Could not find residuals for group ", group, ".")
    }
    
    c(r.g$mean, lavaan::lav_matrix_vech(r.g$cov))
}

# Obtain the diagonal DWLS weight matrix from Gamma.
get.dwls.weight <- function(Gamma) {
    gamma.diag <- diag(Gamma)
    w <- rep(0, length(gamma.diag))
    use <- abs(gamma.diag) > .Machine$double.eps
    w[use] <- 1 / gamma.diag[use]
    W <- diag(w, nrow=length(w))
    dimnames(W) <- dimnames(Gamma)
    W
}

# Obtain sample size from multiply imputed svydesign object.
# In case sample size differs over imputations, takes median over imputations.
# TODO: Does not work with multiple group yet.
get.sample.nobs  <- function(svy.imp.design, group=NULL) {
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
  # Check that Satorra-Bentler or Satterthwaite adjustment is present
  test.options <- lavInspect(lavaan.fit, "options")$test
  if(!any(test.options %in% c("satorra.bentler", "mean.var.adjusted", "Satterthwaite"))) {
    stop("Please refit the model with Satorra-Bentler (MLM) or Satterthwaite (MLMVS) adjustment.") 
  }
  
  UGamma <- lavTech(lavaan.fit, "ugamma")
  real.eigen.values <- Re(eigen(UGamma, only.values = TRUE)$values)

  return(survey::pFsum(x=fitMeasures(lavaan.fit, "chisq"), df=rep(1, length(real.eigen.values)), 
                a=real.eigen.values, ddf=survey::degf(survey.design), lower.tail=FALSE,
                method=method))
}
