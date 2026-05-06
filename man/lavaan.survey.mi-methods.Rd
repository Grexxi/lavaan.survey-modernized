\name{lavaan.survey.mi-methods}
\alias{lavaan.survey.mi-methods}
\alias{parameterEstimates}
\alias{parameterEstimates.lavaan.survey.mi}
\alias{standardizedSolution}
\alias{standardizedSolution.lavaan.survey.mi}
\alias{fitMeasures.lavaan.survey.mi}
\title{Extractor methods for lavaan.survey multiple-imputation fits}
\description{
  Extract pooled parameter estimates, standardized solutions, and
  pooled fit-measure averages from \code{lavaan.survey.mi} objects returned by
  \code{\link{lavaan.survey.ordinal}} when Rubin parameter pooling is used.
}
\usage{
parameterEstimates(object, ...)

standardizedSolution(object, ...)

\method{parameterEstimates}{lavaan.survey.mi}(
  object, se = TRUE, zstat = TRUE, pvalue = TRUE, ci = TRUE,
  standardized = FALSE, level = 0.95, plabel = FALSE,
  cov.std = TRUE, remove.system.eq = TRUE, remove.eq = TRUE,
  remove.ineq = TRUE, remove.def = FALSE, remove.nonfree = FALSE,
  output = c("data.frame", "table"), header = FALSE, ...
)

\method{standardizedSolution}{lavaan.survey.mi}(
  object, type = "std.all", se = TRUE, zstat = TRUE, pvalue = TRUE,
  ci = TRUE, level = 0.95, cov.std = TRUE, remove.eq = TRUE,
  remove.ineq = TRUE, remove.def = FALSE, partable = NULL, GLIST = NULL,
  est = NULL, output = "data.frame"
)

\method{fitMeasures}{lavaan.survey.mi}(
  object, fit.measures = "all", baseline.model = NULL, h1.model = NULL,
  fm.args = list(),
  output = c("vector", "data.frame", "list"), ...
)
}
\arguments{
  \item{object}{
    A \code{lavaan.survey.mi} object.
  }
  \item{se, zstat, pvalue, ci, level, plabel, cov.std, remove.system.eq,
        remove.eq, remove.ineq, remove.def, remove.nonfree, output, header}{
    Arguments following the spirit of \code{\link[lavaan]{parameterEstimates}}.
    For \code{lavaan.survey.mi} objects, test statistics are \code{t}
    statistics with Barnard-Rubin degrees of freedom.
  }
  \item{standardized}{
    Logical. If \code{TRUE}, \code{parameterEstimates()} adds pooled
    \code{std.lv}, \code{std.all}, and \code{std.nox} columns.
  }
  \item{type}{
    Standardization type for \code{standardizedSolution()}, such as
    \code{"std.all"}, \code{"std.lv"}, or \code{"std.nox"}.
  }
  \item{partable, GLIST, est}{
    Included for compatibility with \code{\link[lavaan]{standardizedSolution}}.
    These arguments are ignored for pooled \code{lavaan.survey.mi} objects.
  }
  \item{fit.measures}{
    Character vector of fit measures to extract, or \code{"all"}.
  }
  \item{baseline.model, h1.model, fm.args}{
    Included for compatibility with the \code{\link[lavaan]{fitMeasures}}
    generic. They are currently ignored for pooled \code{lavaan.survey.mi}
    objects, whose fit measures are stored as pooled averages.
  }
  \item{...}{
    Additional arguments. For non-\code{lavaan.survey.mi} objects, the
    exported \code{parameterEstimates()} wrapper delegates to
    \code{\link[lavaan]{parameterEstimates}}.
  }
}
\value{
  \code{parameterEstimates()} returns a data frame with lavaan-style columns
  \code{lhs}, \code{op}, \code{rhs}, \code{est}, and, when requested,
  \code{se}, \code{t}, \code{df}, \code{pvalue}, \code{ci.lower}, and
  \code{ci.upper}. With \code{standardized = TRUE}, it also includes
  \code{std.lv}, \code{std.all}, and \code{std.nox}.
  \code{standardizedSolution()} returns a data frame with pooled
  \code{est.std}; when requested, standard errors and tests are Rubin-pooled
  from the per-imputation \pkg{lavaan} standardized-solution standard errors
  and between-imputation variability. \code{fitMeasures()} returns a named
  vector by default.
}
\details{
  \code{lavaan}'s \code{parameterEstimates()} and
  \code{standardizedSolution()} are not S3 generics. This package therefore
  exports small delegating wrappers: ordinary \code{lavaan} objects are passed
  through to \pkg{lavaan}, while \code{lavaan.survey.mi} objects receive
  Rubin-pooled estimates and Barnard-Rubin degrees of freedom.
}
\seealso{
  \code{\link{lavaan.survey.ordinal}},
  \code{\link[lavaan]{parameterEstimates}},
  \code{\link[lavaan]{standardizedSolution}},
  \code{\link[lavaan]{fitMeasures}}
}
\keyword{survey}
\keyword{models}
