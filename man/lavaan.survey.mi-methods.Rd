\name{lavaan.survey.mi-methods}
\alias{lavaan.survey.mi-methods}
\alias{parameterEstimates}
\alias{parameterEstimates.lavaan.survey.mi}
\alias{fitMeasures.lavaan.survey.mi}
\title{Extractor methods for lavaan.survey multiple-imputation fits}
\description{
  Extract pooled parameter estimates and pooled fit-measure averages from
  \code{lavaan.survey.mi} objects returned by
  \code{\link{lavaan.survey.ordinal}} when Rubin parameter pooling is used.
}
\usage{
parameterEstimates(object, ...)

\method{parameterEstimates}{lavaan.survey.mi}(
  object, se = TRUE, zstat = TRUE, pvalue = TRUE, ci = TRUE,
  standardized = FALSE, level = 0.95, plabel = FALSE,
  remove.system.eq = TRUE, remove.eq = TRUE, remove.ineq = TRUE,
  remove.def = FALSE, remove.nonfree = FALSE,
  output = c("data.frame", "table"), header = FALSE, ...
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
  \item{se, zstat, pvalue, ci, level, plabel, remove.system.eq, remove.eq,
        remove.ineq, remove.def, remove.nonfree, output, header}{
    Arguments following the spirit of \code{\link[lavaan]{parameterEstimates}}.
    For \code{lavaan.survey.mi} objects, test statistics are \code{t}
    statistics with Barnard-Rubin degrees of freedom.
  }
  \item{standardized}{
    Logical. Standardized pooled estimates are not yet implemented for
    \code{lavaan.survey.mi}; if requested, unstandardized estimates are
    returned with an informative message.
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
  \code{ci.upper}. \code{fitMeasures()} returns a named vector by default.
}
\details{
  \code{lavaan}'s \code{parameterEstimates()} is not an S3 generic. This
  package therefore exports a small delegating wrapper: ordinary
  \code{lavaan} objects are passed through to \pkg{lavaan}, while
  \code{lavaan.survey.mi} objects receive Rubin-pooled estimates and
  Barnard-Rubin degrees of freedom.
}
\seealso{
  \code{\link{lavaan.survey.ordinal}},
  \code{\link[lavaan]{parameterEstimates}},
  \code{\link[lavaan]{fitMeasures}}
}
\keyword{survey}
\keyword{models}
