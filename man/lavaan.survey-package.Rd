\name{lavaan.survey-package}
\alias{lavaan.survey-package}
\docType{package}
\title{
	Complex survey analysis of structural equation models (SEM)
}
\description{
	Takes a lavaan fit object and a complex survey design object as input
	and returns a structural equation modeling analysis based on the fit 
	object, where the complex sampling design is taken into account. 
	
	The structural equation model parameter estimates and standard errors
	are design-based. See Satorra and Muth\'en (1995) for details on the 
	procedure. This modernization fork also provides initial support for
	ordinal survey SEM through \code{\link{lavaan.survey.ordinal}}.
}
\details{
\tabular{ll}{
Package: \tab lavaan.survey\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2026-05-03\cr
License: \tab GPL(>=2)\cr
}

}
\author{
    Original author: Daniel Oberski <daniel.oberski@gmail.com>

    Maintainer of this unofficial modernization fork:
    Gregor Waller <gregor.waller@zhaw.ch>
}
\references{
  Bollen, K, Tueller, S, Oberski, DL (2013). Issues in the Structural Equation 
    Modeling of Complex Survey Data. In: Proceedings of the 59th World Statistics 
    Congress 2013 (International Statistical Institute, ed.), Hong Kong.
    \url{http://daob.nl/publications/}

  Oberski, D.L. (2014). lavaan.survey: An R Package for Complex Survey Analysis
    of Structural Equation Models. Journal of Statistical Software, 57(1), 1-27.
    \url{http://www.jstatsoft.org/v57/i01/}.

	Oberski, DL and Saris, WE (2012). A model-based procedure to evaluate
	    the relative effects of different TSE components on structural equation
	    model parameter estimates. Presentation given at the International
	    Total Survey Error Workshop in Santpoort, the Netherlands. 
	    \url{http://daob.nl/publications/}

	Satorra, A, & Bentler, PM (1994). Corrections to test statistics
		and standard errors in covariance structure analysis. 

	Satorra, A, and Muthen, BO (1995). Complex sample data in structural
	   equation modeling. Sociological methodology, 25, 267-316.
}
\keyword{ package }
\seealso{
    \code{\link{lavaan.survey}}
    \code{\link{lavaan.survey.ordinal}}
    
  \code{\link{pval.pFsum}}

  \code{\link{cardinale}}
  \code{\link{ess.dk}}
  \code{\link{ess4.gb}}
  \code{\link{liss}}
  \code{\link{pisa.be.2003}}

  \code{\link{svydesign}}
  \code{\link{svrepdesign}}

  \code{\link{lavaan}}
    \code{\link{cfa}}
    \code{\link{sem}}
}
\examples{
	# See lavaan.survey function.
}
